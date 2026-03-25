function paraMergeTiffFilesToMultiDimStack(inputFolder, outputFilePath, xySpacing, zSpacing, frameInterval, skewDirection)
    tic;
    if ~ischar(outputFilePath)
        outputFilePath = char(outputFilePath);
    end

    fprintf('Listing all TIFF files...\n');
    tiffFiles = dir(fullfile(inputFolder, '*.tif'));
    if isempty(tiffFiles)
        error('No TIFF files found.');
    end

    % --- Metadata parsing ---
    fileNames = {tiffFiles.name}';
    tokens = regexp(fileNames, '.*_T(\d+)_Ch(\d+).tif', 'tokens', 'once');
    metaMatrix = cellfun(@(x) str2double(x), vertcat(tokens{:}));
    [~, sortIdx] = sortrows(metaMatrix);
    sortedFileNames = fileNames(sortIdx);
    sortedMeta = metaMatrix(sortIdx, :);

    numTimePoints = max(sortedMeta(:,1)) + 1;
    numChannels   = max(sortedMeta(:,2)) + 1;

    % --- Pre-read first file for dimensions ---
    doRotate = (nargin > 5) && strcmpi(skewDirection, 'Y');
    firstImg = readTiffSerial(fullfile(inputFolder, sortedFileNames{1}));
    if doRotate
        firstImg = rot90(firstImg, -1);
    end
    [stackSizeY, stackSizeX, stackSizeZ] = size(firstImg);
    clear firstImg;
    fprintf('X:%d, Y:%d, Z:%d, Ch:%d, T:%d\n', stackSizeX, stackSizeY, stackSizeZ, numChannels, numTimePoints);

    % --- Ensure parallel pool is warmed up ---
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local');
    end
    numWorkers = pool.NumWorkers;
    fprintf('Using %d parallel workers.\n', numWorkers);

    % --- TIFF output setup ---
    t = Tiff(outputFilePath, 'w8');

    tagstruct.ImageLength        = stackSizeY;
    tagstruct.ImageWidth         = stackSizeX;
    tagstruct.Photometric        = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample      = 16;
    tagstruct.SamplesPerPixel    = 1;
    tagstruct.RowsPerStrip       = stackSizeY;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software           = 'MATLAB';
    tagstruct.ResolutionUnit     = Tiff.ResolutionUnit.Centimeter;
    tagstruct.XResolution        = 10000 / xySpacing;
    tagstruct.YResolution        = 10000 / xySpacing;
    tagstruct.ImageDescription   = sprintf(['ImageJ=1.52p\nimages=%d\nchannels=%d\nslices=%d\nframes=%d\n' ...
        'hyperstack=true\nmode=grayscale\nloop=false\nspacing=%f\n' ...
        'unit=micron\nfinterval=%f\n'], ...
        stackSizeZ * numChannels * numTimePoints, numChannels, stackSizeZ, numTimePoints, zSpacing, frameInterval);

    % Verify RowsPerStrip is honoured
    t.setTag(tagstruct);
    actualRPS = t.getTag('RowsPerStrip');
    if actualRPS ~= stackSizeY
        fprintf('WARNING: RowsPerStrip requested %d but got %d. Writing will use more strips per slice.\n', stackSizeY, actualRPS);
    else
        fprintf('RowsPerStrip: %d (full frame, good).\n', actualRPS);
    end

    % --- Determine batch size ---
    % With 1-2 channels, batch multiple timepoints so all workers stay busy.
    % e.g. 16 workers / 2 channels = 8 timepoints per batch
    filesPerTimepoint = numChannels;
    batchSize = max(1, floor(numWorkers / filesPerTimepoint));
    fprintf('Batching %d timepoint(s) per parfor call (%d files/batch).\n', batchSize, batchSize * filesPerTimepoint);

    % --- Pre-build full file path list for parfor (avoids broadcast of struct) ---
    allFilePaths = fullfile(inputFolder, sortedFileNames);  % Nx1 cell of full paths

    % --- Main loop: read in T*C batches, write sequentially ---
    totalSlices = stackSizeZ * numChannels * numTimePoints;
    slicesWritten = 0;

    tIdx = 1;
    while tIdx <= numTimePoints
        batchEnd    = min(tIdx + batchSize - 1, numTimePoints);
        batchT      = tIdx:batchEnd;
        numInBatch  = numel(batchT);
        totalFilesInBatch = numInBatch * numChannels;

        % Build flat list of file paths for this batch
        batchPaths = cell(1, totalFilesInBatch);
        for bIdx = 1:numInBatch
            for cIdx = 1:numChannels
                flatIdx = (bIdx-1)*numChannels + cIdx;
                
                % Calculate the global timepoint index (0-based)
                % batchT(bIdx) is the current timepoint number (1 to numTimePoints)
                currentTime = batchT(bIdx) - 1; 
                currentChan = cIdx - 1;
                
                % Map to the 1-based index of allFilePaths
                fileListIdx = (currentTime * numChannels) + (currentChan + 1);
                
                batchPaths{flatIdx} = allFilePaths{fileListIdx};
            end
        end

        % --- Parallel read: one worker per file ---
        batchImgs = cell(1, totalFilesInBatch);
        parfor f = 1:totalFilesInBatch
            img = readTiffSerial(batchPaths{f});
            if doRotate
                img = rot90(img, -1);
            end
            batchImgs{f} = img;
        end

        % --- Sequential write: unpack batch and write in ImageJ order ---
        for bIdx = 1:numInBatch
            for zIdx = 1:stackSizeZ
                for cIdx = 1:numChannels
                    flatIdx = (bIdx-1)*numChannels + cIdx;
                    t.write(batchImgs{flatIdx}(:, :, zIdx));
                    t.writeDirectory();
                    slicesWritten = slicesWritten + 1;
                    % setTag required after writeDirectory, except after the very last slice
                    if slicesWritten < totalSlices
                        t.setTag(tagstruct);
                    end
                end
            end
            fprintf('Processed Timepoint %d/%d\n', batchT(bIdx), numTimePoints);
        end

        % Clear batch from memory before next iteration
        clear batchImgs;
        tIdx = tIdx + batchSize;
    end

    t.close();
    fprintf('Success! Total time: %.2f seconds.\n', toc);
end


% -------------------------------------------------------------------------
% Serial TIFF reader — used inside parfor (no nested parallelism)
% -------------------------------------------------------------------------
function img = readTiffSerial(filePath)
    tObj = Tiff(filePath, 'r');

    firstSlice = tObj.read();
    [h, w] = size(firstSlice);

    % Count Z slices by walking directories
    numSlices = 1;
    while ~tObj.lastDirectory()
        tObj.nextDirectory();
        numSlices = numSlices + 1;
    end

    tObj.setDirectory(1);
    img = zeros(h, w, numSlices, 'uint16');
    img(:, :, 1) = firstSlice;
    for z = 2:numSlices
        tObj.nextDirectory();
        img(:, :, z) = tObj.read();
    end
    tObj.close();
end