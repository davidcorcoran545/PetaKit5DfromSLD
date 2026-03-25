% processes 100 timepoints in 36 seconds compared to 88 seconds for the original.
% this does the reading and writing at the same time 
% but perhaps is more likely to break or have memory issues 
% so probably best not to use this for now

function paraMergeTiffFilesToMultiDimStackVibeClaude4(inputFolder, outputFilePath, xySpacing, zSpacing, frameInterval, skewDirection)
    tic;
    if ~ischar(outputFilePath), outputFilePath = char(outputFilePath); end

    fprintf('Listing all TIFF files...\n');
    tiffFiles = dir(fullfile(inputFolder, '*.tif'));
    if isempty(tiffFiles), error('No TIFF files found.'); end

    % --- Metadata parsing ---
    fileNames = {tiffFiles.name}';
    tokens = regexp(fileNames, '.*_T(\d+)_Ch(\d+).tif', 'tokens', 'once');
    metaMatrix = cellfun(@(x) str2double(x), vertcat(tokens{:}));
    [~, sortIdx] = sortrows(metaMatrix);
    sortedFileNames = fileNames(sortIdx);
    
    numTimePoints = max(metaMatrix(:,1)) + 1;
    numChannels   = max(metaMatrix(:,2)) + 1;

    % --- Dimensions & Rotation ---
    doRotate = (nargin > 5) && strcmpi(skewDirection, 'Y');
    firstImg = readTiffSerial(fullfile(inputFolder, sortedFileNames{1}));
    if doRotate, firstImg = rot90(firstImg, -1); end
    [stackSizeY, stackSizeX, stackSizeZ] = size(firstImg);
    clear firstImg;

    % --- Parallel Pool ---
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local'); end
    numWorkers = pool.NumWorkers;

    % --- TIFF setup ---
    t = Tiff(outputFilePath, 'w8');
    tagstruct.ImageLength = stackSizeY;
    tagstruct.ImageWidth = stackSizeX;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = stackSizeY;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
    tagstruct.XResolution = 10000 / xySpacing;
    tagstruct.YResolution = 10000 / xySpacing;
    tagstruct.ImageDescription = sprintf('ImageJ=1.52p\nimages=%d\nchannels=%d\nslices=%d\nframes=%d\nhyperstack=true\nmode=grayscale\nloop=false\nspacing=%f\nunit=micron\nfinterval=%f\n', ...
        stackSizeZ * numChannels * numTimePoints, numChannels, stackSizeZ, numTimePoints, zSpacing, frameInterval);
    
    t.setTag(tagstruct);

    % --- Pipelining Setup ---
    % We process 1 timepoint per future for maximum granularity
    allFilePaths = fullfile(inputFolder, sortedFileNames);
    totalSlices = stackSizeZ * numChannels * numTimePoints;
    slicesWritten = 0;
    
    % Queue depth: How many timepoints to keep in RAM at once. 
    % Adjust based on RAM (e.g., 2-4 is usually enough to hide latency).
    maxQueueDepth = min(numWorkers, 4);     
    % % Automatically calculate queue depth based on available memory
    % % 1. Calculate size of one timepoint in GB
    % % (Height * Width * Slices * Channels * 2 bytes) / 1024^3
    % bytesPerTimepoint = (stackSizeY * stackSizeX * stackSizeZ * numChannels * 2);
    % timepointGB = bytesPerTimepoint / (1024^3);
    % 
    % % 2. Calculate how many fit in 50% of available RAM
    % mem = memory;
    % availableGB = mem.MaxPossibleArrayBytes / (1024^3);
    % ramLimit = floor((availableGB * 0.5) / timepointGB);
    % 
    % % 3. Set a sensible performance cap
    % % We want enough to hide latency, but not so many that MATLAB stutters
    % perfCap = max(numWorkers * 2, 20); 
    % 
    % % 4. Final Queue Depth
    % maxQueueDepth = min(ramLimit, perfCap);
    % 
    % fprintf('Each timepoint is %.2f MB. Using Queue Depth of %d.\n', ...
    %         bytesPerTimepoint/1024^2, maxQueueDepth);
    
    futures = parallel.FevalFuture.empty(0, numTimePoints);

    fprintf('Starting asynchronous read/write pipeline...\n');

    % Initial queueing
    currT = 1;
    while currT <= min(maxQueueDepth, numTimePoints)
        % Paths for all channels of this timepoint
        indices = (currT-1)*numChannels + (1:numChannels);
        timepointPaths = allFilePaths(indices);
        
        futures(currT) = parfeval(pool, @readTimepointBatch, 1, timepointPaths, doRotate);
        currT = currT + 1;
    end

    % --- Main Loop: Fetch and Replenish ---
    for tIdx = 1:numTimePoints
        % 1. Wait for the specific future in order (maintains T-sequence)
        batchImgs = fetchOutputs(futures(tIdx));
        
        % 2. Write this timepoint to disk
        for zIdx = 1:stackSizeZ
            for cIdx = 1:numChannels
                t.write(batchImgs{cIdx}(:, :, zIdx));
                slicesWritten = slicesWritten + 1;
                if slicesWritten < totalSlices
                    t.writeDirectory();
                    t.setTag(tagstruct);
                end
            end
        end
        fprintf('Written Timepoint %d/%d\n', tIdx, numTimePoints);
        
        % 3. Clear from memory
        batchImgs = []; 
        
        % 4. Queue the next available timepoint to keep the pipeline full
        if currT <= numTimePoints
            indices = (currT-1)*numChannels + (1:numChannels);
            timepointPaths = allFilePaths(indices);
            futures(currT) = parfeval(pool, @readTimepointBatch, 1, timepointPaths, doRotate);
            currT = currT + 1;
        end
    end

    t.close();
    fprintf('Success! Total time: %.2f seconds.\n', toc);
end

% Helper to read all channels for a single timepoint inside a worker
function imgs = readTimepointBatch(paths, doRotate)
    numCh = numel(paths);
    imgs = cell(1, numCh);
    for c = 1:numCh
        img = readTiffSerial(paths{c});
        if doRotate, img = rot90(img, -1); end
        imgs{c} = img;
    end
end

function img = readTiffSerial(filePath)
    tObj = Tiff(filePath, 'r');
    firstSlice = tObj.read();
    [h, w] = size(firstSlice);
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