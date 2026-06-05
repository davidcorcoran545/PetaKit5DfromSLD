function modularPipeline(psfFolder, inputFolder, outputFolder)
% modularPipeline processes microscope data using decon+deskew and/or deskew-only.
% check whether tif files as input are processed correctly
% I may have broken them when I added an option specify the output folder

    % --- UI: Ask the User for Required Folders and Config ---
    % If paths are not provided via arguments, launch the GUI
    if nargin < 3 || isempty(psfFolder) || isempty(inputFolder) || isempty(outputFolder)
        [uiResult, isCanceled] = launchPipelineGUI();
        
        if isCanceled
            disp('Processing canceled by user.');
            return;
        end
        
        psfFolder = uiResult.psfFolder;
        config = uiResult.config;
    else
        % If run programmatically, just use defaults and the provided paths
        config = getDefaultConfig();
        config.inputFolder = inputFolder;
        config.outputFolder = outputFolder;
    end

    %% --- Suppress Bio-Formats debug/info logging ---
    try
        loci.common.DebugTools.setRootLevel('WARN');
    catch
        % If Bio-Formats is not available yet, do nothing
    end

    %% --- Build the PSF file list ---
    psfFiles = dir(fullfile(psfFolder, '*PSF_CH*.tif'));
    if isempty(psfFiles)
        error('No PSF files found in %s', psfFolder);
    end
    
    % Extract the channel number from each filename.
    psfArray = struct('channel', {}, 'fullpath', {});
    for i = 1:length(psfFiles)
        fname = psfFiles(i).name;
        tokens = regexp(fname, 'PSF_CH(\d+)', 'tokens');
        if ~isempty(tokens)
            channelNum = str2double(tokens{1}{1});
            psfArray(end+1).channel = channelNum;  %#ok<AGROW>
            psfArray(end).fullpath = fullfile(psfFiles(i).folder, fname);
        end
    end
    
    if isempty(psfArray)
        error('No valid PSF files found with pattern "PSF_CHX" in %s', psfFolder);
    end
    
    % Sort the PSF files by channel number.
    [~, idx] = sort([psfArray.channel]);
    psfArray = psfArray(idx);
    
    % Update the configuration with PSF file paths and channel patterns.
    config.PSFFullpaths = cell(1, length(psfArray));
    config.ChannelPatterns = cell(1, length(psfArray));
    for i = 1:length(psfArray)
        config.PSFFullpaths{i} = psfArray(i).fullpath;
        % Create a channel pattern (e.g. 'Ch1', 'Ch2', ...).
        config.ChannelPatterns{i} = ['Ch' num2str(psfArray(i).channel)];
    end
    
    % get size Metadata - assuming same for all PSF channels
    r=bfGetReader(psfArray(1).fullpath);    
    psf_metadata = getSizeMetadata(r, 0, config, psfArray(1).fullpath);
    r.close();
    
    fprintf('Selected PSF files:\n');
    disp(config.PSFFullpaths);
    fprintf('Channel Patterns:\n');
    disp(config.ChannelPatterns);
    fprintf('Processing files in folder: %s\n', config.inputFolder);
    fprintf('Saving output files to folder: %s\n', config.outputFolder);
    
    %% --- Process the Input Data ---
    % Determine if the input folder contains .sld, .czi or .tif files.
    sldyFiles = dir(fullfile(config.inputFolder, '*.sldy'));
    sldFiles = dir(fullfile(config.inputFolder, '*.sld'));
    allTifFiles = dir(fullfile(config.inputFolder, '*.tif'));
    allTifFiles = allTifFiles(~[allTifFiles.isdir]);
    cziFiles = dir(fullfile(config.inputFolder, '*.czi'));
    
    if ~isempty(sldFiles)
         % Process each SLD file.
         for i = 1:length(sldFiles)
             sldFullPath = fullfile(sldFiles(i).folder, sldFiles(i).name);
             processSldFile(sldFullPath, config, psf_metadata);
         end
    elseif ~isempty(cziFiles)
        %default czi config
        config = getCziDefaultConfig(config);
        %Process each CZI file
        for i=1:length(cziFiles)
            cziFullPath = fullfile(cziFiles(i).folder, cziFiles(i).name);
            processSldFile(cziFullPath, config, psf_metadata);
        end
    elseif ~isempty(sldyFiles)
        %Process each SLDY file
        for i=1:length(sldyFiles)
            sldyFullPath = fullfile(sldyFiles(i).folder, sldyFiles(i).name);
            processSldFile(sldyFullPath, config, psf_metadata);
        end
    elseif ~isempty(allTifFiles)
         % Pattern to match _T<number>_Ch<number>
         pattern = '_T\d+_Ch\d+';

         % Filter TIF files that match the pattern
         tifFiles = [];
         for i = 1:length(allTifFiles)
            filename = allTifFiles(i).name;
            if ~isempty(regexp(filename, pattern, 'once'))
                tifFiles = [tifFiles; allTifFiles(i)];  % Append matching file
            end
         end

         if ~isempty(tifFiles)
             % Only process matching TIF files
             filePaths = fullfile({tifFiles.folder}, {tifFiles.name});
             disp('Processing series of 3D TIF files...');
             %assume all tifs have same metadata
             seriesResult = processTifFolder(config);

             if isempty(seriesResult)
                 error('No valid TIFF series found in %s', config.inputFolder);
             end

             if strcmp(config.processingMode, 'decon+deskew') || strcmp(config.processingMode, 'both')
                 r=bfGetReader(filePaths{1});                 
                 tif3D_metadata = getSizeMetadata(r, 0, config, filePaths{1});
                 r.close();
                 if (abs(tif3D_metadata.pixelSizeZ - psf_metadata.pixelSizeZ)<1e-3)
                     runDeconDeskewPipeline(seriesResult, config);
                 else
                     warning('Skipping series %d Z spacing does not match PSF.', 0);
                 end
             end
             if strcmp(config.processingMode, 'deskew-only') || strcmp(config.processingMode, 'both')
                 runDeskewOnlyPipeline(seriesResult, config);
             end
             deleteIntermediateFiles(seriesResult.tifDir, config);
         else
             % No matching files, do something else
             % ignore files with ending
             ignoreSuffixes = {'_decon.tif', '_deskew.tif', '_MAX.tif','_decondeskew.tif'};
             disp('Processing Tif files...');
             for i = 1:length(allTifFiles)
                 if ~endsWith(allTifFiles(i).name,ignoreSuffixes,"IgnoreCase",true)
                     filepath = fullfile(allTifFiles(i).folder,allTifFiles(i).name);
                     processSldFile(filepath, config, psf_metadata);
                 end
             end
         end
    else
         error('No .sld or .tif files found in folder %s', config.inputFolder);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: processSldFile
function processSldFile(sldFileName, config, psf_metadata)
    % Open the .sld file using Bio-Formats.
    r = bfGetReader(sldFileName);
    if endsWith(sldFileName, {'.sld','.sldy'}, 'IgnoreCase', true)
        nSeries = r.getSeriesCount();
    else
        nSeries=1;
    end
    
    % Process each series in the .sld file.
    for S = 0:nSeries-1
        if endsWith(sldFileName, {'.sld','.sldy'}, 'IgnoreCase', true)
            r.setSeries(S);
        end
        seriesResult = convertSeriesToTif(r, S, sldFileName, config, psf_metadata);
        if isempty(seriesResult)
            continue;  % Skip series with only one Z-slice.
        end
        
        if strcmp(config.processingMode, 'decon+deskew') || strcmp(config.processingMode, 'both')            
            size_metadata = getSizeMetadata(r, S, config, sldFileName);
            if (abs(size_metadata.pixelSizeZ - psf_metadata.pixelSizeZ)<1e-3)
                runDeconDeskewPipeline(seriesResult, config);
            else
                warning('Skipping series %d Z spacing does not match PSF.', S);
            end         
        end

        if strcmp(config.processingMode, 'deskew-only') || strcmp(config.processingMode, 'both')
            runDeskewOnlyPipeline(seriesResult, config);
        end
        deleteIntermediateFiles(seriesResult.tifDir, config);
    end
    r.close();
end

%% -----------------------------------------------------------------------
%% Local Function: convertSeriesToTif
% Doesn't work well when processing files stored remotely
% The parfor loop overwhelms the network with multiple works trying to
% access the same file. 

% some speed improvements (~4-5 times faster) with parfor and two other changes. 
% It avoids copying the entire array every z-slice. 
% It may run out of memory in some situations
% Maybe bad if the size of each timepoint*number of cpu-cores is bigger than RAM. 
% In future could limit the number of workers based on the size of a timepoint and the available ram

% Also fixed a problem with original code not returning the correct time
% metadata, possibly partly due to a slidebook bug

% also check if it's doing the right thing with the number format,
% previously it was storing each plane as a double in the array, this may
% have changed to storing it as uint16. double may be necessary for later
% maths calculations. 

function seriesResult = convertSeriesToTif(r, seriesIndex, sldFileName, config, psf_metadata)    
    size_metadata = getSizeMetadata(r, seriesIndex, config, sldFileName);
    omeMeta = r.getMetadataStore();
    
    % Calculate the deskewed Z spacing
    deskewedZSpacing = sin(deg2rad(config.skewAngle)) * size_metadata.pixelSizeZ;
    
    % Get time metadata using the main reader (r) before entering parfor    
    % returns the time between the start of one timepoint and the start of the next timepoint 
    % this is not actually completely consistent between timepoints for the 3i LLSM, it is usually within +-0.8 ms of the average
    % I think this is small enough to ignore 
    frameInterval = 0;
    try
        % slidebook or the LLSM has some weird bug where the time gap between the first and second timepoint
        % is sometimes longer than the time gap between the second and third timepoint
        % so ive changed the below to measure the time between 2nd and 3rd timepoint        
        numPlanesInTimepoint = size_metadata.stackSizeZ * size_metadata.stackSizeC;
        if size_metadata.stackSizeT > 2 && numPlanesInTimepoint < omeMeta.getPlaneCount(seriesIndex)
            deltaTsecondTimepoint = omeMeta.getPlaneDeltaT(seriesIndex, numPlanesInTimepoint).value().doubleValue() / 1000;
            deltaTthirdTimepoint = omeMeta.getPlaneDeltaT(seriesIndex, numPlanesInTimepoint*2).value().doubleValue() / 1000;           
            frameInterval = deltaTthirdTimepoint - deltaTsecondTimepoint;
        end
    catch ME
        % Display a descriptive error message including the filename
        fprintf('Error retrieving time metadata, time metadata set to 0 for file: %s (Series %d)\n', sldFileName, seriesIndex);
        fprintf('Reason: %s\n', ME.message);
    end
    
    % Set up paths
    [~, baseFileName, ~] = fileparts(sldFileName);
    if endsWith(sldFileName, {'.sld','.sldy'}, 'IgnoreCase', true)
        seriesName = char(omeMeta.getImageName(seriesIndex));
        seriesNameNoSpaces = strrep(seriesName, ' ', '_');
        currentSeriesFolder = [baseFileName, '_', seriesNameNoSpaces];
    else
        currentSeriesFolder = baseFileName;
    end
    
    currentSeriesPath = fullfile(config.outputFolder, currentSeriesFolder);

    if ~exist(currentSeriesPath, 'dir')
        mkdir(currentSeriesPath);
    end
    
    tifDir = fullfile(currentSeriesPath, 'tifs');
    if ~exist(tifDir, 'dir')
        mkdir(tifDir);
    end

    % Create a directory specifically for unpadded tifs if 'both' mode is selected
    % when doing deskew without decon we don't want the tifs to be padded
    % as it breaks when doing rotate+deskew, and also max intensity projections
    % the deskew only mode never pads the tifs
    isBothMode = strcmp(config.processingMode, 'both');
    tifDirUnpadded = fullfile(currentSeriesPath, 'tifs_unpadded');
    if isBothMode && ~exist(tifDirUnpadded, 'dir')
        mkdir(tifDirUnpadded);
    end

    % 1. Extract dimensions into local variables for parallel broadcasting
    stackSizeT = size_metadata.stackSizeT;
    stackSizeC = size_metadata.stackSizeC;
    stackSizeZ = size_metadata.stackSizeZ;


    % 2. Pre-fetch one plane to get dimensions/type for workers
    sample_plane = bfGetPlane(r, 1);
    [actualRows, actualCols] = size(sample_plane);
    native_class = class(sample_plane);

    % Evaluate padding requirement outside the loop
    applyPaddingFlag = strcmp(config.processingMode, 'decon+deskew') || strcmp(config.processingMode, 'both');

    % Pre-compute Z and X/Y padding once per series.
    % Padding is applied only to deconvolution inputs, not deskew-only inputs.
    xyPadInfo = [];
    zPadInfo = [];

    if applyPaddingFlag
        zPadInfo = getSymmetricZGoodPaddingInfo(stackSizeZ, config);

        fprintf(['Z padding settings: method=%s. ', ...
            'Z %d -> %d, pad each side = %d\n'], ...
            zPadInfo.method, ...
            zPadInfo.originalZ, zPadInfo.targetZ, zPadInfo.padZ);

        xyPadInfo = getSymmetricXYGoodPaddingInfo(actualRows, actualCols, config);

        fprintf(['XY mirror padding settings: X=%d, Y=%d. ', ...
            'Y %d -> %d, pad each side = %d; ', ...
            'X %d -> %d, pad each side = %d\n'], ...
            xyPadInfo.padXEnabled, xyPadInfo.padYEnabled, ...
            xyPadInfo.originalY, xyPadInfo.targetY, xyPadInfo.padY, ...
            xyPadInfo.originalX, xyPadInfo.targetX, xyPadInfo.padX);
    end

    fprintf('Starting parallel conversion of %d timepoints...\n', stackSizeT);

    % --- BEGIN PARALLEL PROCESSING ---

    % Initialize ONE reader per worker to prevent network lock collisions 
    % and massive metadata parsing overhead.
    readerFactory = @() initWorkerReader(sldFileName, seriesIndex);
    readerCleanup = @(r) safeCloseReader(r);
    workerReaderConst = parallel.pool.Constant(readerFactory, readerCleanup);
        
    % needs testing
    % Dynamically determine safe worker count based on current hardware
    % Possibly will prevent network I/O and RAM exhaustion
    maxWorkers = determineOptimalWorkers(); 
    fprintf('Dynamically allocated %d workers for this machine.\n', maxWorkers);
    
    parfor (T = 0:stackSizeT-1, maxWorkers)
        % Suppress Bio-Formats debug/info logging
        try
            loci.common.DebugTools.setRootLevel('WARN');
        catch
            % If Bio-Formats is not available yet, do nothing
        end       
        
        % Access the persistent reader assigned to this specific worker
        worker_r = workerReaderConst.Value;
        
        for C = 0:stackSizeC-1
            % Memory Guard Logic
            % May need improving in future, doesn't do much here
            local_array = [];             
            try                
                % Try to preallocate full 3D stack (zeros)
                % This is fast and should succeed for almost all LLSM data sizes and computer memory availability 
                % If this fails catch will also probably fail later on
                local_array = zeros(actualRows, actualCols, stackSizeZ, native_class);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:nomem') || strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded')
                   warning('Could not preallocate full 3D stack. Continuing with dynamic array growth; this may still fail if memory is insufficient.');
                   local_array = [];     
                else
                    rethrow(ME);
                end
            end
            
            % Data Extraction
            for Z = 0:stackSizeZ-1
                idx = worker_r.getIndex(Z, C, T) + 1;
                plane = bfGetPlane(worker_r, idx);
                
                % Insertion
                local_array(:, :, Z + 1) = plane;
            end
            
            % Output Processing
            baseOutputArray = uint16(local_array(:, :, 1:stackSizeZ));
            
            % Filename Construction
            strS = num2str(seriesIndex);
            strT = pad(num2str(T), 4, 'left', '0');
            strC = num2str(C);
            tifFileName = fullfile(tifDir, sprintf('%s_S%s_T%s_Ch%s.tif', baseFileName, strS, strT, strC));

            % Check rotation once to save repetition
            rotateY = isfield(config,'skewDirection') && strcmpi(config.skewDirection, 'Y');

            % Write and optionally pad the tif file depending on the processing mode
            if strcmp(config.processingMode, 'deskew-only')
                if rotateY, baseOutputArray = rot90(baseOutputArray, 1); end
                parallelWriteTiff(tifFileName, baseOutputArray);

            elseif strcmp(config.processingMode, 'decon+deskew')
                paddedArray = applyZPadding(baseOutputArray, config, zPadInfo);

                % New: mirror-pad X/Y after Z padding
                paddedArray = applyXYSymmetricMirrorPadding(paddedArray, xyPadInfo);

                if rotateY
                    paddedArray = rot90(paddedArray, 1);
                end

                parallelWriteTiff(tifFileName, paddedArray);

            elseif isBothMode
                % Create both versions in memory
                paddedArray = applyZPadding(baseOutputArray, config, zPadInfo);

                % New: mirror-pad X/Y after Z padding for the decon input only
                paddedArray = applyXYSymmetricMirrorPadding(paddedArray, xyPadInfo);

                if rotateY
                    baseOutputArray = rot90(baseOutputArray, 1);
                    paddedArray = rot90(paddedArray, 1);
                end
                
                % Write padded version to main 'tifs' folder (for decon+deskew)
                parallelWriteTiff(tifFileName, paddedArray);
                
                % Write unpadded version to 'tifs_unpadded' (for deskew-only)
                tifFileNameUnpadded = fullfile(tifDirUnpadded, sprintf('%s_S%s_T%s_Ch%s.tif', baseFileName, strS, strT, strC));
                parallelWriteTiff(tifFileNameUnpadded, baseOutputArray);
            end            

        end
        
    end
    % --- END PARALLEL PROCESSING ---

    seriesResult.tifDir = tifDir;
    seriesResult.tifDirUnpadded = tifDirUnpadded; 
    seriesResult.currentSeriesFolder = currentSeriesFolder;
    seriesResult.currentSeriesPath = currentSeriesPath;
    seriesResult.frameInterval = frameInterval;
    seriesResult.pixelSizeX = size_metadata.pixelSizeX;
    seriesResult.pixelSizeY = size_metadata.pixelSizeY;
    seriesResult.deskewedZSpacing = deskewedZSpacing;
    seriesResult.xyPadInfo = xyPadInfo;
    seriesResult.zPadInfo = zPadInfo;
    
end

function size_metadata = getSizeMetadata(r, seriesIndex, config, sourceName)
    % Robust metadata reader with simple warnings and safe fallbacks

    if nargin < 3 || ~isstruct(config)
        config = getDefaultConfig();
    end
    if nargin < 4 || isempty(sourceName)
        sourceName = '<unknown source>';
    end

    % Safe fallback values
    defaultConfig = getDefaultConfig();
    if ~isfield(config, 'xyPixelSize') || ~isscalar(config.xyPixelSize) || ~isfinite(config.xyPixelSize) || config.xyPixelSize <= 0
        config.xyPixelSize = defaultConfig.xyPixelSize;
    end
    if ~isfield(config, 'dz') || ~isscalar(config.dz) || ~isfinite(config.dz) || config.dz <= 0
        config.dz = defaultConfig.dz;
    end

    omeMeta = r.getMetadataStore();

    % Dimensions (Cast to double to prevent Java object math errors downstream)
    size_metadata.stackSizeX = double(omeMeta.getPixelsSizeX(seriesIndex).getValue());
    size_metadata.stackSizeY = double(omeMeta.getPixelsSizeY(seriesIndex).getValue());
    size_metadata.stackSizeZ = double(omeMeta.getPixelsSizeZ(seriesIndex).getValue());
    size_metadata.stackSizeC = double(omeMeta.getPixelsSizeC(seriesIndex).getValue());
    size_metadata.stackSizeT = double(omeMeta.getPixelsSizeT(seriesIndex).getValue());

    % Physical sizes with fallback
    size_metadata.pixelSizeX = readPhysicalSize( ...
        @() omeMeta.getPixelsPhysicalSizeX(seriesIndex), ...
        config.xyPixelSize, 'X', sourceName, seriesIndex);

    size_metadata.pixelSizeY = readPhysicalSize( ...
        @() omeMeta.getPixelsPhysicalSizeY(seriesIndex), ...
        config.xyPixelSize, 'Y', sourceName, seriesIndex);

    size_metadata.pixelSizeZ = readPhysicalSize( ...
        @() omeMeta.getPixelsPhysicalSizeZ(seriesIndex), ...
        config.dz, 'Z', sourceName, seriesIndex);
end

function val = readPhysicalSize(getter, fallback, axisName, sourceName, seriesIndex)
    val = fallback;
    try
        obj = getter();
        if ~isempty(obj)
            tmp = double(obj.value());
            if isfinite(tmp) && tmp > 0
                val = tmp;
                return;
            end
        end
    catch
        % ignore, use fallback below
    end

    warning('modularPipeline:MetadataFallback', ...
        'Physical pixel size %s metadata unavailable or invalid for "%s" (series %d). Using fallback value %.6g um.', ...
        axisName, sourceName, seriesIndex, fallback);
end

%% -----------------------------------------------------------------------
%% Local Function: processTifFolder
% For a folder of existing 3D TIFF files, copy the usable raw TIFFs into
% the output folder and run the downstream pipeline from there.
%
% This avoids writing intermediate results into the input folder and
% prevents deleteIntermediateFiles() from deleting original input data.
function seriesResult = processTifFolder(config)
    inputTifDir = config.inputFolder;
    
    [~, currentSeriesFolder, ~] = fileparts(inputTifDir);
    if isempty(currentSeriesFolder)
        currentSeriesFolder = 'raw_tif_series';
    end
    
    currentSeriesPath = fullfile(config.outputFolder, currentSeriesFolder);
    
    if ~exist(currentSeriesPath, 'dir')
        mkdir(currentSeriesPath);
    end
    
    tifDir = fullfile(currentSeriesPath, 'tifs');
    
    if ~exist(tifDir, 'dir')
        mkdir(tifDir);
    end
    
    % Copy only likely raw input TIFFs into the output working folder.
    % Avoid copying already processed output files.
    ignoreSuffixes = {'_decon.tif', '_deskew.tif', '_MAX.tif', '_decondeskew.tif'};
    
    allTifs = dir(fullfile(inputTifDir, '*.tif'));
    allTifs = allTifs(~[allTifs.isdir]);
    
    for i = 1:length(allTifs)
        srcName = allTifs(i).name;
    
        if endsWith(srcName, ignoreSuffixes, 'IgnoreCase', true)
            continue;
        end
    
        srcPath = fullfile(allTifs(i).folder, srcName);
        dstPath = fullfile(tifDir, srcName);
    
        if ~exist(dstPath, 'file')
            copyfile(srcPath, dstPath);
        end
    end
    
    seriesResult.tifDir = tifDir;
    seriesResult.tifDirUnpadded = fullfile(currentSeriesPath, 'tifs_unpadded');
    seriesResult.currentSeriesFolder = currentSeriesFolder;
    seriesResult.currentSeriesPath = currentSeriesPath;
    seriesResult.frameInterval = 0;  % Default if metadata is unavailable.
    seriesResult.pixelSizeX = config.xyPixelSize;
    seriesResult.pixelSizeY = config.xyPixelSize;
    seriesResult.deskewedZSpacing = sin(deg2rad(config.skewAngle)) * config.dz;
    seriesResult.xyPadInfo = [];
end





%% -----------------------------------------------------------------------
%% Local Function: applyZPadding
function paddedArray = applyZPadding(array, config, zPadInfo)

if nargin < 3 || isempty(zPadInfo)
    zPadInfo = getSymmetricZGoodPaddingInfo(size(array, 3), config);
end

zPad = zPadInfo.padZ;

switch config.z_edge_padding
    case 'none'
        paddedArray = array;

    case 'zero'
        if zPad == 0
            paddedArray = array;
        else
            paddedArray = padarray(array, [0, 0, zPad], 0, 'both');
        end

    case 'mirror'
        % Mirrors the skewed data, then skews the padded areas for realism,
        % Does not duplicate the first or last real slice.
        % Needs testing more.

        [ny, nx, nz] = size(array);

        if zPad == 0
            paddedArray = array;
            return;
        end

        if nz < 2
            error('Cannot mirror-pad Z when the stack has fewer than 2 slices.');
        end

        % --- SLOPE CALCULATION ---
        % Current empirical working estimate.
        % In future, this could be calculated from dz, xyPixelSize and skewAngle.
        baseSlope = -3.55;

        % --- INITIALIZE PADDED ARRAY ---
        totalZ = nz + 2 * zPad;
        paddedArray = zeros(ny, nx, totalZ, 'like', array);

        % 1. PLACE ORIGINAL DATA IN THE CENTER
        paddedArray(:, :, zPad + 1 : zPad + nz) = array;

        % 2. BUILD Z REFLECTION INDICES WITHOUT DUPLICATING EDGE SLICES
        zIdxWithPrePad  = mirrorIndexVector(nz, zPad, 0);
        zIdxWithPostPad = mirrorIndexVector(nz, 0, zPad);

        bottomSourceZ = zIdxWithPrePad(1:zPad);
        topSourceZ    = zIdxWithPostPad(nz + 1 : nz + zPad);

        % 3. BOTTOM Z PADDING
        for i = 1:zPad
            targetZ = i;
            sourceZ_in_orig = bottomSourceZ(i);

            zDistance = targetZ - (zPad + sourceZ_in_orig);
            shiftX = zDistance * baseSlope;

            translated_slice = imtranslate( ...
                array(:, :, sourceZ_in_orig), ...
                [shiftX, 0], ...
                'Method', 'cubic', ...
                'FillValues', 0);

            colsToReplace = min(nx, ceil(abs(shiftX)) + 1);

            if colsToReplace > 0
                noise_cols = cast( ...
                    config.gaussian_mean + config.gaussian_std .* randn(ny, colsToReplace), ...
                    'like', array);

                if shiftX > 0
                    translated_slice(:, 1:colsToReplace) = noise_cols;
                elseif shiftX < 0
                    translated_slice(:, nx - colsToReplace + 1 : nx) = noise_cols;
                end
            end

            paddedArray(:, :, targetZ) = translated_slice;
        end

        % 4. TOP Z PADDING
        for i = 1:zPad
            targetZ = zPad + nz + i;
            sourceZ_in_orig = topSourceZ(i);

            zDistance = targetZ - (zPad + sourceZ_in_orig);
            shiftX = zDistance * baseSlope;

            translated_slice = imtranslate( ...
                array(:, :, sourceZ_in_orig), ...
                [shiftX, 0], ...
                'Method', 'cubic', ...
                'FillValues', 0);

            colsToReplace = min(nx, ceil(abs(shiftX)) + 1);

            if colsToReplace > 0
                noise_cols = cast( ...
                    config.gaussian_mean + config.gaussian_std .* randn(ny, colsToReplace), ...
                    'like', array);

                if shiftX > 0
                    translated_slice(:, 1:colsToReplace) = noise_cols;
                elseif shiftX < 0
                    translated_slice(:, nx - colsToReplace + 1 : nx) = noise_cols;
                end
            end

            paddedArray(:, :, targetZ) = translated_slice;
        end

    case 'gaussian'
        if zPad == 0
            paddedArray = array;
        else
            frontPad = config.gaussian_mean + config.gaussian_std .* ...
                randn(size(array, 1), size(array, 2), zPad);
            backPad = config.gaussian_mean + config.gaussian_std .* ...
                randn(size(array, 1), size(array, 2), zPad);

            frontPad = cast(frontPad, 'like', array);
            backPad  = cast(backPad,  'like', array);

            paddedArray = cat(3, frontPad, array, backPad);
        end

    case 'fixed'
        if zPad == 0
            paddedArray = array;
        else
            frontPad = cast(config.fixed_value * ...
                ones(size(array, 1), size(array, 2), zPad), 'like', array);
            backPad = cast(config.fixed_value * ...
                ones(size(array, 1), size(array, 2), zPad), 'like', array);

            paddedArray = cat(3, frontPad, array, backPad);
        end

    otherwise
        error('Invalid z_edge_padding option: %s', config.z_edge_padding);
end

end

function padInfo = getSymmetricZGoodPaddingInfo(nz, config)
% Calculates symmetric Z padding so the final Z size is good for FFT/GPU use.
%
% Requirements:
%   1. At least config.z_padding slices on each side, if Z padding is enabled
%   2. Final Z size is good for cuFFT-style FFT performance
%   3. Good sizes factor only into 2, 3, 5, and 7
%
% Assumes array layout:
%   array(Y, X, Z)
if ~isfield(config, 'z_edge_padding') || isempty(config.z_edge_padding)
    config.z_edge_padding = 'none';
end

if ~isfield(config, 'z_padding') || isempty(config.z_padding)
    config.z_padding = 0;
end

padInfo.method = config.z_edge_padding;
padInfo.originalZ = nz;
padInfo.padZ = 0;
padInfo.targetZ = nz;

% If globally disabled, return no padding.
if strcmpi(config.z_edge_padding, 'none')
    return;
end

minPad = max(0, round(config.z_padding));

% If requested Z padding is zero, leave unchanged.
if minPad == 0
    return;
end

padInfo.padZ = getOptimalSymmetricPadAmount(nz, minPad);
padInfo.targetZ = nz + 2 * padInfo.padZ;
end

%% -----------------------------------------------------------------------
%% Local Function: getSymmetricXYGoodPaddingInfo
function padInfo = getSymmetricXYGoodPaddingInfo(ny, nx, config)
% Calculates symmetric mirror-padding amounts for X and Y.
%
% Requirements:
%   1. At least config.xy_padding pixels on each enabled side
%   2. Final enabled X/Y sizes are independently good for cuFFT
%   3. Good sizes factor only into 2, 3, 5, and 7
%
% Assumes array layout:
%   array(Y, X, Z)

if ~isfield(config, 'xy_edge_padding') || isempty(config.xy_edge_padding)
    config.xy_edge_padding = 'mirror';
end

if ~isfield(config, 'xy_padding') || isempty(config.xy_padding)
    config.xy_padding = 32;
end

if ~isfield(config, 'xy_pad_x') || isempty(config.xy_pad_x)
    config.xy_pad_x = true;
end

if ~isfield(config, 'xy_pad_y') || isempty(config.xy_pad_y)
    config.xy_pad_y = true;
end

padInfo.method = config.xy_edge_padding;
padInfo.originalY = ny;
padInfo.originalX = nx;
padInfo.padXEnabled = logical(config.xy_pad_x);
padInfo.padYEnabled = logical(config.xy_pad_y);

% Default: no padding
padInfo.padY = 0;
padInfo.padX = 0;
padInfo.targetY = ny;
padInfo.targetX = nx;

% If globally disabled, return no padding
if strcmpi(config.xy_edge_padding, 'none')
    return;
end

if ~strcmpi(config.xy_edge_padding, 'mirror')
    error('Unsupported xy_edge_padding option: %s', config.xy_edge_padding);
end

% Independently enable Y padding
if padInfo.padYEnabled
    padInfo.padY = getOptimalSymmetricPadAmount(ny, config.xy_padding);
    padInfo.targetY = ny + 2 * padInfo.padY;
end

% Independently enable X padding
if padInfo.padXEnabled
    padInfo.padX = getOptimalSymmetricPadAmount(nx, config.xy_padding);
    padInfo.targetX = nx + 2 * padInfo.padX;
end
end

%% -----------------------------------------------------------------------
%% Local Function: getOptimalSymmetricPadAmount
function padAmount = getOptimalSymmetricPadAmount(origSize, minPad)
% Finds the smallest symmetric padding amount such that:
%
%   targetSize = origSize + 2 * padAmount
%
% is a good FFT size.
%
% Good FFT sizes here are those whose prime factors are only:
%   2, 3, 5, 7

allowedPrimes = [2, 3, 5, 7];

targetSize = origSize + 2 * minPad;

while true
    f = factor(targetSize);

    if isempty(f) || all(ismember(f, allowedPrimes))
        break;
    end

    % Increment by 2 so the padding remains symmetric
    targetSize = targetSize + 2;
end

padAmount = (targetSize - origSize) / 2;
end

%% -----------------------------------------------------------------------
%% Local Function: applyXYSymmetricMirrorPadding
function paddedArray = applyXYSymmetricMirrorPadding(array, padInfo)
% Applies symmetric mirror padding in Y and X.
%
% This uses reflection without duplicating the edge pixel.
%
% Assumes array layout:
%   array(Y, X, Z)

if isempty(padInfo) || strcmpi(padInfo.method, 'none')
    paddedArray = array;
    return;
end

if ~strcmpi(padInfo.method, 'mirror')
    error('Unsupported XY padding method: %s', padInfo.method);
end

idxY = mirrorIndexVector(size(array, 1), padInfo.padY, padInfo.padY);
idxX = mirrorIndexVector(size(array, 2), padInfo.padX, padInfo.padX);

paddedArray = array(idxY, idxX, :);
end

%% -----------------------------------------------------------------------
%% Local Function: removeXYSymmetricPadding
function imgOut = removeXYSymmetricPadding(imgIn, padInfo, rotateY)
% Removes symmetric X/Y padding after deconvolution.
%
% If skewDirection == 'Y', your pipeline applies:
%
%   rot90(paddedArray, 1)
%
% before writing the decon input TIFF.
%
% Therefore, after deconvolution, the padded X dimension maps to rows,
% and the padded Y dimension maps to columns.

if isempty(padInfo) || strcmpi(padInfo.method, 'none')
    imgOut = imgIn;
    return;
end

if rotateY
    % After rot90(..., 1):
    %   rows correspond to original X
    %   columns correspond to original Y
    padRows = padInfo.padX;
    padCols = padInfo.padY;
else
    % Normal orientation:
    %   rows correspond to Y
    %   columns correspond to X
    padRows = padInfo.padY;
    padCols = padInfo.padX;
end

[szY, szX, ~] = size(imgIn);

if szY <= 2 * padRows || szX <= 2 * padCols
    error('Not enough X/Y size for padding removal. Image size is [%d, %d], padRows=%d, padCols=%d.', ...
        szY, szX, padRows, padCols);
end

imgOut = imgIn( ...
    (padRows + 1):(szY - padRows), ...
    (padCols + 1):(szX - padCols), ...
    :);
end

%% -----------------------------------------------------------------------
%% Local Function: mirrorIndexVector
function idx = mirrorIndexVector(n, padPre, padPost)
% Creates reflected indices for mirror padding without duplicating edge pixels.
%
% Example:
%   Original:
%       [1 2 3 4]
%
%   padPre = 2, padPost = 3 gives index vector:
%       [3 2 1 2 3 4 3 2 1]
%
% This is reflection padding rather than edge replication.

if padPre == 0 && padPost == 0
    idx = 1:n;
    return;
end

if n < 2
    error('Cannot mirror-pad a dimension of size less than 2.');
end

positions = (1 - padPre):(n + padPost);

period = 2 * n - 2;

idx = mod(positions - 1, period) + 1;

over = idx > n;
idx(over) = period - idx(over) + 2;
end

%% -----------------------------------------------------------------------
%% Local Function: runDeconDeskewPipeline
function runDeconDeskewPipeline(seriesResult, config)
    fprintf('Running deconvolution+deskew pipeline for series: %s\n', seriesResult.currentSeriesFolder);
    
    % --- Deconvolution step ---
    if isfield(config, 'deconAlgorithm') && strcmp(config.deconAlgorithm, 'RLGC')
        % RLGC Branch
        deconDir = fullfile(seriesResult.tifDir, config.resultDirName);
        if ~exist(deconDir, 'dir')
            mkdir(deconDir);
        end
        
        fileList = dir(fullfile(seriesResult.tifDir, '*.tif'));
        for fIdx = 1:length(fileList)
            fileName = fileList(fIdx).name;
            filePath = fullfile(fileList(fIdx).folder, fileName);
            
            % Extract channel number from filename to match correct PSF
            tokens = regexp(fileName, '_Ch(\d+)', 'tokens');
            if ~isempty(tokens)
                chStr = ['Ch' tokens{1}{1}];
                psfIdx = find(strcmp(config.ChannelPatterns, chStr), 1);
                
                if ~isempty(psfIdx)
                    psfPath = config.PSFFullpaths{psfIdx};
                    outPath = fullfile(deconDir, fileName);
                    RLGC(filePath, psfPath, outPath, config.Background, config.Save16bit);
                else
                    warning('Could not find matching PSF for channel %s. Skipping %s.', chStr, fileName);
                end
            end
        end
    else
        % Original Branch
        % Test whether this could be faster, doesn't seem to use many cores
        XR_decon_data_wrapper(seriesResult.tifDir, 'resultDirName', config.resultDirName, 'xyPixelSize', config.xyPixelSize, ...
                    'dz', config.dz, 'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'PSFFullpaths', config.PSFFullpaths, ...
                    'dzPSF', config.dzPSF, 'parseSettingFile', config.parseSettingFile, 'RLmethod', config.RLmethod, ...
                    'wienerAlpha', config.wienerAlpha, 'OTFCumThresh', config.OTFCumThresh, 'skewed', config.skewed, ...
                    'Background', config.Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', config.DeconIter, ...
                    'fixIter', config.fixIter, 'EdgeErosion', config.EdgeErosion, 'Save16bit', config.Save16bit, ...
                    'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'parseCluster', config.parseCluster, ...
                    'largeFile', config.largeFile, 'GPUJob', config.GPUJob, 'debug', config.debug, 'cpusPerTask', config.cpusPerTask, ...
                    'ConfigFile', config.ConfigFile, 'GPUConfigFile', config.GPUConfigFile, 'mccMode', config.mccMode);
    end
    
    if config.GPUJob && gpuDeviceCount('available') > 0
         reset(gpuDevice);
    end


    % Remove X/Y and Z padding from the decon results, before we deskew.
    deconDir = fullfile(seriesResult.tifDir, config.resultDirName);
    removePaddingFromDir(deconDir, config, seriesResult);


    % --- Deskew step ---.
    dataPath_exps = fullfile(seriesResult.tifDir, config.resultDirName);
    XR_deskew_rotate_data_wrapper(dataPath_exps, 'skewAngle', config.skewAngle, 'flipZstack', config.flipZstack, ...
        'DSRCombined', config.DSRCombined, 'rotate', config.rotate, 'xyPixelSize', config.xyPixelSize, 'dz', config.dz, ...
        'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'largeFile', config.largeFile, ...
        'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'Save16bit', config.Save16bit, 'parseCluster', config.parseCluster, ...
        'masterCompute', config.masterCompute, 'configFile', config.configFile, 'mccMode', config.mccMode);
    
    % Merge the deconvolved+deskewed images.
    outputTiffFile = fullfile(config.outputFolder, [seriesResult.currentSeriesFolder, '_decondeskew.tif']);
    
    % If rotation is enabled then merge the 'DSR' folder, otherwise 'DS'    
    if config.rotate
        deconDSDir = fullfile(dataPath_exps, 'DSR');
    else
        deconDSDir = fullfile(dataPath_exps, 'DS');
    end
    
    if isfield(config,'skewDirection')
        paraMergeTiffFilesToMultiDimStack(deconDSDir, outputTiffFile, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeTiffFilesToMultiDimStack(deconDSDir, outputTiffFile, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval);
    end

    outputTiffFileMax = fullfile(config.outputFolder, [seriesResult.currentSeriesFolder, '_decondeskew_MAX.tif']);
    inputToMergeMax = fullfile(deconDSDir, 'MIPs');
    if isfield(config,'skewDirection')
        paraMergeMaxToStack(inputToMergeMax, outputTiffFileMax, seriesResult.pixelSizeX, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeMaxToStack(inputToMergeMax, outputTiffFileMax, seriesResult.pixelSizeX, seriesResult.frameInterval);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: runDeskewOnlyPipeline
function runDeskewOnlyPipeline(seriesResult, config)
    fprintf('Running deskew-only pipeline for series: %s\n', seriesResult.currentSeriesFolder);
    
    % Use the unpadded directory if we used the 'both' deskew and deskew+decon mode, otherwise use the standard directory
    if strcmp(config.processingMode, 'both')
        inputDir = seriesResult.tifDirUnpadded;
    else
        inputDir = seriesResult.tifDir;
    end

    XR_deskew_rotate_data_wrapper(inputDir, 'resultDirName', config.resultDirNameDeskew, 'skewAngle', config.skewAngle, 'flipZstack', config.flipZstack, ...
        'DSRCombined', config.DSRCombined, 'rotate', config.rotate, 'xyPixelSize', config.xyPixelSize, 'dz', config.dz, ...
        'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'largeFile', config.largeFile, ...
        'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'Save16bit', config.Save16bit, 'parseCluster', config.parseCluster, ...
        'masterCompute', config.masterCompute, 'configFile', config.configFile, 'mccMode', config.mccMode);
    
    % Merge the deskew-only images.
    outputTiffFileDeskew = fullfile(config.outputFolder, [seriesResult.currentSeriesFolder, '_deskew.tif']);
    
    % If rotation is enabled then merge the 'DSR' folder, otherwise 'DS'
    if config.rotate
        deskewDSDir = fullfile(inputDir, [config.resultDirNameDeskew, 'R']); % e.g., 'DSR'
    else
        deskewDSDir = fullfile(inputDir, config.resultDirNameDeskew); % e.g., 'DS'
    end
    
    if isfield(config,'skewDirection')
        paraMergeTiffFilesToMultiDimStack(deskewDSDir, outputTiffFileDeskew, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeTiffFilesToMultiDimStack(deskewDSDir, outputTiffFileDeskew, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval);
    end
    
    outputTiffFileDeskewMax = fullfile(config.outputFolder, [seriesResult.currentSeriesFolder, '_deskew_MAX.tif']);
    inputToMergeDeskewMax = fullfile(deskewDSDir, 'MIPs');
    if isfield(config,'skewDirection')
        paraMergeMaxToStack(inputToMergeDeskewMax, outputTiffFileDeskewMax, seriesResult.pixelSizeX, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeMaxToStack(inputToMergeDeskewMax, outputTiffFileDeskewMax, seriesResult.pixelSizeX, seriesResult.frameInterval);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: removePaddingFromDir
function removePaddingFromDir(targetDir, config, seriesResult)

    fileList = dir(fullfile(targetDir, '*.tif'));

    hasXYPadInfo = nargin >= 3 && ...
                   isstruct(seriesResult) && ...
                   isfield(seriesResult, 'xyPadInfo') && ...
                   ~isempty(seriesResult.xyPadInfo);

    rotateY = isfield(config, 'skewDirection') && strcmpi(config.skewDirection, 'Y');

    for i = 1:length(fileList)
        filePath = fullfile(targetDir, fileList(i).name);
        img = parallelReadTiff(filePath);

        % Remove X/Y padding first
        if hasXYPadInfo
            img = removeXYSymmetricPadding(img, seriesResult.xyPadInfo, rotateY);
        end

        % Remove Z padding second.
        % Use the actual calculated GPU-friendly Z pad, not just config.z_padding.
        zPadToRemove = 0;

        if nargin >= 3 && ...
                isstruct(seriesResult) && ...
                isfield(seriesResult, 'zPadInfo') && ...
                ~isempty(seriesResult.zPadInfo) && ...
                isfield(seriesResult.zPadInfo, 'padZ')

            zPadToRemove = seriesResult.zPadInfo.padZ;

        elseif isfield(config, 'z_edge_padding') && ...
                ~strcmp(config.z_edge_padding, 'none') && ...
                isfield(config, 'z_padding')

            % Backwards-compatible fallback.
            zPadToRemove = config.z_padding;
        end

        if zPadToRemove > 0
            if size(img, 3) > 2 * zPadToRemove
                img = img(:, :, (zPadToRemove + 1):(end - zPadToRemove));
            else
                warning('Not enough Z depth for padding removal in file: %s', fileList(i).name);
            end
        end

        parallelWriteTiff(filePath, img);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: deleteIntermediateFiles
function deleteIntermediateFiles(tifDir, config)
    if config.deleteRawTif
         deleteFilesInDir(tifDir);
         % Also clean up unpadded tifs if 'both' mode was used
         deleteFilesInDir(fullfile(tifDir, '..', 'tifs_unpadded')); 
    end
    
    if config.deleteDeconTif
         % Base Deconvolution Directory
         deconDir = fullfile(tifDir, config.resultDirName);
         
         % Define all possible intermediate folders to ensure a clean sweep.
         % This catches both standard (DS) and rotated (DSR) outputs.
         foldersToClean = {
             deconDir, ...
             fullfile(deconDir, 'MIPs'), ...
             fullfile(deconDir, 'DS'), ...
             fullfile(deconDir, 'DS', 'MIPs'), ...
             fullfile(deconDir, 'DSR'), ...
             fullfile(deconDir, 'DSR', 'MIPs'), ...
             fullfile(tifDir, config.resultDirNameDeskew), ...
             fullfile(tifDir, config.resultDirNameDeskew, 'MIPs'), ...
             fullfile(tifDir, [config.resultDirNameDeskew, 'R']), ... % e.g., 'DSR'
             fullfile(tifDir, [config.resultDirNameDeskew, 'R'], 'MIPs')
         };
         
         % Loop through and delete .tif files safely
         for i = 1:length(foldersToClean)
             deleteFilesInDir(foldersToClean{i});
         end
    end
end

%% Helper Function: deleteFilesInDir
% Checks if the directory exists and deletes all .tif files inside to prevent repetitive code.
function deleteFilesInDir(targetDir)
    if exist(targetDir, 'dir')
        files = dir(fullfile(targetDir, '*.tif'));
        for k = 1:length(files)
            delete(fullfile(files(k).folder, files(k).name));
        end
    end
end

%% -----------------------------------------------------------------------
%% Local Function: initWorkerReader
function r = initWorkerReader(fileName, seriesIndex)
try
    % 1. Force the MATLAB worker to load the Bio-Formats Java library
    % The '1' argument prevents it from prompting the user in the console
    bfCheckJavaPath(1);

    % 2. Initialize the reader without the Memoizer
    baseReader = loci.formats.ImageReader();
    filler = loci.formats.ChannelFiller(baseReader);
    r = loci.formats.ChannelSeparator(filler);
    r.setId(fileName);
    r.setSeries(seriesIndex);
catch ME
    error('Failed to initialise Bio-Formats reader on worker for file "%s", series %d. Reason: %s', ...
        fileName, seriesIndex, ME.message);
end
end

%% -----------------------------------------------------------------------
%% Local Function: safeCloseReader
function safeCloseReader(r)
% Defensively close the reader to prevent "Dot indexing" crashes 
% if the reader failed to initialize properly.
if ~isempty(r) && (isobject(r) || isjava(r))
    try
        r.close();
    catch
        % Ignore errors during cleanup
    end
end
end

% move this in future
function config = getCziDefaultConfig(config)
    config.xyPixelSize = 0.1449922;
    config.skewAngle = 30.0;
    config.skewDirection = 'Y';
end

%% -----------------------------------------------------------------------
%% Local Function: determineOptimalWorkers
function optimalWorkers = determineOptimalWorkers()
% Dynamically calculates safe worker limits based on the computer's hardware

% 1. Get the number of physical CPU cores (not logical threads)
physicalCores = feature('numcores');

% 2. Calculate a safe baseline (Use half of available physical cores)
% A 4-core laptop gets 2 workers. Your 26-core machine gets 13 workers.
baseWorkers = floor(physicalCores / 2);

% 3. Set a strict ceiling to protect the hard drive from I/O bottlenecking
% Even on a massive supercomputer, 12 workers reading TIFFs simultaneously
% is usually the absolute limit for standard NVMe SSDs.
ioCeiling = 12; 

% 4. Choose the safest number
optimalWorkers = min(baseWorkers, ioCeiling);

% Ensure it never returns less than 1
optimalWorkers = max(1, optimalWorkers);
end