function modularPipeline(psfFolder, inputFolder)
% modularPipeline prompts the user for:
%    1. A folder containing PSF files. Assumes PSF filenames contain the text
%       'PSF_CHX' (where X is a number) and uses those to build the PSF paths.
%    2. A folder containing the files to deconvolve (either .sld or .tif).
%
% The code then processes the data using one or both pipelines (decon+deskew and/or
% deskew-only), with z–axis padding applied during conversion and removed before merging.

    %% --- UI: Ask the User for Required Folders ---
    % Allow test code to pass paths directly
    if nargin < 1 || isempty(psfFolder)
        psfFolder = uigetdir([], 'Select the folder containing your PSF files');
        if psfFolder == 0
            error('No PSF folder selected.');
        end
    end

    if nargin < 2 || isempty(inputFolder)
        inputFolder = uigetdir([], 'Select the folder containing files to deconvolve (SLD or TIFF)');
        if inputFolder == 0
            error('No input folder selected.');
        end
    end
    
    %% --- Get Default Configuration and Update with User Choices ---
    config = getDefaultConfig();
    config.inputFolder = inputFolder;
    
    % Build the PSF file list by scanning the PSF folder.
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
    psf_metadata = getSizeMetadata(r,0);
    r.close();
    fprintf('Selected PSF files:\n');
    disp(config.PSFFullpaths);
    fprintf('Channel Patterns:\n');
    disp(config.ChannelPatterns);
    fprintf('Processing files in folder: %s\n', config.inputFolder);
    
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
                    tif3D_metadata = getSizeMetadata(r,0);
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
            size_metadata = getSizeMetadata(r,S);    
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
    size_metadata = getSizeMetadata(r, seriesIndex);
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
    
    currentSeriesPath = fullfile(config.inputFolder, currentSeriesFolder);
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

    fprintf('Starting parallel conversion of %d timepoints...\n', stackSizeT);

    % --- BEGIN PARALLEL PROCESSING ---
    parfor T = 0:stackSizeT-1
        % Create a worker-specific reader for this timepoint
        worker_r = bfGetReader(sldFileName);
        worker_r.setSeries(seriesIndex);
        
        for C = 0:stackSizeC-1
            % Memory Guard Logic
            use_fast_mode = true;
            local_array = []; 
            
            try
                % Attempt fast preallocation (zeros)
                local_array = zeros(actualRows, actualCols, stackSizeZ, native_class);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:nomem') || strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded')
                    use_fast_mode = false; % Fallback to safe growth
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
            rotateY = isfield(config,'skewDirection') && (config.skewDirection == 'Y');

            % Write and optionally pad the tif file depending on the processing mode 
            if strcmp(config.processingMode, 'deskew-only')
                if rotateY, baseOutputArray = rot90(baseOutputArray, 1); end
                parallelWriteTiff(tifFileName, baseOutputArray);
                
            elseif strcmp(config.processingMode, 'decon+deskew')
                paddedArray = applyZPadding(baseOutputArray, config);
                if rotateY, paddedArray = rot90(paddedArray, 1); end
                parallelWriteTiff(tifFileName, paddedArray);
                
            elseif isBothMode
                % Create both versions in memory
                paddedArray = applyZPadding(baseOutputArray, config);
                
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
        
        % Close the worker reader to release the file lock
        worker_r.close();
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
    
end

function size_metadata = getSizeMetadata(r, seriesIndex)
    omeMeta = r.getMetadataStore();
    % Extract image dimensions.
    size_metadata.stackSizeX = omeMeta.getPixelsSizeX(seriesIndex).getValue();
    size_metadata.stackSizeY = omeMeta.getPixelsSizeY(seriesIndex).getValue();
    size_metadata.stackSizeZ = omeMeta.getPixelsSizeZ(seriesIndex).getValue();
    size_metadata.stackSizeC = omeMeta.getPixelsSizeC(seriesIndex).getValue();
    size_metadata.stackSizeT = omeMeta.getPixelsSizeT(seriesIndex).getValue();
    
    % Get physical pixel sizes. If unavailable (or NaN) then substitute defaults from config.
    pixelSizeX_obj = omeMeta.getPixelsPhysicalSizeX(seriesIndex);
    if isempty(pixelSizeX_obj)
        size_metadata.pixelSizeX = config.xyPixelSize;
    else
        size_metadata.pixelSizeX = double(pixelSizeX_obj.value());
        if isnan(size_metadata.pixelSizeX)
            size_metadata.pixelSizeX = config.xyPixelSize;
        end
    end
    
    pixelSizeY_obj = omeMeta.getPixelsPhysicalSizeY(seriesIndex);
    if isempty(pixelSizeY_obj)
        size_metadata.pixelSizeY = config.xyPixelSize;
    else
        size_metadata.pixelSizeY = double(pixelSizeY_obj.value());
        if isnan(size_metadata.pixelSizeY)
            size_metadata.pixelSizeY = config.xyPixelSize;
        end
    end
    
    pixelSizeZ_obj = omeMeta.getPixelsPhysicalSizeZ(seriesIndex);
    if isempty(pixelSizeZ_obj)
        size_metadata.pixelSizeZ = config.dz;
    else
        size_metadata.pixelSizeZ = double(pixelSizeZ_obj.value());
        if isnan(size_metadata.pixelSizeZ)
            size_metadata.pixelSizeZ = config.dz;
        end
    end
end

%% -----------------------------------------------------------------------
%% Local Function: processTifFolder
function seriesResult = processTifFolder(config)
    % For a folder of TIFF files, assume that the folder itself contains the raw images.
    tifDir = config.inputFolder;
    [~, currentSeriesFolder, ~] = fileparts(tifDir);
    if isempty(currentSeriesFolder)
       currentSeriesFolder = 'raw_tif_series';
    end
    
    seriesResult.tifDir = tifDir;
    seriesResult.currentSeriesFolder = currentSeriesFolder;
    seriesResult.currentSeriesPath = tifDir;
    seriesResult.frameInterval = 0;  % Default if metadata is unavailable.
    seriesResult.pixelSizeX = config.xyPixelSize;
    seriesResult.deskewedZSpacing = sin(deg2rad(config.skewAngle)) * config.dz;
end

%% -----------------------------------------------------------------------
%% Local Function: applyZPadding
function paddedArray = applyZPadding(array, config)
switch config.z_edge_padding
    case 'none'
        paddedArray = array;
    case 'zero'
        paddedArray = padarray(array, [0, 0, config.z_padding], 0, 'both');
    case 'mirror'
        % Mirrors the skewed data, skews the padded areas for realism, 
        % --- SETUP ---
        [ny, nx, nz] = size(array);
        zPad = config.z_padding;
        
        % --- SLOPE CALCULATION ---
        % z_step (0.5) / x_pixel (0.104) * tan(32.8)
        % Theoretical slope is ~-3.0978 pixels. Using -3.55 as a working estimate.
        % In future replace with calculated version depending on z-step xpixel and angle, also allow zeiss/czi
        baseSlope = -3.55;         
        
        % --- INITIALIZE PADDED ARRAY ---
        totalZ = nz + (2 * zPad);
        paddedArray = zeros(ny, nx, totalZ, 'like', array);
        
        % 1. PLACE ORIGINAL DATA IN THE CENTER
        paddedArray(:, :, zPad+1 : zPad+nz) = array;
        
        % 2. BOTTOM PADDING (Slices 1 to zPad)
        for i = 1:zPad
            targetZ = i;
            sourceZ_in_orig = zPad - i + 1;
            zDistance = targetZ - (zPad + sourceZ_in_orig);
            shiftX = zDistance * baseSlope;
            
            % Apply geometric transformation
            translated_slice = imtranslate(array(:,:,sourceZ_in_orig), [shiftX, 0], 'Method', 'cubic', 'FillValues', 0);
            
            % The previous transformation creates empty regions filled with zeros, now we will fill these with a gaussian noise to better simulate camera noise.
            % Calculate affected columns + 1 extra column
            colsToReplace = ceil(abs(shiftX)) + 1;
            
            if colsToReplace > 0
                % Generate Gaussian noise
                noise_cols = cast(config.gaussian_mean + config.gaussian_std .* randn(ny, colsToReplace), 'like', array);
                
                if shiftX > 0
                    % Image shifted right -> empty space is on the left
                    translated_slice(:, 1:colsToReplace) = noise_cols;
                else
                    % Image shifted left -> empty space is on the right
                    translated_slice(:, (nx - colsToReplace + 1):nx) = noise_cols;
                end
            end
            
            paddedArray(:,:,targetZ) = translated_slice;
        end        
        % 3. TOP PADDING (Slices zPad+nz+1 to totalZ)
        for i = 1:zPad
            targetZ = zPad + nz + i;
            sourceZ_in_orig = nz - i + 1;
            zDistance = targetZ - (zPad + sourceZ_in_orig);
            shiftX = zDistance * baseSlope;
            
            % Apply geometric transformation
            translated_slice = imtranslate(array(:,:,sourceZ_in_orig), [shiftX, 0], 'Method', 'cubic', 'FillValues', 0);
            
            % The previous transformation creates empty regions filled with zeros, now we will fill these with a gaussian noise to better simulate camera noise.
            % Calculate affected columns + 1 extra column
            colsToReplace = ceil(abs(shiftX)) + 1;            
            if colsToReplace > 0
                % Generate Gaussian noise
                noise_cols = cast(config.gaussian_mean + config.gaussian_std .* randn(ny, colsToReplace), 'like', array);                
                if shiftX > 0
                    % Image shifted right -> empty space is on the left
                    translated_slice(:, 1:colsToReplace) = noise_cols;
                else
                    % Image shifted left -> empty space is on the right
                    translated_slice(:, (nx - colsToReplace + 1):nx) = noise_cols;
                end
            end            
            paddedArray(:,:,targetZ) = translated_slice;
        end       
    case 'gaussian'
        frontPad = config.gaussian_mean + config.gaussian_std .* randn(size(array,1), size(array,2), config.z_padding);
        backPad  = config.gaussian_mean + config.gaussian_std .* randn(size(array,1), size(array,2), config.z_padding);
        paddedArray = cat(3, frontPad, array, backPad);
    case 'fixed'
        frontPad = config.fixed_value * ones(size(array,1), size(array,2), config.z_padding);
        backPad  = config.fixed_value * ones(size(array,1), size(array,2), config.z_padding);
        paddedArray = cat(3, frontPad, array, backPad);
    otherwise
        error('Invalid z_edge_padding option: %s', config.z_edge_padding);
end

end

%% -----------------------------------------------------------------------
%% Local Function: runDeconDeskewPipeline
function runDeconDeskewPipeline(seriesResult, config)
    fprintf('Running deconvolution+deskew pipeline for series: %s\n', seriesResult.currentSeriesFolder);
    
    % Deconvolution step.
    XR_decon_data_wrapper(seriesResult.tifDir, 'resultDirName', config.resultDirName, 'xyPixelSize', config.xyPixelSize, ...
                'dz', config.dz, 'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'PSFFullpaths', config.PSFFullpaths, ...
                'dzPSF', config.dzPSF, 'parseSettingFile', config.parseSettingFile, 'RLmethod', config.RLmethod, ...
                'wienerAlpha', config.wienerAlpha, 'OTFCumThresh', config.OTFCumThresh, 'skewed', config.skewed, ...
                'Background', config.Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', config.DeconIter, ...
                'fixIter', config.fixIter, 'EdgeErosion', config.EdgeErosion, 'Save16bit', config.Save16bit, ...
                'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'parseCluster', config.parseCluster, ...
                'largeFile', config.largeFile, 'GPUJob', config.GPUJob, 'debug', config.debug, 'cpusPerTask', config.cpusPerTask, ...
                'ConfigFile', config.ConfigFile, 'GPUConfigFile', config.GPUConfigFile, 'mccMode', config.mccMode);
            
    if config.GPUJob && gpuDeviceCount('available') > 0
         reset(gpuDevice);
    end
    
    % % Remove z-padding from the decon results, before we deskew.
    deconDir = fullfile(seriesResult.tifDir, config.resultDirName);
    removePaddingFromDir(deconDir, config);

    % Deskew step.
    dataPath_exps = fullfile(seriesResult.tifDir, config.resultDirName);
    XR_deskew_rotate_data_wrapper(dataPath_exps, 'skewAngle', config.skewAngle, 'flipZstack', config.flipZstack, ...
        'DSRCombined', config.DSRCombined, 'rotate', config.rotate, 'xyPixelSize', config.xyPixelSize, 'dz', config.dz, ...
        'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'largeFile', config.largeFile, ...
        'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'Save16bit', config.Save16bit, 'parseCluster', config.parseCluster, ...
        'masterCompute', config.masterCompute, 'configFile', config.configFile, 'mccMode', config.mccMode);
    
    % Merge the deconvolved+deskewed images.
    outputTiffFile = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_decondeskew.tif']);
    deconDSDir = fullfile(dataPath_exps, 'DS');
    if isfield(config,'skewDirection')
        paraMergeTiffFilesToMultiDimStack(deconDSDir, outputTiffFile, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeTiffFilesToMultiDimStack(deconDSDir, outputTiffFile, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval);
    end
    
    outputTiffFileMax = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_decondeskew_MAX.tif']);
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

    XR_deskew_rotate_data_wrapper(seriesResult.tifDir, 'resultDirName', config.resultDirNameDeskew, 'skewAngle', config.skewAngle, 'flipZstack', config.flipZstack, ...
        'DSRCombined', config.DSRCombined, 'rotate', config.rotate, 'xyPixelSize', config.xyPixelSize, 'dz', config.dz, ...
        'Reverse', config.Reverse, 'ChannelPatterns', config.ChannelPatterns, 'largeFile', config.largeFile, ...
        'zarrFile', config.zarrFile, 'saveZarr', config.saveZarr, 'Save16bit', config.Save16bit, 'parseCluster', config.parseCluster, ...
        'masterCompute', config.masterCompute, 'configFile', config.configFile, 'mccMode', config.mccMode);
    
    % Merge the deskew-only images.
    outputTiffFileDeskew = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_deskew.tif']);
    deskewDSDir = fullfile(seriesResult.tifDir, config.resultDirNameDeskew);
    if isfield(config,'skewDirection')
        paraMergeTiffFilesToMultiDimStack(deskewDSDir, outputTiffFileDeskew, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeTiffFilesToMultiDimStack(deskewDSDir, outputTiffFileDeskew, seriesResult.pixelSizeX, seriesResult.deskewedZSpacing, seriesResult.frameInterval);
    end
    
    outputTiffFileDeskewMax = fullfile(config.inputFolder, [seriesResult.currentSeriesFolder, '_deskew_MAX.tif']);
    inputToMergeDeskewMax = fullfile(deskewDSDir, 'MIPs');
    if isfield(config,'skewDirection')
        paraMergeMaxToStack(inputToMergeDeskewMax, outputTiffFileDeskewMax, seriesResult.pixelSizeX, seriesResult.frameInterval, config.skewDirection);
    else
        paraMergeMaxToStack(inputToMergeDeskewMax, outputTiffFileDeskewMax, seriesResult.pixelSizeX, seriesResult.frameInterval);
    end
end

%% -----------------------------------------------------------------------
%% Local Function: removePaddingFromDir
function removePaddingFromDir(targetDir, config)
    fileList = dir(fullfile(targetDir, '*.tif'));
    if ~strcmp(config.z_edge_padding, 'none')
        for i = 1:length(fileList)
            filePath = fullfile(targetDir, fileList(i).name);
            img = parallelReadTiff(filePath);
            if size(img, 3) > 2 * config.z_padding
                img_no_padding = img(:, :, (config.z_padding+1):(end-config.z_padding));
                parallelWriteTiff(filePath, img_no_padding);
            else
                warning('Not enough depth for padding removal in file: %s', fileList(i).name);
            end
        end
    end
end

%% -----------------------------------------------------------------------
%% Local Function: deleteIntermediateFiles
function deleteIntermediateFiles(tifDir, config)
    if config.deleteRawTif
         rawFiles = dir(fullfile(tifDir, '*.tif'));
         for k = 1:length(rawFiles)
             delete(fullfile(rawFiles(k).folder, rawFiles(k).name));
         end
    end
    if config.deleteDeconTif
         deconTifDir = fullfile(tifDir, config.resultDirName);
         deconFiles = dir(fullfile(deconTifDir, '*.tif'));
         for k = 1:length(deconFiles)
             delete(fullfile(deconFiles(k).folder, deconFiles(k).name));
         end
    end
end

%% -----------------------------------------------------------------------
%% Local Function: getDefaultConfig
function config = getDefaultConfig()
    % Folder settings (these will be updated by the user selections).
    config.inputFolder = '';   % The folder to deconvolve (set via UI)
    
    % PSF related settings will be updated via UI:
    config.PSFFullpaths = {};  % Cell array of PSF file paths.
    config.ChannelPatterns = {}; % Patterns (e.g., 'Ch1','Ch2', ...) for each channel.
    
    % Imaging parameters.
    config.dz = 0.5;
    config.xyPixelSize = 0.104;
    
    % z–axis padding settings.
    config.z_edge_padding = 'mirror';   % Options: 'none', 'zero', 'mirror', 'gaussian', 'fixed'
    config.z_padding = 20;  % this needs to come with a warning about not padding more than half the size of the stack
    config.gaussian_mean = 102.27;
    config.gaussian_std = 3.17;
    config.fixed_value = 100;
    
    % Flags for deleting intermediate files.
    config.deleteRawTif = false;
    config.deleteDeconTif = false;
    
    % Deconvolution parameters.
    config.RLmethod = 'simplified';
    config.DeconIter = 1;
    config.wienerAlpha = 0.05;
    
    % Acquisition and PSF parameters.
    config.Reverse = true;
    config.dzPSF = 0.5;
    config.parseSettingFile = false;
    % The ChannelPatterns and PSFFullpaths will be updated from the PSF folder.
    config.OTFCumThresh = 0.9;
    config.skewed = true;
    config.Background = 100;
    config.fixIter = true;
    config.EdgeErosion = 0;
    config.Save16bit = true;
    config.zarrFile = false;
    config.saveZarr = false;
    config.cpusPerTask = 4;
    config.parseCluster = false;
    config.largeFile = false;
    config.GPUJob = true;
    config.debug = false;
    config.ConfigFile = '';
    config.GPUConfigFile = '';
    config.mccMode = false;
    
    % Deskew and rotation parameters.
    config.rotate = false;
    config.skewAngle = 32.8;
    config.flipZstack = true;
    config.DSRCombined = false;
    config.masterCompute = true;
    config.configFile = '';
    
    % Processing mode: Choose 'deskew-only', 'decon+deskew', or 'both'.
    config.processingMode = 'decon+deskew';
    
    % Output folder names.
    config.resultDirName = 'deconvolved';       % For deconvolution+deskew branch.
    config.resultDirNameDeskew = 'DS';   % For deskew-only branch.
end

function config = getCziDefaultConfig(config)
    config.xyPixelSize = 0.1449922;
    config.skewAngle = 30.0;
    config.skewDirection = 'Y';
end
