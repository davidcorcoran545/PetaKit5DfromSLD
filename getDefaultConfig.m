%% getDefaultConfig
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
    config.z_padding = 30;  % this needs to come with a warning about not padding more than half the size of the stack
    config.gaussian_mean = 102.27;
    config.gaussian_std = 3.17;
    config.fixed_value = 100;
    
    % Flags for deleting intermediate files.
    config.deleteRawTif = true;
    config.deleteDeconTif = true;
    
    % Deconvolution parameters.
    config.deconAlgorithm = 'PetaKit5D'; % Options: 'petaKit5D', 'RLGC'
    config.RLmethod = 'simplified';
    config.DeconIter = 10;
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

