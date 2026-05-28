%% launchPipelineGUI
function [uiResult, canceled] = launchPipelineGUI()
    % Setup defaults and initialize outputs
    canceled = true;
    uiResult = struct();

    % Load last used UI state safely
    uiState = loadLastUIState();

    defaultPSF   = uiState.psfFolder;
    defaultInput = uiState.inputFolder;
    cfg          = uiState.config;
    
    % Try to load output folder safely, otherwise default to input folder
    if isfield(uiState, 'config') && isfield(uiState.config, 'outputFolder')
        defaultOutput = uiState.config.outputFolder;
    else
        defaultOutput = defaultInput;
    end

    cfg          = uiState.config;

    % Increase the below if you add more features to give more room
    fig = uifigure('Name', 'Modular Pipeline Configuration', 'Position', [100, 100, 520, 700]);

    % --- 1. Folders ---
    yPos = 670; 
    uilabel(fig, 'Position', [20, yPos, 400, 22], 'Text', '1. Select Folders', 'FontWeight', 'bold');
    
    yPos = yPos - 30;
    btnInput = uibutton(fig, 'Position', [20, yPos, 120, 22], 'Text', 'Select Input Folder');
    lblInput = uilabel(fig, 'Position', [150, yPos, 310, 22], 'Text', defaultInput, 'Interpreter', 'none');
    btnInput.ButtonPushedFcn = @(btn,event) folderSelect(lblInput, 'Select Input Folder');

    yPos = yPos - 30;
    btnOutput = uibutton(fig, 'Position', [20, yPos, 120, 22], 'Text', 'Select Output Folder');
    lblOutput = uilabel(fig, 'Position', [150, yPos, 310, 22], 'Text', defaultOutput, 'Interpreter', 'none');
    btnOutput.ButtonPushedFcn = @(btn,event) folderSelect(lblOutput, 'Select Output Folder');

    yPos = yPos - 30;
    btnPSF = uibutton(fig, 'Position', [20, yPos, 120, 22], 'Text', 'Select PSF Folder');
    lblPSF = uilabel(fig, 'Position', [150, yPos, 310, 22], 'Text', defaultPSF, 'Interpreter', 'none');
    btnPSF.ButtonPushedFcn = @(btn,event) folderSelect(lblPSF, 'Select PSF Folder');

    % --- 2. Parameters ---
    yPos = yPos - 40;
    uilabel(fig, 'Position', [20, yPos, 400, 22], 'Text', '2. Metadata', 'FontWeight', 'bold');

    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 200, 22], 'Text', 'Image XY pixel size (microns):');
    editXY = uieditfield(fig, 'numeric', 'Position', [250, yPos, 80, 22], 'Value', cfg.xyPixelSize);

    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 200, 22], 'Text', 'Image Z-step size (microns):');
    editDz = uieditfield(fig, 'numeric', 'Position', [250, yPos, 80, 22], 'Value', cfg.dz);
    
    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 200, 22], 'Text', 'PSF Z-step size (microns):');
    editDzPSF = uieditfield(fig, 'numeric', 'Position', [250, yPos, 80, 22], 'Value', cfg.dzPSF);

    % --- 3. Processing Mode ---
    yPos = yPos - 40;
    uilabel(fig, 'Position', [20, yPos, 400, 22], ...
        'Text', '3. Processing Mode', 'FontWeight', 'bold');

    yPos = yPos - 25;
    chkDeskew = uicheckbox(fig, 'Position', [30, yPos, 150, 22], ...
        'Text', 'Deskew', ...
        'Value', strcmp(cfg.processingMode, 'deskew-only') || strcmp(cfg.processingMode, 'both'));

    yPos = yPos - 25;
    chkDecon = uicheckbox(fig, 'Position', [30, yPos, 220, 22], ...
        'Text', 'Deconvolve and then deskew', ...
        'Value', strcmp(cfg.processingMode, 'decon+deskew') || strcmp(cfg.processingMode, 'both'));

    yPos = yPos - 25;
    chkRotate = uicheckbox(fig, 'Position', [30, yPos, 250, 22], ...
        'Text', 'Coverslip correction (rotation)', ...
        'Value', cfg.rotate);

    % --- 4. Deconvolution Settings ---
    yPos = yPos - 40;
    uilabel(fig, 'Position', [20, yPos, 400, 22], 'Text', '4. Deconvolution Settings', 'FontWeight', 'bold');

    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 220, 22], 'Text', 'Deconvolution method:');
    ddAlgorithm = uidropdown(fig, 'Position', [250, yPos, 100, 22], ...
        'Items', {'PetaKit5D', 'RLGC'}, 'Value', cfg.deconAlgorithm);

    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 220, 22], 'Text', 'Number of Deconvolution Iterations:');
    editIter = uieditfield(fig, 'numeric', 'Position', [250, yPos, 80, 22], 'Value', cfg.DeconIter);
    
    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 220, 22], 'Text', 'Background value to subtract:');
    editBG = uieditfield(fig, 'numeric', 'Position', [250, yPos, 80, 22], 'Value', cfg.Background);
   
    yPos = yPos - 30;
    chkZPad = uicheckbox(fig, 'Position', [30, yPos, 220, 22], ...
        'Text', 'Pad image in the z-axis', ...
        'Value', ~strcmpi(cfg.z_edge_padding, 'none'));

    yPos = yPos - 30;
    uilabel(fig, 'Position', [20, yPos, 220, 22], ...
        'Text', 'Number of slices to pad:');

    editZPad = uieditfield(fig, 'numeric', ...
        'Position', [250, yPos, 80, 22], ...
        'Value', cfg.z_padding, ...
        'Limits', [0 Inf], ...
        'RoundFractionalValues', true);
    editZPad.Enable = matlab.lang.OnOffSwitchState(chkZPad.Value);
    chkZPad.ValueChangedFcn = @(src,event) set(editZPad, 'Enable', matlab.lang.OnOffSwitchState(src.Value))

    % --- 5. File Management ---
    yPos = yPos - 40;
    uilabel(fig, 'Position', [20, yPos, 400, 22], ...
        'Text', '5. File Management', 'FontWeight', 'bold');

    yPos = yPos - 25;
    chkDelInter = uicheckbox(fig, 'Position', [30, yPos, 320, 22], ...
        'Text', 'Delete intermediate processed TIFs', ...
        'Value', cfg.deleteDeconTif || cfg.deleteRawTif);

    % --- Run Button ---
    btnRun = uibutton(fig, 'Position', [190, 20, 100, 35], 'Text', 'Run Pipeline', ...
        'ButtonPushedFcn', @runPipeline);

    drawnow;    
    uiwait(fig);

    % --- Callbacks ---
    function folderSelect(lblTarget, promptTitle)
        startPath = lblTarget.Text;
        if ~isfolder(startPath), startPath = pwd; end
        selectedFolder = uigetdir(startPath, promptTitle);
        if selectedFolder ~= 0
            lblTarget.Text = selectedFolder;
        end
    end

    function runPipeline(~, ~)
        baseConfig = getDefaultConfig();        
        baseConfig.inputFolder = lblInput.Text;
        baseConfig.outputFolder = lblOutput.Text;
        uiResult.psfFolder     = lblPSF.Text;

        % Map GUI values
        baseConfig.deconAlgorithm = ddAlgorithm.Value;
        baseConfig.dz             = editDz.Value;
        baseConfig.xyPixelSize    = editXY.Value;
        baseConfig.dzPSF          = editDzPSF.Value;
        baseConfig.DeconIter      = editIter.Value;
        baseConfig.Background     = editBG.Value;

        % z-axis padding settings
        baseConfig.z_padding = editZPad.Value;
        if chkZPad.Value
            baseConfig.z_edge_padding = 'mirror';
        else
            baseConfig.z_edge_padding = 'none';
        end

        % File deletion logic:
        % User clicking delete intermediates will delete both the initial
        % TIFs converted from the sld/czi files and the intermediate
        % deconvolved and deskewed tifs
        baseConfig.deleteDeconTif = chkDelInter.Value;
        baseConfig.deleteRawTif   = chkDelInter.Value;


        if chkDeskew.Value && chkDecon.Value
            baseConfig.processingMode = 'both';
        elseif chkDecon.Value
            baseConfig.processingMode = 'decon+deskew';
        else
            baseConfig.processingMode = 'deskew-only';
        end

        baseConfig.rotate = chkRotate.Value;

        uiResult.config = baseConfig;

        % Save UI state safely
        saveLastUIState(uiResult.psfFolder, lblInput.Text, baseConfig);

        canceled = false;
        delete(fig);
    end
end

function uiState = loadLastUIState()
    % Always start from known-good defaults
    uiState = struct();
    uiState.psfFolder   = pwd;
    uiState.inputFolder = pwd;
    uiState.config      = getDefaultConfig();

    try
        if ispref('ModPipeline', 'lastUIState')
            saved = getpref('ModPipeline', 'lastUIState');

            if isstruct(saved)
                % Folders
                if isfield(saved, 'psfFolder') && isTextScalar(saved.psfFolder) && isfolder(saved.psfFolder)
                    uiState.psfFolder = char(saved.psfFolder);
                end

                if isfield(saved, 'inputFolder') && isTextScalar(saved.inputFolder) && isfolder(saved.inputFolder)
                    uiState.inputFolder = char(saved.inputFolder);
                    uiState.config.inputFolder = char(saved.inputFolder);
                end

                % Config
                if isfield(saved, 'config') && isstruct(saved.config)
                    uiState.config = mergeSavedConfig(uiState.config, saved.config);
                end
            end
        else
            % Backward compatibility with your existing saved folder prefs
            try
                if ispref('ModPipeline', 'lastPSF')
                    tmp = getpref('ModPipeline', 'lastPSF');
                    if isTextScalar(tmp) && isfolder(tmp)
                        uiState.psfFolder = char(tmp);
                    end
                end
                if ispref('ModPipeline', 'lastInput')
                    tmp = getpref('ModPipeline', 'lastInput');
                    if isTextScalar(tmp) && isfolder(tmp)
                        uiState.inputFolder = char(tmp);
                        uiState.config.inputFolder = char(tmp);
                    end
                end
            catch
                % Ignore and keep defaults
            end
        end
    catch
        % If anything goes wrong, silently keep defaults
        uiState = struct();
        uiState.psfFolder   = pwd;
        uiState.inputFolder = pwd;
        uiState.config      = getDefaultConfig();
    end

    % Normalise dropdown value so it always matches UI items exactly
    if strcmpi(uiState.config.deconAlgorithm, 'petakit5d')
        uiState.config.deconAlgorithm = 'PetaKit5D';
    elseif strcmpi(uiState.config.deconAlgorithm, 'rlgc')
        uiState.config.deconAlgorithm = 'RLGC';
    else
        uiState.config.deconAlgorithm = 'PetaKit5D';
    end
end

function saveLastUIState(psfFolder, inputFolder, config)
    try
        saved = struct();
        saved.psfFolder   = char(psfFolder);
        saved.inputFolder = char(inputFolder);
        saved.config      = sanitizeConfigForSave(config);

        setpref('ModPipeline', 'lastUIState', saved);

        % Optional: keep old prefs too for compatibility/debugging
        setpref('ModPipeline', 'lastPSF', saved.psfFolder);
        setpref('ModPipeline', 'lastInput', saved.inputFolder);
    catch
        % Intentionally do nothing:
        % failure to save preferences should never block processing
    end
end

function cfg = mergeSavedConfig(defaultCfg, savedCfg)
    cfg = defaultCfg;

    % Numeric scalars
    cfg.dz          = pickNumeric(savedCfg, 'dz',          defaultCfg.dz,          @(x) isfinite(x) && x > 0);
    cfg.xyPixelSize = pickNumeric(savedCfg, 'xyPixelSize', defaultCfg.xyPixelSize, @(x) isfinite(x) && x > 0);
    cfg.dzPSF       = pickNumeric(savedCfg, 'dzPSF',       defaultCfg.dzPSF,       @(x) isfinite(x) && x > 0);
    cfg.DeconIter   = pickNumeric(savedCfg, 'DeconIter',   defaultCfg.DeconIter,   @(x) isfinite(x) && x >= 1 && mod(x,1)==0);
    cfg.Background  = pickNumeric(savedCfg, 'Background',  defaultCfg.Background,  @(x) isfinite(x) && x >= 0);
    cfg.z_padding   = pickNumeric(savedCfg, 'z_padding',   defaultCfg.z_padding,   @(x) isfinite(x) && x >= 0 && mod(x,1)==0);

    % Logicals
    cfg.deleteRawTif   = pickLogical(savedCfg, 'deleteRawTif',   defaultCfg.deleteRawTif);
    cfg.deleteDeconTif = pickLogical(savedCfg, 'deleteDeconTif', defaultCfg.deleteDeconTif);
    cfg.rotate         = pickLogical(savedCfg, 'rotate',         defaultCfg.rotate);

    % Enums / strings
    cfg.processingMode = pickEnum(savedCfg, 'processingMode', ...
        {'deskew-only','decon+deskew','both'}, defaultCfg.processingMode);

    cfg.deconAlgorithm = pickEnum(savedCfg, 'deconAlgorithm', ...
        {'PetaKit5D','RLGC','petaKit5D'}, defaultCfg.deconAlgorithm);

    cfg.z_edge_padding = pickEnum(savedCfg, 'z_edge_padding', ...
        {'none','zero','mirror','gaussian','fixed'}, defaultCfg.z_edge_padding);

    % Normalise
    if strcmpi(cfg.deconAlgorithm, 'petakit5d')
        cfg.deconAlgorithm = 'PetaKit5D';
    end
end

function out = sanitizeConfigForSave(cfg)
    % Start from defaults, then merge valid values from cfg
    out = mergeSavedConfig(getDefaultConfig(), cfg);
end

function tf = isTextScalar(x)
    tf = (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x));
end

function val = pickNumeric(s, fieldName, defaultVal, validator)
    val = defaultVal;
    try
        if isfield(s, fieldName)
            x = s.(fieldName);
            if isnumeric(x) && isscalar(x) && validator(x)
                val = double(x);
            end
        end
    catch
    end
end

function val = pickLogical(s, fieldName, defaultVal)
    val = defaultVal;
    try
        if isfield(s, fieldName)
            x = s.(fieldName);
            if islogical(x) && isscalar(x)
                val = x;
            elseif isnumeric(x) && isscalar(x) && any(x == [0 1])
                val = logical(x);
            end
        end
    catch
    end
end

function val = pickEnum(s, fieldName, validValues, defaultVal)
    val = defaultVal;
    try
        if isfield(s, fieldName)
            x = s.(fieldName);
            if ischar(x) || (isstring(x) && isscalar(x))
                x = char(x);
                if any(strcmpi(x, validValues))
                    % Return canonical spelling from validValues if possible
                    idx = find(strcmpi(x, validValues), 1, 'first');
                    val = validValues{idx};
                end
            end
        end
    catch
    end
end