function RLGC(image_filename, psf_filename, output_filename, bg_val, save16bit, seed)
% RLGC, Richardson-Lucy deconvolution using gradient consensus to stop iterations locally.
%
% Not tested or verified! 
% May not be an accurate representation of the original python version
% MATLAB adaptation by David Corcoran, 2026
%
% Original by James Manton, 2023
% jmanton@mrc-lmb.cam.ac.uk
% Developed in collaboration with Andy York (Calico), Jan Becker (Oxford) and Craig Russell (EMBL EBI)
% https://zenodo.org/records/10278919
% https://github.com/jdmanton/rlgc/blob/master/rlgc.py
%
% Performance modifications specific to matlab:
%   - Exact grouped scalar-n binomial split:
%       split1(x) ~ Binomial(floor(img(x)), 0.5)
%       split2    = img - split1
%       Instead of using binornd(img_int, 0.5) as heterogeneous array-valued
%       n can be dramatically slower than repeated exact scalar-n calls grouped by n.%
%   - GPU-aware RNG seeding with gpurng
%   - Reused normalized/logged Hu for KL divergence
%   - Branchless update to avoid masked indexed assignment
%
% Inputs:
%   image_filename  - input TIFF stack
%   psf_filename    - PSF TIFF stack (needs to be processed)
%   output_filename - output deconvolved TIFF stack
%   bg_val          - scalar background to subtract before deconvolution (default 0)
%   save16bit       - true to save uint16, false to save single (default true)
%   seed            - RNG seed (default 42)

if nargin < 6; seed = 42; end
if nargin < 5; save16bit = true; end
if nargin < 4; bg_val = 0; end

%% ----------------------------
%  Load image
% -----------------------------
disp(['Loading data for RLGC: ', image_filename]);

img_info = imfinfo(image_filename);
nImgZ = numel(img_info);

img = zeros(img_info(1).Height, img_info(1).Width, nImgZ, 'single');
for k = 1:nImgZ
    img(:, :, k) = single(imread(image_filename, k));
end

if bg_val > 0
    img = max(img - single(bg_val), 0);
end

% Convert MATLAB (Y, X, Z) -> internal (Z, Y, X)
img = permute(img, [3, 1, 2]);

%% ----------------------------
%  Load PSF
% -----------------------------
psf_info = imfinfo(psf_filename);
nPsfZ = numel(psf_info);

psf_temp = zeros(psf_info(1).Height, psf_info(1).Width, nPsfZ, 'single');
for k = 1:nPsfZ
    psf_temp(:, :, k) = single(imread(psf_filename, k));
end

% Convert MATLAB (Y, X, Z) -> internal (Z, Y, X)
psf_temp = permute(psf_temp, [3, 1, 2]);

% Pad PSF to image size
psf = zeros(size(img), 'single');
sz_p = size(psf_temp);
psf(1:sz_p(1), 1:sz_p(2), 1:sz_p(3)) = psf_temp;

% Center the PSF to mimic original Python/NumPy behaviour
for i = 1:ndims(psf)
    psf = circshift(psf, floor(size(psf, i) / 2), i);
end
for i = 1:ndims(psf_temp)
    psf = circshift(psf, -floor(size(psf_temp, i) / 2), i);
end

% Inverse FFT shift and normalize
psf = ifftshift(psf);
psf = psf / sum(psf, 'all');

%% ----------------------------
%  GPU setup
% -----------------------------
use_gpu = canUseGPU();

% Seed CPU RNG in case any CPU-side randomness is used
rng(seed);

if use_gpu
    img = gpuArray(img);
    psf = gpuArray(psf);

    % GPU RNG is separate from CPU RNG
    gpurng(seed, "Philox");
end

%% ----------------------------
%  Precompute OTFs
% -----------------------------
otf = fftn(psf);
otfT = conj(otf);
otfotfT = otf .* otfT;

clear psf psf_temp;

%% ----------------------------
%  Precompute exact split ingredients
% -----------------------------
% Preserve original split semantics:
%   split1 = binornd(floor(img), 0.5)
%   split2 = img - split1
img_int = floor(img);
img_frac = img - img_int;

if isgpuarray(img_int)
    img_int_cpu = gather(img_int);
else
    img_int_cpu = img_int;
end

lin_nonzero = find(img_int_cpu > 0);
vals_nonzero = double(img_int_cpu(lin_nonzero));
max_val = max(vals_nonzero);

idx_groups_all = accumarray( ...
    vals_nonzero, ...
    lin_nonzero, ...
    [max_val, 1], ...
    @(x){x}, ...
    {[]} );

nonempty = ~cellfun(@isempty, idx_groups_all);
group_n = find(nonempty);
idx_groups = idx_groups_all(nonempty);

clear img_int_cpu;

%% ----------------------------
%  Initialize reconstruction
% -----------------------------
recon = mean(img, 'all') .* ones(size(img), 'like', img);
previous_recon = recon;

num_iters = 0;
prev_kld1 = inf;
prev_kld2 = inf;
start_time = tic;

%% ----------------------------
%  Main RLGC loop
% -----------------------------
while true
    iter_start_time = tic;

    % Exact grouped scalar-n binomial split
    split1 = grouped_binomial_split(img_int, idx_groups, group_n, use_gpu);
    split2 = (img_int - split1) + img_frac;

    % Forward model
    Hu = fftconv(recon, otf);

    % Reuse normalized/logged Hu for all KLD evaluations
    [Hu_norm, logHu] = prepare_kld_p(Hu);
    kld1 = kldiv_prepared(Hu_norm, logHu, split1);
    kld2 = kldiv_prepared(Hu_norm, logHu, split2);

    % Stop when both split-image KLDs have increased
    if (kld1 > prev_kld1) && (kld2 > prev_kld2)
        recon = previous_recon;
        fprintf('Optimum result obtained after %d iterations with a total time of %1.1f seconds.\n', ...
            num_iters - 1, toc(start_time));
        break;
    end

    prev_kld1 = kld1;
    prev_kld2 = kld2;

    % Backprojection terms
    denom = 0.5 * (Hu + 1e-12);
    HTratio1 = fftconv(split1 ./ denom, otfT);
    HTratio2 = fftconv(split2 ./ denom, otfT);
    HTratio  = HTratio1 + HTratio2;

    % Save previous estimate
    previous_recon = recon;

    % Gradient consensus update (branchless)
    consensus_map = fftconv((HTratio1 - 1) .* (HTratio2 - 1), otfotfT);
    mask = cast(consensus_map >= 0, 'like', recon);
    recon = recon .* (1 + (HTratio - 1) .* mask);

    calc_time = toc(iter_start_time);
    fprintf('Iteration %03d completed in %1.3f s.\n', num_iters + 1, calc_time);
    num_iters = num_iters + 1;

    % Optional cleanup of temporaries
    clear split1 split2 Hu Hu_norm logHu denom HTratio1 HTratio2 HTratio consensus_map mask;
end

%% ----------------------------
%  Gather and save result
% -----------------------------
if isgpuarray(recon)
    recon = gather(recon);
end

% Convert internal (Z, Y, X) -> MATLAB save order (Y, X, Z)
recon_save = permute(recon, [2, 3, 1]);

t = Tiff(output_filename, 'w8');

tags.ImageLength = size(recon_save, 1);
tags.ImageWidth = size(recon_save, 2);
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.SamplesPerPixel = 1;
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

if save16bit
    recon_save = uint16(recon_save);
    tags.BitsPerSample = 16;
    tags.SampleFormat = Tiff.SampleFormat.UInt;
else
    tags.BitsPerSample = 32;
    tags.SampleFormat = Tiff.SampleFormat.IEEEFP;
end

for k = 1:size(recon_save, 3)
    t.setTag(tags);
    t.write(recon_save(:, :, k));
    if k < size(recon_save, 3)
        t.writeDirectory();
    end
end
t.close();
end

%% ============================================================
% Local helper functions
% =============================================================

function split1 = grouped_binomial_split(img_int, idx_groups, group_n, use_gpu)
% Exact grouped scalar-n binomial split.
%
% For each unique n > 0:
%   vals ~ Binomial(n, 0.5), size = [count_of_pixels_with_that_n, 1]
%   scatter vals back into split1
%
% Preserves exact Binomial(n, 0.5) marginals for each voxel.

split1 = zeros(size(img_int), 'like', img_int);

for i = 1:numel(group_n)
    n = group_n(i);
    c = numel(idx_groups{i});

    if c == 0
        continue;
    end

    if use_gpu
        % Give binornd a gpuArray data input so it executes on the GPU
        vals = single(binornd(gpuArray(single(n)), 0.5, c, 1));
    else
        vals = single(binornd(single(n), 0.5, c, 1));
    end

    split1(idx_groups{i}) = vals;
end
end

function out = fftconv(x, H)
% FFT-based circular convolution
out = real(ifftn(fftn(x) .* H));
end

function [p_norm, logp] = prepare_kld_p(p)
% Prepare normalized p and log(p) once for reuse in multiple KLD evaluations
p = p + 1e-4;
p_norm = p / sum(p, 'all');
logp = log(p_norm);
end

function kld = kldiv_prepared(p_norm, logp, q)
% KLD using precomputed normalized/logged p
q = q + 1e-4;
q = q / sum(q, 'all');
logq = log(q);

kld_map = p_norm .* (logp - logq);
kld_map(isnan(kld_map)) = 0;
kld = sum(kld_map, 'all');
end