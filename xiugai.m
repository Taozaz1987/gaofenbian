clear; clc; close all;

%% 1. Signal Generation
fs = 1000; 
dt = 1/fs; 
T = 1; 
t = 0:dt:T-dt; 
nt = length(t); 

% Quality Factor and Freq
Q = 80;
f0 = 30; % 30 Hz

% Generate Minimum Phase Ricker Wavelet
wavelet_length = 0.2; 
nw = round(wavelet_length * fs);
tw = (-nw/2:nw/2) * dt;
ricker_zero = (1 - 2*(pi*f0*tw).^2) .* exp(-(pi*f0*tw).^2);
ricker_min = real(hilbert(ricker_zero));
% Normalize wavelet energy
ricker_min = ricker_min / max(abs(ricker_min));

% Place wavelet in center for stationary comparison
w = zeros(1, nt);
center_idx = floor(length(t)/2);
start_idx = center_idx - floor(length(ricker_min)/2);
w(start_idx : start_idx+length(ricker_min)-1) = ricker_min;

% Generate Reflectivity
rng(42); % Fixed seed for reproducibility
reflection_coef = randn(1, nt);
% Sparse reflectivity
reflection_coef(abs(reflection_coef) < 1.0) = 0; 
window_size = 5;
smooth_filter = ones(1, window_size)/window_size;
% reflection_coef = conv(reflection_coef, smooth_filter, 'same'); % Optional smoothing
reflection_coef = reflection_coef / max(abs(reflection_coef)) * 0.5;

%% 2. Generate Non-stationary Signal (Q-attenuation)
nfft = 2^nextpow2(nt);
f_fft = (0:nfft-1) * (fs/nfft);
f_pos = f_fft(1:nfft/2+1); 

% Vectorized construction of non-stationary signal
S_nonstationary = zeros(1, nfft);
W_pos = fft(ricker_min, nfft); W_pos = W_pos(1:nfft/2+1); % FFT of wavelet kernel

fprintf('Generating non-stationary signal...\n');
for tau_idx = 1:nt
    if reflection_coef(tau_idx) == 0, continue; end
    tau = t(tau_idx);
    
    % Constant Q transfer function: exp(-pi*f*t/Q) * exp(i*phase)
    amp_atten = exp(-pi * f_pos * tau / Q);
    
    % Minimum phase calculation for Q filter (Futterman relation)
    % A simple approximation for constant Q phase drift
    ref_f = f0; % reference freq
    phase_lag = 2*pi*f_pos * tau .* (1/pi/Q) .* log((f_pos+eps)/ref_f); 
    
    Q_filter = amp_atten .* exp(-1i * phase_lag);
    
    % Add reflection * wavelet * Q-filter to spectrum
    S_comp = W_pos .* Q_filter * reflection_coef(tau_idx);
    
    % Shift to time tau (linear phase shift)
    S_shifted = S_comp .* exp(-1i * 2 * pi * f_pos * tau); % Note: Wavelet is centered at 0 locally
    
    % Accumulate
    S_nonstationary(1:length(f_pos)) = S_nonstationary(1:length(f_pos)) + S_shifted;
end

% Symmetry for IFFT
if mod(nfft, 2) == 0
    S_nonstationary(nfft/2+2:end) = conj(S_nonstationary(nfft/2:-1:2));
else
    S_nonstationary((nfft+1)/2+1:end) = conj(S_nonstationary((nfft+1)/2:-1:2));
end

s_nonstationary = real(ifft(S_nonstationary, nfft));
s_nonstationary = s_nonstationary(1:nt);
s_nonstationary = s_nonstationary / max(abs(s_nonstationary));

% Stationary Reference
s_stationary = conv(reflection_coef, ricker_min, 'same');
s_stationary = s_stationary / max(abs(s_stationary));

% Plot Input Signals
figure('Position', [100, 100, 800, 500]);
subplot(3,1,1); plot(t, reflection_coef, 'k'); title('Reflectivity'); grid on;
subplot(3,1,2); plot(t, s_stationary, 'b'); title('Stationary Synthetic (No Q)'); grid on;
subplot(3,1,3); plot(t, s_nonstationary, 'r'); title(['Non-stationary (Q=' num2str(Q) ')']); grid on;
xlabel('Time (s)');

%% 3. Adaptive Local Frequency (ALF) Estimation
fprintf('Estimating Local Frequency...\n');
f_vec = linspace(0, 100, 150); % Frequency range for T-F analysis
k_win = 4; % Window scaling factor

% Calculate Adaptive Local Frequency
f_optloc = optloc(s_nonstationary, t); 
% Smooth local frequency for stability in T-F
f_optloc = smooth(f_optloc, 20)'; 

%% 4. ALFMT (Forward Transform)
fprintf('Computing ALFMT...\n');
ALFMT = zeros(length(t), length(f_vec));
sigma_store = zeros(length(t), length(f_vec));

for tau_idx = 1:length(t)
    tau = t(tau_idx);
    f_loc_val = abs(f_optloc(tau_idx)); % Ensure positive freq
    
    % Calculate Window Widths for all frequencies at this time step
    % sigma = k / (f_loc + |f_loc - f|)
    sigma_inv_vec = (f_loc_val + abs(f_loc_val - f_vec)) / k_win; 
    
    % Optimized Inner Loop: Compute Fourier integral for this time step
    % Use a truncated Gaussian window for speed
    for f_idx = 1:length(f_vec)
        cur_f = f_vec(f_idx);
        inv_sig = sigma_inv_vec(f_idx);
        
        % Define effective window range (3 sigma)
        eff_width = 3.5 / (inv_sig + eps); 
        t_start = tau - eff_width;
        t_end = tau + eff_width;
        
        % Indices for window
        idx_range = (t >= t_start) & (t <= t_end);
        if ~any(idx_range), continue; end
        
        t_sub = t(idx_range);
        s_sub = s_nonstationary(idx_range);
        
        % Window function
        win = (inv_sig / sqrt(2*pi)) * exp( -((t_sub - tau).^2) * (inv_sig^2) / 2 );
        
        % Fourier Transform integral
        integrand = s_sub .* win .* exp(-1i * 2 * pi * cur_f * t_sub);
        ALFMT(tau_idx, f_idx) = sum(integrand) * dt; % Trapezoidal approx via sum*dt
    end
end

%% 5. Dynamic Deconvolution
fprintf('Performing Dynamic Deconvolution...\n');

% A. Smoothing the Amplitude Spectrum (|W|)
magnitude_ALFMT = abs(ALFMT);
% Use 2D Gaussian smoothing instead of complex hyperbolic logic for robustness
% Sigma: Time (samples) approx 100ms, Freq (samples) approx 5Hz
sigma_t = round(0.1 / dt); 
sigma_f = round(5 / (f_vec(2)-f_vec(1))); 
try
    ALFMT_s_smooth = imgaussfilt(magnitude_ALFMT, [sigma_t, sigma_f]);
catch
    % Fallback if image toolbox not present
    h_gauss = fspecial_gauss([sigma_t*2, sigma_f*2], sigma_t/2);
    ALFMT_s_smooth = conv2(magnitude_ALFMT, h_gauss, 'same');
end

% B. Minimum Phase Estimation via Hilbert
% Hilbert along frequency axis (dimension 2)
% log(|W|) -> Hilbert -> Imag -> Phase
log_spec = log(ALFMT_s_smooth + eps);
% Note: hilbert in MATLAB operates down columns. We need rows (Freqs).
% Transpose, Hilbert, Transpose.
analytic = hilbert(log_spec.').'; 
phi = imag(analytic);

% C. Deconvolution Division
mu = 0.1; % Pre-whitening / Stabilization factor
A_max = max(ALFMT_s_smooth(:));
Denominator = ALFMT_s_smooth + mu * A_max;

% V_r = V_s / W_est
V_ALFMT_r = ALFMT ./ Denominator .* exp(-1i * phi);

%% 6. Inverse ALFMT (Reconstruction)
fprintf('Reconstructing Signal...\n');

% Collapse Time Axis (Integration over tau)
% R(f) = Integral V_r(tau, f) dtau
R_f_spectrum = sum(V_ALFMT_r, 1) * dt; 

% Inverse Fourier Transform (Integration over f)
% r(t) = Integral R(f) exp(i 2 pi f t) df
% Since we have positive frequencies only, r(t) = 2 * Real(...)
r_reconstruct = zeros(1, nt);
df = f_vec(2) - f_vec(1);

for t_idx = 1:nt
    % Discrete integration over frequency
    integrand = R_f_spectrum .* exp(1i * 2 * pi * f_vec * t(t_idx));
    val = sum(integrand) * df;
    r_reconstruct(t_idx) = 2 * real(val); % Factor of 2 for one-sided spectrum
end

%% 7. Comparison and Plotting
% Convolve recovered reflectivity with stationary wavelet to see improvement
s_recon = conv(r_reconstruct, ricker_min, 'same');
% Normalize
s_recon = s_recon / max(abs(s_recon));

figure('Position', [100, 100, 1000, 600], 'Color', 'w');

subplot(4,1,1);
imagesc(t, f_vec, abs(ALFMT)'); axis xy; colormap(jet);
title('ALFMT Spectrum (Original Attenuated)'); ylabel('Freq (Hz)');

subplot(4,1,2);
imagesc(t, f_vec, abs(V_ALFMT_r)'); axis xy; colormap(jet);
title('ALFMT Spectrum (Deconvolved)'); ylabel('Freq (Hz)');

subplot(4,1,3);
plot(t, f_optloc, 'k', 'LineWidth', 1.5);
title('Estimated Adaptive Local Frequency'); ylabel('Hz');

subplot(4,1,4);
plot(t, s_stationary, 'b:', 'LineWidth', 1.0); hold on;
plot(t, s_nonstationary, 'k-', 'LineWidth', 1.0);
plot(t, s_recon, 'r', 'LineWidth', 1.3);
legend('Stationary (Ideal)', 'Attenuated (Input)', 'Restored (Output)');
title('Trace Comparison'); grid on; xlabel('Time (s)');


%% --- Helper Functions ---

function f_optloc = optloc(x, t)
    % Calculates Adaptive Local Frequency
    % Solves the regularization problem: f = (I + xi^2 P'P)^-1 S* n
    
    dt = t(2) - t(1);
    N = length(x);
    
    % Parameters
    xi = 20;        % Smoothing parameter (higher = smoother)
    epsilon = 1e-6; % Stabilization
    
    % Complex Trace & Envelope
    x_c = hilbert(x);
    h = imag(x_c);
    env = abs(x_c);
    
    % Derivative of phase (Instantaneous Frequency precursor)
    % n(t) = (x h' - x' h) / (2pi)
    x_grad = gradient(x(:), dt);
    h_grad = gradient(h(:), dt);
    
    % Numerator n (unscaled frequency)
    n = (x(:) .* h_grad - x_grad .* h(:)) / (2*pi);
    
    % Denominator d (Envelope squared)
    d = env(:).^2;
    
    % Construct Regularization Matrix P based on Local Energy
    % P weighs the smoothing. High energy -> less smoothing needed?
    % Or P enforces smoothness relative to energy changes.
    
    % Calculate local average energy for P construction
    win = gausswin(50); win = win/sum(win);
    A2_avg = conv(d, win, 'same') + epsilon;
    
    % Build Sparse P matrix (First order difference weighted)
    % P(i, i) = 1/sqrt(E_i), P(i, i+1) = -1/sqrt(E_j)
    inv_root_E = 1 ./ sqrt(A2_avg);
    
    % Using sparse diagonals for efficiency
    % Standard Laplacian smoothing weighted by data energy
    e = ones(N, 1);
    % Simplification: Standard Shaping Regularization
    % Solve: (D + xi^2 D_2' D_2) f = n
    % where D is data weight (envelope), D_2 is roughness.
    
    % Implementing logic similar to Equation 15/17 in paper:
    % f = (lambda^2 I + S*(D - lambda^2 I)) \ S* n
    % Simplified logic for robustness: Weighted Least Squares
    % Minimize || d*f - n ||^2 + xi^2 || L f ||^2
    
    L = spdiags([-e e], [0 1], N-1, N); % Derivative matrix
    W_data = spdiags(d + epsilon, 0, N, N); % Data weights
    
    % System: (W_data^2 + xi^2 * L'*L) f = W_data * n
    % Note: Use W_data (envelope) to trust high energy areas more
    
    A_sys = W_data + xi^2 * (L' * L);
    rhs = n;
    
    f_optloc = A_sys \ rhs;
    
    % Bound frequencies to reasonable physical limits
    f_optloc(f_optloc < 0) = 0;
    f_optloc(f_optloc > 200) = 200; % Max plausible freq
end

function h = fspecial_gauss(hsize, sigma)
    % Fallback for Gaussian kernel if Image Toolbox is missing
    siz = (hsize-1)/2;
    [x,y] = meshgrid(-siz(2):siz(2), -siz(1):siz(1));
    arg = -(x.*x + y.*y)/(2*sigma*sigma);
    h = exp(arg);
    h(h<eps*max(h(:))) = 0;
    sumh = sum(h(:));
    if sumh ~= 0, h = h/sumh; end
end