function seismic_alfmt_demo()
% SEISMIC_ALFMT_DEMO Reproduces the ALFMT method and Dynamic Deconvolution
%
% This script implements:
% 1. Adaptive Local Frequency (ALF) estimation.
% 2. Adaptive Local Frequency Modulation Transform (ALFMT).
% 3. Dynamic Deconvolution based on ALFMT.
% 4. Synthetic data tests (FM signal and Seismic trace).
%
% Reference: Huang et al., "Dynamic deconvolution based on adaptive local 
% frequency modulation transform", Petroleum Science, 2025.

    %% 1. Synthetic FM Signal Test (Verify ALFMT)
    fprintf('Running FM Signal Test...\n');
    fs = 1000; % Sampling frequency (Hz)
    T = 1;     % Duration (s)
    t = 0:1/fs:T-1/fs;
    N = length(t);
    
    % Eq. 33: Nonlinear frequency modulation signal
    % s(t) = 0.6 * sin(2*pi*(120*t + 2.5*sin(20*t)))
    % Note: Phase is integral of freq. The formula in paper is likely the phase argument directly.
    % If freq is 120 + ..., phase is 120t + ... 
    % We strictly follow Eq 33 expression as the signal S.
    sig_fm = 0.6 * sin(2*pi * (120*t + 2.5*sin(20*t)));
    
    % Compute Ground Truth Instantaneous Frequency (for comparison)
    % f(t) = 1/(2pi) * d/dt (phase) = 120 + 2.5*20*cos(20t) = 120 + 50*cos(20t)
    f_true = 120 + 50 * cos(20*t);

    % Estimate Adaptive Local Frequency
    % Parameters for regularization (approximate based on paper visual)
    lambda_params.epsilon = 1e-6; % Damping
    lambda_params.xi = 10;        % Smoothing factor (xi)
    f_loc = estimate_local_frequency(sig_fm, fs, lambda_params);
    
    % Perform ALFMT
    freq_vec = 0:1:200; % Frequency range
    k_alfmt = 1;        % Scalar factor k (usually 1)
    [spec_alfmt, ~] = alfmt_forward(sig_fm, t, freq_vec, f_loc, k_alfmt);
    
    % Plot FM results
    figure('Name', 'ALFMT FM Signal Test', 'Color', 'w');
    subplot(3,1,1); plot(t, sig_fm); title('Original Signal'); xlabel('Time (s)');
    subplot(3,1,2); 
    plot(t, f_true, 'r--', 'LineWidth', 1.5); hold on;
    plot(t, f_loc, 'b', 'LineWidth', 1);
    legend('True IF', 'Adaptive Local Freq'); title('Frequency Estimation'); xlabel('Time (s)'); ylabel('Hz');
    subplot(3,1,3);
    imagesc(t, freq_vec, abs(spec_alfmt)); axis xy; colormap(jet);
    title('ALFMT Amplitude Spectrum'); xlabel('Time (s)'); ylabel('Frequency (Hz)');
    
    %% 2. Synthetic Seismic Dynamic Deconvolution
    fprintf('Running Dynamic Deconvolution Test...\n');
    
    % Simulation Parameters
    dt = 0.002;      % 2ms sampling (500Hz)
    len_t = 1.0;     % 1s record
    t_seis = 0:dt:len_t-dt;
    Q = 60;          % Quality factor
    f_dom = 30;      % Ricker dominant freq
    
    % Generate Reflectivity (Sparse spikes)
    r = zeros(size(t_seis));
    r(100) = 0.5; r(200) = -0.3; r(300) = 0.8; r(400) = -0.5;
    
    % Generate Non-stationary Seismic Signal (Forward Q Model)
    seis_att = generate_nonstationary_seismic(r, t_seis, f_dom, Q);
    
    % Apply ALFMT Dynamic Deconvolution
    % 1. Estimate Local Frequency of seismic trace
    f_loc_seis = estimate_local_frequency(seis_att, 1/dt, lambda_params);
    
    % 2. Dynamic Deconvolution
    % Parameters: mu (stabilization), smoothing window sizes
    decon_params.mu = 0.01; 
    decon_params.freq_vec = 0:1:150;
    decon_params.k = 1;
    
    [r_est, seis_rec, spec_orig, spec_restored] = alfmt_dynamic_deconvolution(seis_att, t_seis, f_loc_seis, decon_params);
    
    % Plot Deconvolution Results
    figure('Name', 'Dynamic Deconvolution Results', 'Color', 'w');
    subplot(3,1,1); plot(t_seis, r); title('Reflectivity');
    subplot(3,1,2); plot(t_seis, seis_att); title('Attenuated Seismic Signal (Input)');
    subplot(3,1,3); plot(t_seis, seis_rec, 'r'); title('Compensated Signal (Output)');
    xlabel('Time (s)');
    
    figure('Name', 'Time-Frequency Comparison', 'Color', 'w');
    subplot(2,1,1); imagesc(t_seis, decon_params.freq_vec, abs(spec_orig)); axis xy; colormap(jet);
    title('Original Spectrum'); ylabel('Hz');
    subplot(2,1,2); imagesc(t_seis, decon_params.freq_vec, abs(spec_restored)); axis xy; colormap(jet);
    title('Restored Spectrum'); ylabel('Hz'); xlabel('Time (s)');

end

%% Core Algorithms

function [f_opt_loc] = estimate_local_frequency(x, fs, params)
    % ESTIMATE_LOCAL_FREQUENCY Computes adaptive local frequency (Eqs. 9-17)
    % Inputs:
    %   x: Signal vector
    %   fs: Sampling frequency
    %   params: Struct with xi (smoothness), epsilon (damping)
    
    N = length(x);
    x = x(:);
    
    % 1. Complex trace and Envelope (Eq. 6-8)
    h_x = imag(hilbert(x)); % Hilbert transform (imaginary part)
    env_sq = x.^2 + h_x.^2; % Envelope squared (A^2)
    
    % 2. Instantaneous Frequency (Eq. 9)
    % Discrete derivative approximation for inst freq
    % f = (1/2pi) * (x * h' - x' * h) / (x^2 + h^2)
    % Use central differences for derivatives
    dx = gradient(x) * fs;
    dh = gradient(h_x) * fs;
    
    num = x .* dh - dx .* h_x;
    den = env_sq + params.epsilon; % Add epsilon to avoid division by zero
    
    % Unsmoothed instantaneous frequency vectors n and d (Eq. 10)
    n_vec = num / (2*pi);
    d_vec = den; % This corresponds to the weight D in the paper
    
    % 3. Construct Regularization Operator S* (Eq. 13)
    % P is difference matrix based on local energy (Eq. 14)
    % We compute P as a sparse matrix
    
    % Local energy average <A^2>
    % Use a small smoothing window to estimate local energy average for P calculation
    win_len = round(N/20); if win_len<3, win_len=3; end
    avg_E = smooth(env_sq, win_len); 
    
    % Construct diagonals of P
    % P(i,i) = - sqrt(E(i)/E(i)) = -1
    % P(i,i+1) = sqrt(E(i+1)/E(i))
    ratio = sqrt([avg_E(2:end); avg_E(end)] ./ (avg_E + params.epsilon));
    
    % P is (N-1)xN matrix, approximated here as NxN for simplicity in regularization
    % Standard first order difference with weights
    e = ones(N,1);
    P_diag0 = -ones(N,1);
    P_diag1 = ratio; % super-diagonal
    
    P = spdiags([P_diag0, P_diag1], [0, 1], N, N);
    % Remove last row artifact if needed or keep square
    P(N, :) = 0; 
    
    % S* = (I + xi^2 P'P)^-1
    xi = params.xi;
    I = speye(N);
    S_inv = I + xi^2 * (P' * P);
    
    % 4. Local Adaptive Function lambda(t) (Eq. 16)
    % Local rms of data
    lambda_t = sqrt(smooth(x.^2 + h_x.^2, round(fs*0.1)) + params.epsilon); 
    % Note: Paper defines lambda based on sum. We assume normalized.
    Lambda = spdiags(lambda_t.^2, 0, N, N);
    
    % 5. Solve for f_opt_loc (Eq. 17)
    % f_opt = (Lambda^2*I + S*(D - Lambda^2*I))^-1 * S* * n
    % This equation involves S* explicitly. 
    % To implement efficiently:
    % Let y = S* * n  => Solve (I + xi^2 P'P) y = n
    % Let M = S*(D - Lambda) => Column-wise solve?
    %
    % Simplification: The paper formulates this to balance data fitting (D) and smoothness.
    % Standard Shaping Regularization: f = S [ n/d ] is naive.
    % The formula Eq 17 incorporates the data weight D and the adaptive lambda.
    % Rearranging Eq 17 is complex. We use a Conjugate Gradient or direct solve if N is small.
    %
    % Let's assume the standard Formel shaping logic:
    % Model: D f ~ n
    % Regularization: f is smooth (S applied).
    % Inversion: (D + R) f = n. 
    % Here we solve the specific linear algebra defined.
    
    % Calculate y = S* * n
    y = S_inv \ n_vec;
    
    % Calculate Operator A = Lambda + S*(D - Lambda)
    % Since S* is dense, we don't form it. We assume Eq 17 implies:
    % f = (Lambda + S*(D - Lambda)) \ y ?? 
    % Actually Eq 17: f = (Lambda + S*(D - Lambda))^-1 * S* * n
    % Let RHS = y.
    % LHS matrix A. We need A * f = y.
    % A = Lambda + S* * (D - Lambda).
    % We can solve A f = y using iterative method (e.g. gmres)
    % matvec operation: v -> Lambda*v + S_inv \ ((D - Lambda)*v)
    
    D_mat = spdiags(d_vec, 0, N, N);
    D_minus_L = D_mat - Lambda;
    
    afun = @(v) Lambda*v + (S_inv \ (D_minus_L * v));
    
    [f_opt_loc, ~] = gmres(afun, y, [], 1e-6, 50);
    
    f_opt_loc = abs(f_opt_loc); % Frequencies are positive
end

function [spec, t, f] = alfmt_forward(x, t, f_vec, f_loc, k)
    % ALFMT_FORWARD computes the transform (Eq. 4)
    % x: signal, t: time, f_vec: frequency vector, f_loc: local frequency
    
    dt = t(2)-t(1);
    N = length(t);
    M = length(f_vec);
    spec = zeros(M, N);
    
    % To vectorize partially, loop over frequencies
    % For each frequency f, the window width depends on time tau via f_loc(tau)
    % h(t, f; tau) = C * exp( - (t-tau)^2 * A^2 / 2k^2 )
    % where A = f_loc(tau) + |f_loc(tau) - f|
    
    for j = 1:M
        freq = f_vec(j);
        
        % Compute width parameter A for all tau at this freq
        % Delta f(tau) assumed to be |f_loc(tau) - freq|
        A_vec = f_loc + abs(f_loc - freq);
        
        % Sigma vector: sigma = k ./ A_vec
        % Window width varies with tau.
        % Convolution is non-stationary.
        
        % Implementation: 
        % V(tau, f) = integral x(t) w(t-tau, tau) e^-i2pft dt
        % Let y(t) = x(t) * e^-i2pft
        % V(tau, f) = integral y(t) w(t-tau, tau) dt
        % This is a summation. 
        % Since w shape depends on tau, we cannot use FFT convolution.
        % We calculate discrete sum. To speed up, limit integration to +/- 3 sigma.
        
        y = x .* exp(-1i * 2 * pi * freq * t(:));
        
        % We iterate over tau (time) to compute the value at each point
        % Optimization: Gaussian window decays fast.
        
        vals = zeros(N, 1);
        for i = 1:N
            tau = t(i);
            A = A_vec(i);
            
            % Determine effective window length (e.g., 3 standard deviations)
            % sigma = k / A;
            % 3*sigma = 3k/A.
            % Indices range
            sigma = k / (A + 1e-6);
            half_width = 3 * sigma; 
            
            % Define range of t indices
            t_min = tau - half_width;
            t_max = tau + half_width;
            
            idx_start = max(1, floor((t_min - t(1))/dt) + 1);
            idx_end = min(N, ceil((t_max - t(1))/dt) + 1);
            idx_range = idx_start:idx_end;
            
            if isempty(idx_range), continue; end
            
            t_sub = t(idx_range).';
            win = (A / (sqrt(2*pi)*k)) * exp( - ((t_sub - tau).^2 * A^2) / (2*k^2) );
            
            vals(i) = sum(y(idx_range) .* win) * dt;
        end
        spec(j, :) = vals;
    end
end

function [x_rec] = alfmt_inverse(spec, f_vec)
    % ALFMT_INVERSE (Eq. 5 / 31)
    % x(t) = Integral ( Integral ALFMT(tau, f) dtau ) e^i2pft df
    % Inner integral is sum over time axis (stacking)
    
    % Sum over tau (time)
    % Note: This assumes the window partition of unity property roughly holds
    % or that we reconstruct from the "voice" at each frequency.
    % Based on S-transform inversion: Integral S(tau, f) e^i2pft df ? No.
    % Standard Inverse: sum over frequency of S(t, f) e^i2pft?
    % Eq 5 says: Integral dtau Integral df S(tau, f) e^i2pft
    % Let X_f(f) = Integral S(tau, f) dtau.
    % Then x(t) = Integral X_f(f) e^i2pft df = Inverse Fourier of X_f.
    
    % 1. Collapse time axis
    X_f = sum(spec, 2); % Sum over columns (time)
    
    % 2. Inverse FFT?
    % Since f_vec might not be uniform or satisfy FFT grid, we use discrete sum.
    % x(t) = sum X_f(j) * exp(i 2 pi f_j t) df
    
    df = f_vec(2) - f_vec(1);
    [M, N] = size(spec);
    
    % Create t vector if not passed (assuming normalized 0..N-1 or passed in context)
    % We compute for N samples.
    % Matrix multiplication: x = exp(i 2 pi f t)' * (X_f * df)
    % Wait, Eq 5 implies x(t).
    
    % Using IFFT approximation if f starts at 0 and is linear
    % Construct full spectrum for IFFT
    % Here we just use direct summation for clarity
    t_idx = 0:N-1; 
    % Assumes t starts at 0, dt=1 for reconstruction grid, scaled later.
    % Actually we need to match original phase.
    % Let's use the explicit sum formula.
    
    % Note: X_f includes the dt factor from the forward integral sum? 
    % In forward: sum * dt. So X_f is approx integral.
    
    % Sum over frequencies
    % x(n) = sum_j X_f(j) * exp(i 2 pi f_j t_n) * df
    % Note: Inverse ST usually doesn't need sum over tau?
    % Standard ST Inverse: x(t) = integral S(t, f) e^i2pft df? No.
    % Standard ST Inverse: x(t) = integral (integral S(tau, f) dtau) e^i2pft df.
    % Yes, Eq 5 matches standard S-transform inverse.
    
    X_f = sum(spec, 2); % This is sum, need * dt_sampling? 
    % Forward was sum * dt. So X_f is roughly Integral(S dtau).
    
    % Reconstruct
    % We need to multiply by exp(i 2 pi f t). 
    % Since we don't have t here, we assume standard reconstruction length.
    
    % Use simple inverse Fourier Logic if uniform
    % x = ifft(X_f) ?
    % Only if X_f represents the Fourier spectrum.
    % alfmt_forward gives Time-Frequency distribution.
    % Summing over time gives Fourier Spectrum.
    % So yes, IFFT of X_f recovers x.
    
    % However, f_vec is 0 to Fmax. We need two-sided or real signal assumption.
    % We assume real signal x.
    x_rec_complex = zeros(1, N);
    t_vec = (0:N-1) * (1/(2*max(f_vec))); % Placeholder t
    % We need original t. But let's do discrete sum:
    
    % To make it robust:
    % x(n) = Real ( 2 * Sum_{f>0} X_f(f) e^{i 2 pi f t_n} df ) + X_f(0)
    % Assuming sampling dt from context (passed implicitly or normalized).
    % In this demo, we do loop reconstruction.
    
    x_rec = zeros(1, N);
    % Assuming t goes 0 to (N-1)*dt
    % We need dt. But let's assume f_vec is in Hz and t in s.
    % We will pass t in context or assume dt from f_vec?
    % Usually df * dt * N = 1.
    
    % Simplified reconstruction:
    % x_rec = real(ifft(X_f_full))
    % But let's stick to the summation integral Eq 5.
    
    % For reproduction, we use the `alfmt_dynamic_deconvolution` loop which handles recovery.
    x_rec = 0; % Placeholder
end

function [r_est, seis_rec, spec, spec_rec] = alfmt_dynamic_deconvolution(seis, t, f_loc, params)
    dt = t(2)-t(1);
    f_vec = params.freq_vec;
    k = params.k;
    
    % 1. Forward ALFMT
    spec = alfmt_forward(seis, t, f_vec, f_loc, k);
    
    % 2. Estimate Wavelet Spectrum (Hyperbolic/Smooth Smoothing)
    % Approximation: Smooth the magnitude spectrum
    % We apply a 2D Gaussian filter to |S(tau, f)|
    % Window size: 200ms x 10Hz approx
    sigma_t = 0.1 / dt; % 100ms
    sigma_f = 5 / (f_vec(2)-f_vec(1)); % 5Hz
    
    amp_spec = abs(spec);
    % Use imgaussfilt for smoothing (requires Image Proc Toolbox)
    % Fallback if unavailable: simple convolution
    try
        W_amp = imgaussfilt(amp_spec, [sigma_f, sigma_t]);
    catch
        W_amp = smooth2d(amp_spec, sigma_f, sigma_t);
    end
    
    % 3. Estimate Phase (Minimum Phase Assumption) Eq. 29
    % Phi(f) = Hilbert( ln(W_amp) ) over frequency axis
    ln_A = log(W_amp + 1e-6);
    Phase = zeros(size(ln_A));
    for i = 1:length(t)
        Phase(:, i) = imag(hilbert(ln_A(:, i)));
    end
    
    % 4. Estimate Reflectivity Spectrum (Eq. 30)
    mu = params.mu;
    A_max = max(W_amp(:));
    
    % Deconvolution in TF domain
    % R_est = S / ( |W| + mu*Amax ) * exp(-i phi)
    % Note: S has phase of W + phase of R. We remove W phase.
    % Wait, S is complex. S = |S| exp(i phi_total).
    % We want R. R ~ S / W.
    % W_est = |W|_smooth * exp(i Phase).
    % R_est = S / W_est.
    
    W_est = W_amp .* exp(1i * Phase);
    % Stabilized division
    Denominator = W_est;
    % Using the paper's formula style (magnitude stabilization)
    % R = S * conj(W) / (|W|^2 + mu) ?
    % Paper Eq 30: V_r = V_s / (|V_s|_hs + mu*Amax) * exp(-i phi)
    % This is effectively S / (|W| + eps) * exp(-i phi_w)
    % = |S|/|W| * exp(i(phi_s - phi_w)) = |R| exp(i phi_r)
    
    spec_rec = spec ./ (W_amp + mu * A_max) .* exp(-1i * Phase);
    
    % 5. Inverse ALFMT to get reflectivity r(t)
    % Reconstruct r(t)
    % Sum over tau (collapse to Fourier)
    R_f = sum(spec_rec, 2) * dt; % Scale by dt for integration
    
    % Inverse Fourier Transform
    % We construct a 2-sided spectrum for IFFT or sum cosines
    % x(t) = 2 * Real( Sum R_f * exp(i 2 pi f t) df )
    df = f_vec(2) - f_vec(1);
    
    % Vectorized Inverse Summation
    % r_est(t)
    E = exp(1i * 2 * pi * f_vec(:) * t); % M x N matrix
    r_est = 2 * real(sum(R_f .* E, 1) * df); % Factor 2 for one-sided
    
    % 6. Convolve with Desired Wavelet (e.g. Ricker 30Hz) Eq. 32
    % Generate target wavelet
    [wav, tw] = ricker_wavelet(30, dt);
    seis_rec = conv(r_est, wav, 'same');
    
end

function [w, t] = ricker_wavelet(f, dt)
    nw = 2/f/dt;
    t = -floor(nw/2)*dt : dt : floor(nw/2)*dt;
    s = (1 - 2*pi^2*f^2*t.^2) .* exp(-pi^2*f^2*t.^2);
    w = s;
end

function out = smooth2d(img, sig_r, sig_c)
    % Basic 2D smoothing using separable convolution
    gw_r = gausswin(round(6*sig_r)+1); gw_r = gw_r/sum(gw_r);
    gw_c = gausswin(round(6*sig_c)+1); gw_c = gw_c/sum(gw_c);
    out = conv2(gw_r, gw_c, img, 'same');
end

function seis = generate_nonstationary_seismic(r, t, f_dom, Q)
    % Generates seismic trace with Q attenuation
    % S(t) = Sum r(tau) * w(t-tau, tau)
    % w(t, tau) is evolved wavelet
    
    seis = zeros(size(t));
    dt = t(2)-t(1);
    
    % Base Ricker spectrum
    Nt = length(t);
    Nf = 2^nextpow2(Nt);
    f = (0:Nf/2)' * (1/dt/Nf);
    
    % Ricker freq domain
    om = 2*pi*f;
    W0 = (2/sqrt(pi)) * (f.^2 / f_dom^2) .* exp(-(f.^2)/f_dom^2);
    % Proper scaling/phase for Ricker
    % Typically Ricker in time: (1-2x^2)exp(-x^2). Spectrum is above.
    
    % Constant Q Attenuation Filter
    % A(f, tau) = exp(- pi f tau / Q)
    % Phi(f, tau) = Minimum phase via Hilbert
    
    for i = 1:length(t)
        if r(i) == 0, continue; end
        
        tau = t(i);
        
        % Attenuation
        Atten = exp(-pi * f * tau / Q);
        
        % Minimum Phase
        % H(ln A). A is real exp. ln A = -pi f tau / Q.
        % Hilbert of linear freq (-k*f)
        % For constant Q, use standard phase dispersion term:
        % exp( -i * 2 * pi * f * tau * (1/pi/Q) * ln(f/f_ref) )
        % We use f_dom as reference
        f_ref = f_dom;
        Disp = exp( -1i * 2 * pi * f * tau .* (1/pi/Q) .* log( (f+1e-6)/f_ref ) );
        
        % Total Wavelet Spectrum at tau
        W_tau = W0 .* Atten .* Disp;
        
        % Time domain wavelet
        w_t = real(ifft(W_tau, Nf));
        w_t = [w_t(Nf/2+1:end); w_t(1:Nf/2)]; % fftshift
        w_t = w_t * (2*Nf); % Scale fix approximation
        
        % Shift to tau
        % We generated w centered at 0 (after shift).
        % Add to seismic at position i
        % Simple valid convolution accumulation
        
        % Cut wavelet to reasonable length
        Lw = round(0.15/dt); % 150ms half width
        center = Nf/2;
        w_valid = w_t(center-Lw : center+Lw);
        
        i_start = i - Lw;
        i_end = i + Lw;
        w_start = 1; 
        w_end = length(w_valid);
        
        % Boundary checks
        if i_start < 1
            w_start = 1 + (1 - i_start);
            i_start = 1;
        end
        if i_end > Nt
            w_end = w_end - (i_end - Nt);
            i_end = Nt;
        end
        
        if w_start <= w_end
            seis(i_start:i_end) = seis(i_start:i_end) + r(i) * w_valid(w_start:w_end).';
        end
    end
end