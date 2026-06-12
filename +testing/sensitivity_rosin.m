% SENSITIVITY_ROSIN  Probe rosinThreshold stability and failure modes.

rand('state', 7); randn('state', 7);
H = 480; W = 640;

function score_map = make_map(noise_mean, sig_lo, sig_hi, n_noise, n_sig)
    H = 480; W = 640;
    score_map = zeros(H, W);
    idx = randperm(H * W, n_noise);
    score_map(idx) = -noise_mean * log(rand(1, n_noise));
    free = find(score_map == 0);
    sidx = free(randperm(numel(free), n_sig));
    score_map(sidx) = sig_lo + (sig_hi - sig_lo) * rand(1, n_sig);
end

base = make_map(0.004, 0.05, 0.55, 135000, 7000);

fprintf('--- Binning / smoothing sensitivity (same map) ---\n');
for nb = [64 128 256 512 1024]
    th = rosinThreshold(base, 'NumBins', nb);
    fprintf('  NumBins = %4d                  -> th = %.5f\n', nb, th);
end
for sw = [1 3 5 11 21]
    th = rosinThreshold(base, 'SmoothWidth', sw);
    fprintf('  SmoothWidth = %2d               -> th = %.5f\n', sw, th);
end

fprintf('--- Signal/noise overlap (signal lower edge moves down) ---\n');
for lo = [0.05 0.03 0.02 0.01 0.005]
    m  = make_map(0.004, lo, 0.55, 135000, 7000);
    th = rosinThreshold(m);
    n_idx = m > 0 & m < lo;  % approx: pixels certainly noise
    fprintf('  signal starts at %.3f          -> th = %.5f\n', lo, th);
end

fprintf('--- Signal fraction (7000 -> 700 -> 70 signal pixels) ---\n');
for ns = [7000 2000 700 200 70]
    m  = make_map(0.004, 0.05, 0.55, 135000, ns);
    th = rosinThreshold(m);
    fprintf('  n_signal = %5d               -> th = %.5f\n', ns, th);
end

fprintf('--- Pathological: noise only, NO signal ---\n');
m_pure = make_map(0.004, 0.5, 0.55, 135000, 0 + 1);  % 1 token pixel
m_pure(m_pure >= 0.5) = 0;                            % strip it: pure noise
th = rosinThreshold(m_pure);
frac_removed = sum(m_pure(:) > 0 & m_pure(:) < th) / sum(m_pure(:) > 0);
fprintf('  pure exponential noise         -> th = %.5f (removes %.1f%%)\n', ...
    th, 100 * frac_removed);
