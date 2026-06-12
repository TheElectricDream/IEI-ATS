% DEMO_ROSINTHRESHOLD  Synthetic behavior check for rosinThreshold.
%
% Builds a 640x480 score map mimicking the observed filter_mask
% distribution: ~135k noise pixels with exponential scores driven
% toward zero, plus ~7k signal pixels spread over [0.05, 0.7]
% (clustered spatially, though spatial structure is irrelevant to
% a histogram method), plus a handful of extreme stragglers.

rng_seed = 42;
rand('state', rng_seed); randn('state', rng_seed);

H = 480; W = 640;
score_map = zeros(H, W);

% --- Noise population: exponential, mean 0.004 (the spike) ---------
n_noise = 135000;
noise_idx = randperm(H * W, n_noise);
score_map(noise_idx) = -0.004 * log(rand(1, n_noise));

% --- Signal population: broad low tail, ~7k pixels -----------------
n_sig = 7000;
free = find(score_map == 0);
sig_idx = free(randperm(numel(free), n_sig));
% Mixture: most signal between 0.05 and 0.3, some up to 0.55
sig_scores = 0.05 + 0.25 * rand(1, n_sig);
hi = rand(1, n_sig) < 0.08;
sig_scores(hi) = 0.3 + 0.25 * rand(1, sum(hi));
score_map(sig_idx) = sig_scores;

% --- A few extreme stragglers (surviving hot pixels) ---------------
strag = free(end-4:end);
score_map(strag) = 0.65 + 0.05 * rand(1, 5);

% --- Run: Rosin original anchor and robust-quantile anchor ---------
[th_orig, d1] = rosinThreshold(score_map);
[th_q999, d2] = rosinThreshold(score_map, 'TailQuantile', 0.999);
[th_fine, d3] = rosinThreshold(score_map, 'NumBins', 512);

fprintf('Threshold (Rosin original anchor):      %.5f\n', th_orig);
fprintf('Threshold (TailQuantile = 0.999):       %.5f\n', th_q999);
fprintf('Threshold (512 bins, original anchor):  %.5f\n', th_fine);

% Ground-truth-aware sanity numbers
noise_scores = score_map(noise_idx);
fprintf('Noise pixels above th_orig:  %d / %d (%.3f%%)\n', ...
    sum(noise_scores > th_orig), n_noise, 100*sum(noise_scores > th_orig)/n_noise);
fprintf('Signal pixels below th_orig: %d / %d (%.3f%%)\n', ...
    sum(sig_scores < th_orig), n_sig, 100*sum(sig_scores < th_orig)/n_sig);

% --- Diagnostic plot ------------------------------------------------
figure('visible', 'off', 'position', [0 0 1400 500]);

subplot(1, 2, 1);
bar(d1.bin_centers, d1.counts, 'histc'); hold on;
plot(d1.chord_x, d1.chord_y, 'r-', 'linewidth', 2);
plot([th_orig th_orig], ylim, 'k--', 'linewidth', 2);
plot([th_q999 th_q999], ylim, 'm:', 'linewidth', 2);
title('Rosin threshold, linear count axis');
xlabel('Score'); ylabel('Count');
legend('Histogram', 'Chord', 'th (original)', 'th (q=0.999)');

subplot(1, 2, 2);
semilogy(d1.bin_centers, d1.counts_smooth + 1, 'b-'); hold on;
plot([th_orig th_orig], ylim, 'k--', 'linewidth', 2);
plot([th_q999 th_q999], ylim, 'm:', 'linewidth', 2);
title('Smoothed histogram, log count axis');
xlabel('Score'); ylabel('Count + 1');

print('rosin_demo.png', '-dpng', '-r110');
