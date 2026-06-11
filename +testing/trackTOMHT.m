function results = trackTOMHT(x, y, t, opts)
%TRACKTOMHT  Track-Oriented Multiple Hypothesis Tracking for a 2-D point cloud.
%
%   results = TRACKTOMHT(x, y, t) processes a point cloud given as vectors
%   x, y, t.  Each triple (x(i), y(i), t(i)) is one measurement: a 2-D
%   position observed at time t(i).  Measurements sharing the same time value
%   form one "scan" (slice).  Within a scan, some points lie on a smooth,
%   slowly-curving trajectory ("real") and the rest are clutter ("noise").
%   Detections may be missing on any scan -- including many scans in a row,
%   after which the target may reappear.
%
%   results = TRACKTOMHT(x, y, t, opts) overrides defaults via struct OPTS.
%
%   TRACKTOMHT() with no arguments runs a built-in synthetic demo (two
%   curving tracks, clutter, and a deliberate long dropout) and plots the
%   result -- a quick way to confirm the function behaves as expected.
%
%   This is a self-contained implementation of Track-Oriented MHT.  It uses
%   NO MATLAB toolboxes.  Core ingredients:
%       * a constant-velocity (CV) Kalman filter per track; the bounded
%         process noise q absorbs the "slight curvature" of the paths
%         (Bar-Shalom, Li & Kirubarajan, "Estimation with Applications...",
%         2001);
%       * track trees that BRANCH on every gated measurement plus a
%         missed-detection hypothesis (Reid, IEEE TAC 24(6), 1979;
%         Kurien, 1990);
%       * a log-likelihood-ratio (LLR) track score (Sittler, IEEE TMC 1964;
%         Blackman, IEEE AES Mag 19(1), 2004; Bar-Shalom, Willett & Tian,
%         "Tracking and Data Fusion", 2011, ch. on track scoring);
%       * N-scan pruning -- the deferred-decision mechanism that lets a
%         track coast through a long dropout and re-acquire the target
%         afterwards (Blackman 2004);
%       * a global hypothesis taken as a maximum-weight set of mutually
%         compatible tracks (no two tracks share a measurement).  The exact
%         problem is a maximum-weight independent set (NP-hard) that full
%         TOMHT solves with an integer program; here it is approximated
%         greedily, which is the standard practical fallback.
%
%   OUTPUT struct RESULTS:
%       .tracks    confirmed tracks (struct array), each with fields
%                  .est     2 x nScan estimated [x;y] (NaN before birth)
%                  .detFlag 1 x nScan, true where a detection was used,
%                           false where the scan was coasted (a miss)
%                  .hist    1 x nScan global measurement index, or NaN
%                  .score   final LLR track score
%                  .len     number of detections in the track
%       .scanTime  1 x nScan time stamp of each scan
%       .diag      per-scan diagnostics: .nHyp (#live hypotheses),
%                  .bestScore (best LLR)
%
%   Example:
%       trackTOMHT;                       % run the demo
%       r = trackTOMHT(x, y, t);          % run on your own data
%       o.Pd = 0.85; o.maxCoast = 12;     % tune, then
%       r = trackTOMHT(x, y, t, o);

% ===================== 0. Demo mode =======================================
if nargin == 0
    [x, y, t, truth] = local_demo_data();
    opts = struct();
    opts.truth = truth;                  % pass ground truth for overlay plots
end
if nargin < 4 || isempty(opts), opts = struct(); end

x = x(:); y = y(:); t = t(:);
N = numel(x);

% ===================== 1. Options / defaults ==============================
Pd         = optget(opts, 'Pd',         0.1);  % detection probability
sigma      = optget(opts, 'sigma',      5.80);  % meas. noise std (per axis)
q          = optget(opts, 'q',          0.03);  % CV process-noise PSD (accel)
gate       = optget(opts, 'gate',      11.83);  % chi-sq gate, 2 DOF (~0.997).
                                                % Hardcoded to avoid chi2inv.
Nscan      = optget(opts, 'Nscan',      5);     % N-scan pruning depth
maxCoast   = optget(opts, 'maxCoast',  12);     % max consecutive misses
maxBranch  = optget(opts, 'maxBranch', 10);     % best leaves kept per track tree
maxTracks  = optget(opts, 'maxTracks', 300);    % cap on live hypotheses
confirm    = optget(opts, 'confirmThresh', 8);  % LLR needed to confirm
minLen     = optget(opts, 'minLen',     6);     % min detections to confirm
v0var      = optget(opts, 'v0var',      4);     % init velocity variance
newScore0  = optget(opts, 'newScore0',  0);     % LLR of a brand-new track
dtBin      = optget(opts, 'dtBin',     []);     % optional time-bin width
doPlot     = optget(opts, 'plot',    true);

% ===================== 2. Group measurements into scans (slices) ==========
% Points are grouped by equal time value.  For (near-)continuous t (e.g. an
% event stream) pass opts.dtBin to bin time into slices first.
if isempty(dtBin)
    key = t;
else
    key = round(t / dtBin);
end
[~, ~, scanIdx] = unique(key);          % scanIdx: each point's scan, sorted
nScan = max(scanIdx);
scanTime    = zeros(1, nScan);
scanMembers = cell(1, nScan);
for k = 1:nScan
    m = find(scanIdx == k);
    scanMembers{k} = m;
    scanTime(k) = mean(t(m));
end

% ===================== 3. Static model pieces =============================
H  = [1 0 0 0; 0 0 1 0];                % observe position only
R  = sigma^2 * eye(2);                  % measurement covariance
I4 = eye(4);

% Clutter spatial density beta_FT (false alarms / unit area / scan), used in
% the LLR.  Estimated from the data unless supplied.
bb     = [min(x) max(x) min(y) max(y)];
area   = max(bb(2)-bb(1), 1) * max(bb(4)-bb(3), 1);
betaFT = optget(opts, 'betaFT', N / (nScan * area));
logBeta = log(max(betaFT, realmin));

% ===================== 4. TOMHT recursion =================================
tracks   = [];                          % struct array of live track leaves
nextTree = 1;                           % unique id per track tree (root)
diagNHyp = zeros(1, nScan);
diagBest = -inf(1, nScan);

for k = 1:nScan
    mIdx = scanMembers{k};              % global indices of this scan's points
    Z    = [x(mIdx) y(mIdx)]';          % 2 x m measurements
    m    = numel(mIdx);

    if k == 1
        dt = 1;                         % no prediction on first scan
    else
        dt = scanTime(k) - scanTime(k-1);
        if dt <= 0, dt = 1; end
    end
    [F, Q] = local_cvmodel(dt, q);

    newTracks      = {};                % branches produced this scan
    gatedByExisting = false(1, m);      % was meas j claimed by any track?

    % ---- 4a. extend every existing track (branching) --------------------
    for it = 1:numel(tracks)
        tr = tracks(it);

        xp = F * tr.state;              % predicted state
        Pp = F * tr.cov * F' + Q;       % predicted covariance
        S  = H * Pp * H' + R;           % innovation covariance
        Sinv    = inv(S);
        logdetS = log(max(det(S), realmin));

        % --- missed-detection branch (the target was not seen this scan) --
        mb         = tr;
        mb.state   = xp;
        mb.cov     = Pp;
        mb.score   = tr.score + log(1 - Pd);          % miss penalty
        mb.hist(k) = NaN;
        mb.est(:,k)= H * xp;
        mb.detFlag(k) = false;
        mb.missRun = tr.missRun + 1;
        if mb.missRun <= maxCoast       % delete after too many misses
            newTracks{end+1} = mb;
        end

        % --- one detection branch per gated measurement -------------------
        for j = 1:m
            nu = Z(:,j) - H * xp;       % innovation
            d2 = nu' * Sinv * nu;       % squared Mahalanobis distance
            if d2 <= gate
                gatedByExisting(j) = true;
                Kg = Pp * H' * Sinv;    % Kalman gain
                xu = xp + Kg * nu;
                Pu = (I4 - Kg * H) * Pp;

                db          = tr;
                db.state    = xu;
                db.cov      = Pu;
                % LLR increment for a detection (target vs. clutter):
                %   dL = ln(Pd) + ln N(z; Hxp,S) - ln(beta_FT)
                db.score    = tr.score + log(Pd) - logBeta ...
                              - 0.5*(2*log(2*pi) + logdetS + d2);
                db.hist(k)  = mIdx(j);
                db.est(:,k) = H * xu;
                db.detFlag(k) = true;
                db.missRun  = 0;
                db.len      = tr.len + 1;
                newTracks{end+1} = db;
            end
        end
    end

    % ---- 4b. seed a new track tree from each unclaimed measurement -------
    % (Only measurements outside every existing gate start new trees; this
    %  keeps the hypothesis count bounded and lets genuinely new / re-born
    %  targets appear.)
    for j = 1:m
        if ~gatedByExisting(j)
            s          = local_make_track(nScan);
            s.state    = [Z(1,j); 0; Z(2,j); 0];      % one-point init, v=0
            s.cov      = diag([sigma^2, v0var, sigma^2, v0var]);
            s.score    = newScore0;
            s.hist(k)  = mIdx(j);
            s.est(:,k) = Z(:,j);
            s.detFlag(k) = true;
            s.treeId   = nextTree;  nextTree = nextTree + 1;
            s.startScan = k;
            s.missRun  = 0;
            s.len      = 1;
            newTracks{end+1} = s;
        end
    end

    if isempty(newTracks)
        tracks = [];
        diagNHyp(k) = 0;
        continue
    end
    tracks = [newTracks{:}];            % cell -> struct array

    % ---- 4c. delete un-established seeds, then prune within each tree -----
    % A seed measurement that never gathered a 2nd detection (still len==1
    % after a few scans) is treated as clutter and dropped.  Pruning is done
    % PER TREE (keep the maxBranch best leaves of each tree) so that a weaker
    % target -- e.g. one penalised by a dropout -- is never deleted just
    % because some other target's track scores higher.
    estd   = ~([tracks.len] == 1 & [tracks.missRun] >= 3);
    tracks = tracks(estd);
    tracks = local_pertree_prune(tracks, maxBranch);

    % ---- 4d. N-scan pruning (commit decisions N scans back) --------------
    tracks = local_nscan_prune(tracks, k, Nscan);

    % ---- 4e. hard cap on number of live hypotheses -----------------------
    if numel(tracks) > maxTracks
        sc = [tracks.score];
        [~, ord] = sort(sc, 'descend');
        tracks = tracks(ord(1:maxTracks));
    end

    diagNHyp(k) = numel(tracks);
    diagBest(k) = max([tracks.score]);
end

% ===================== 5. Global hypothesis + confirmation ================
% Greedy maximum-weight set of mutually compatible tracks: take tracks in
% descending score and accept one only if it shares no measurement with an
% already-accepted track.
results.scanTime = scanTime;
results.diag     = struct('nHyp', diagNHyp, 'bestScore', diagBest);

empty = local_make_track(nScan); empty(1) = [];   % 0x0 typed struct array
if isempty(tracks)
    results.tracks = empty;
    warning('trackTOMHT:noTracks', 'No tracks survived.');
else
    sc = [tracks.score];
    [~, ord] = sort(sc, 'descend');
    used = false(1, N);
    sel  = [];
    for ii = ord
        idxs = tracks(ii).hist(~isnan(tracks(ii).hist));
        if isempty(idxs),        continue, end
        if any(used(idxs)),      continue, end       % conflict
        used(idxs) = true;
        sel(end+1) = ii;
    end
    chosen = tracks(sel);
    keep   = ([chosen.score] >= confirm) & ([chosen.len] >= minLen);
    if any(keep)
        results.tracks = chosen(keep);
    else
        results.tracks = empty;
    end
end

% ===================== 6. Plots ===========================================
if doPlot
    local_plots(results, x, y, t, opts);
end
end % ====================== end main function ==============================


% -------------------------------------------------------------------------
function [F, Q] = local_cvmodel(dt, q)
%LOCAL_CVMODEL  Constant-velocity transition F and discrete white-noise-
%   acceleration process covariance Q for state [px; vx; py; vy].
F  = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
Qx = q * [dt^4/4, dt^3/2; dt^3/2, dt^2];
Q  = [Qx, zeros(2); zeros(2), Qx];
end


% -------------------------------------------------------------------------
function s = local_make_track(nScan)
%LOCAL_MAKE_TRACK  Allocate a track struct with a fixed field order so that
%   struct arrays concatenate cleanly.
s.state    = zeros(4,1);
s.cov      = eye(4);
s.score    = 0;
s.hist     = nan(1, nScan);
s.est      = nan(2, nScan);
s.detFlag  = false(1, nScan);
s.treeId   = 0;
s.startScan = 1;
s.missRun  = 0;
s.len      = 0;
end


% -------------------------------------------------------------------------
function tracks = local_nscan_prune(tracks, k, N)
%LOCAL_NSCAN_PRUNE  Within each track tree, keep only branches that agree
%   with the best-scoring branch at scan (k-N+1).  This realises the MHT
%   "deferred decision": the assignment N scans back is committed to the
%   currently most likely branch, while recent scans stay ambiguous.
col = k - N + 1;
if col < 1, return, end
trees = [tracks.treeId];
ut    = unique(trees);
keepMask = true(1, numel(tracks));
for i = 1:numel(ut)
    idx = find(trees == ut(i));
    if numel(idx) <= 1, continue, end
    [~, b] = max([tracks(idx).score]);
    bv = tracks(idx(b)).hist(col);                  % best branch's assignment
    for jj = 1:numel(idx)
        hv = tracks(idx(jj)).hist(col);
        agree = (isnan(hv) && isnan(bv)) || (hv == bv);
        if ~agree, keepMask(idx(jj)) = false; end
    end
end
tracks = tracks(keepMask);
end


% -------------------------------------------------------------------------
function tracks = local_pertree_prune(tracks, K)
%LOCAL_PERTREE_PRUNE  Keep the K best-scoring leaves within each track tree.
if isempty(tracks), return, end
trees = [tracks.treeId];
ut    = unique(trees);
keep  = false(1, numel(tracks));
for i = 1:numel(ut)
    idx = find(trees == ut(i));
    if numel(idx) <= K
        keep(idx) = true;
    else
        [~, ord] = sort([tracks(idx).score], 'descend');
        keep(idx(ord(1:K))) = true;
    end
end
tracks = tracks(keep);
end


% -------------------------------------------------------------------------
function v = optget(opts, f, d)
%OPTGET  Return opts.(f) if present & non-empty, else default d.
if isfield(opts, f) && ~isempty(opts.(f)), v = opts.(f); else, v = d; end
end


% -------------------------------------------------------------------------
function local_plots(results, x, y, t, opts)
%LOCAL_PLOTS  Three validation figures.
trk = results.tracks;
nT  = numel(trk);
cols = lines(max(nT,1));
hasTruth = isfield(opts,'truth') && ~isempty(opts.truth);

% --- Fig 1: top-down x-y view ---
figure('Name','TOMHT - XY');
hold on; box on; grid on;
scatter(x, y, 12, [0.72 0.72 0.72], 'filled');      % all points (gray)
if hasTruth
    for c = 1:numel(opts.truth.paths)
        p = opts.truth.paths{c};
        plot(p(:,1), p(:,2), 'k--', 'LineWidth', 1);
    end
end
for c = 1:nT
    e  = trk(c).est;
    al = ~all(isnan(e), 1);                         % scans where track alive
    plot(e(1,al), e(2,al), '-', 'Color', cols(c,:), 'LineWidth', 1.6);
    df = al & trk(c).detFlag;                        % detections
    co = al & ~trk(c).detFlag;                       % coasted (misses)
    plot(e(1,df), e(2,df), 'o', 'Color', cols(c,:), ...
         'MarkerFaceColor', cols(c,:), 'MarkerSize', 4);
    plot(e(1,co), e(2,co), 'x', 'Color', cols(c,:), ...
         'MarkerSize', 8, 'LineWidth', 1.3);
end
xlabel('x'); ylabel('y'); axis equal;
title(sprintf('TOMHT: %d confirmed track(s).  o = detection,  x = coasted,  -- = truth', nT));

% --- Fig 2: x-y-t (shows the slices and continuity through dropouts) ---
figure('Name','TOMHT - XYT');
hold on; box on; grid on;
scatter3(x, y, t, 10, [0.72 0.72 0.72], 'filled');
for c = 1:nT
    e  = trk(c).est;
    al = ~all(isnan(e), 1);
    plot3(e(1,al), e(2,al), results.scanTime(al), '-', ...
          'Color', cols(c,:), 'LineWidth', 1.8);
end
xlabel('x'); ylabel('y'); zlabel('t'); view(40, 25);
title('Point cloud (gray) and confirmed tracks through time');

% --- Fig 3: diagnostics ---
figure('Name','TOMHT - diagnostics');
subplot(2,1,1);
plot(1:numel(results.diag.nHyp), results.diag.nHyp, '-o', 'MarkerSize', 3);
ylabel('# live hypotheses'); grid on; box on;
title('Track-tree branching vs. pruning');
subplot(2,1,2);
bs = results.diag.bestScore; bs(~isfinite(bs)) = NaN;
plot(1:numel(bs), bs, '-', 'LineWidth', 1.2);
xlabel('scan'); ylabel('best LLR score'); grid on; box on;
title('Best track score (log-likelihood ratio)');
end


% -------------------------------------------------------------------------
function [x, y, t, truth] = local_demo_data()
%LOCAL_DEMO_DATA  Two smoothly-curving tracks + clutter + a long dropout.
rng(1);
nScan = 60;
k  = (1:nScan)';
Ax = 5 + 1.25*k;            Ay = 28 + 10*sin(2*pi*k/110);   % gentle curve
Bx = 6 + 1.05*k;            By = 60 - 0.55*k + 0.011*k.^2;  % gentle parabola
A  = [Ax Ay];   B = [Bx By];

PdTrue = 0.90;  sig = 0.80;  lamC = 2.0;
gapA = (k >= 27 & k <= 34);                     % 8-scan dropout for track A

x = []; y = []; t = [];
for kk = 1:nScan
    if rand < PdTrue && ~gapA(kk)               % track A detection
        x(end+1) = A(kk,1) + sig*randn;  y(end+1) = A(kk,2) + sig*randn;  t(end+1) = kk;
    end
    if rand < PdTrue                            % track B detection
        x(end+1) = B(kk,1) + sig*randn;  y(end+1) = B(kk,2) + sig*randn;  t(end+1) = kk;
    end
    nc = local_poisson(lamC);                   % clutter
    for c = 1:nc
        x(end+1) = 2 + 85*rand;  y(end+1) = 10 + 60*rand;  t(end+1) = kk;
    end
end
x = x(:); y = y(:); t = t(:);
truth.paths = {A, B};
truth.gapA  = find(gapA);
end


% -------------------------------------------------------------------------
function n = local_poisson(lam)
%LOCAL_POISSON  Knuth's Poisson sampler (no Statistics Toolbox needed).
L = exp(-lam);  k = 0;  p = 1;
while true
    k = k + 1;  p = p * rand;
    if p <= L, break, end
end
n = k - 1;
end