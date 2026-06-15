function [iei, mean_map, var_map, valid_map] = updateWindowedIEI( ...
        iei, ev_x, ev_y, ev_t, imgSz, K)
% COMPUTEWINDOWEDIEI  Per-pixel inter-event-interval (IEI) mean and variance
% for one time window, backfilled from past windows where the window is
% too sparse to estimate reliably.
%
% For every pixel active in the window:
%   1. count its events this window (c);
%   2. if c < K, borrow the newest (K - c) real timestamps retained from
%      previous windows so the estimate is formed from up to K points;
%   3. compute the IEI mean and Bessel-corrected variance over that
%      timestamp list, using the same closed-form identities as the
%      original computeNeighborhoodStats.
%
% Borrowed timestamps are used ONLY in the statistics. They are retained
% solely because they are among the newest K real timestamps at that pixel;
% they are never emitted as events.
%
%   K = 1  -> no pixel is ever sparse (every active pixel has >= 1 event),
%             so no backfill occurs and mean_map/var_map are exactly the IEI
%             maps of the window alone. Sanity check: mean_map equals the
%             old t_mean_diff and var_map equals t_std_diff.^2.
%   K >= 2 -> sparse pixels are topped up to K points.
%
% Persistent state (ONLY the ring buffer survives between calls):
%   iei.buf  [nPix x K]  newest timestamp in the last column, NaN = empty.
%
% Outputs are FRESH every window (zero / false at inactive pixels), matching
% the convention of computeNeighborhoodStats. Any temporal holding or
% smoothing of the per-pixel tau belongs in a downstream EMA, not here.
%
% Coordinates: ev_x = row, ev_y = col. Inputs need not be pre-sorted.

    nPix = prod(imgSz);
    if ~isfield(iei,'buf') || isempty(iei.buf) || size(iei.buf,2) ~= K
        iei.buf = nan(nPix, K);
    end

    mean_map  = zeros(imgSz);
    var_map   = zeros(imgSz);
    valid_map = false(imgSz);
    if isempty(ev_t), return; end

    ev_x = ev_x(:);  ev_y = ev_y(:);  ev_t = ev_t(:);

    % ----------------------------------------------------------------------
    % 1. Index the window's events by active pixel.
    % ----------------------------------------------------------------------
    lin = sub2ind(imgSz, ev_x, ev_y);
    uid = unique(lin, 'stable');             % active pixels (linear index)
    A   = numel(uid);
    row_of_pix      = zeros(nPix, 1);
    row_of_pix(uid) = (1:A)';
    win_row = row_of_pix(lin);               % active-pixel row of each event
    counts  = accumarray(win_row, 1, [A 1]); % events this window per pixel

    % ----------------------------------------------------------------------
    % 2. Pull backfill from the ring buffer: newest (K - counts) per pixel.
    %    need = 0 for any pixel with counts >= K, so dense pixels borrow none.
    % ----------------------------------------------------------------------
    need = max(0, K - counts);

    oldbuf  = iei.buf(uid, :);               % [A x K], newest in last column
    [br, bc] = find(~isnan(oldbuf));         % rows/cols of retained timestamps
    br = br(:);  bc = bc(:);
    old_ts  = oldbuf(sub2ind([A K], br, bc));
    old_ts  = old_ts(:);                     % row when A==1; force column

    % buffer is right-packed, so newest entry is the largest column:
    rank_from_new = K - bc + 1;              % 1 = newest
    borrow = rank_from_new <= need(br);      % take the newest 'need' of each
    bx_borrow = br(borrow);
    bt_borrow = old_ts(borrow);

    % ----------------------------------------------------------------------
    % 3. Statistics stream = window events + borrowed timestamps, grouped by
    %    pixel and sorted in time (borrowed are older, so they fall first).
    %    Then apply the computeNeighborhoodStats IEI identities.
    % ----------------------------------------------------------------------
    R = [win_row;  bx_borrow];
    T = [ev_t;     bt_borrow];
    M = sortrows([R, T], [1 2]);
    R = M(:,1);  T = M(:,2);

    isStart = [true; R(2:end) ~= R(1:end-1)];
    gstart  = find(isStart);
    gend    = [gstart(2:end) - 1; numel(R)];
    n       = gend - gstart + 1;             % points used per pixel
    nIEI    = n - 1;                         % intervals per pixel

    % IEI mean via the telescoping identity mean(diff) = (t_end - t_first)/nIEI
    t_first = T(gstart);
    t_last  = T(gend);
    mu_d        = (t_last - t_first) ./ max(nIEI, 1);
    mu_d(n < 2) = 0;

    % IEI variance: sum of squared within-pixel diffs via accumarray,
    % excluding the diffs that straddle a pixel boundary.
    d      = diff(T);
    within = true(numel(T) - 1, 1);
    within(gend(1:end-1)) = false;
    glab   = cumsum(isStart);  glab = glab(1:end-1);   % group label per diff
    sum_d2 = accumarray(glab(within), d(within).^2, [A 1]);

    var_d        = (sum_d2 - nIEI .* mu_d.^2) ./ max(nIEI - 1, 1);
    var_d        = max(var_d, 0);            % clamp float noise
    var_d(n < 3) = 0;                        % need >= 2 intervals for a std

    valid = n >= K;                 % emit only once K events are assembled
    mu_d(~valid)  = NaN;            % NaN, not 0 — a 0 mean detonates CV downstream
    var_d(~valid) = NaN;
    
    mean_map(uid)  = mu_d;
    var_map(uid)   = var_d;
    valid_map(uid) = valid;

    % ----------------------------------------------------------------------
    % 4. Refresh the ring buffer: newest K real timestamps per active pixel
    %    (old retained timestamps + this window's events; older ones drop).
    % ----------------------------------------------------------------------
    Rb = [br;       win_row];
    Tb = [old_ts;   ev_t];
    Mb = sortrows([Rb, Tb], [1 2]);
    Rb = Mb(:,1);  Tb = Mb(:,2);

    sB    = [true; Rb(2:end) ~= Rb(1:end-1)];
    gsB   = find(sB);
    grpB  = cumsum(sB);
    nB    = accumarray(grpB, 1, [A 1]);
    posB  = (1:numel(Rb))' - gsB(grpB);      % 0-based position within pixel
    rankB = nB(grpB) - posB;                 % 1 = newest

    keep   = rankB <= K;
    col    = K - rankB(keep) + 1;            % newest -> column K
    newbuf = nan(A, K);
    newbuf(sub2ind([A K], grpB(keep), col)) = Tb(keep);
    iei.buf(uid, :) = newbuf;
end