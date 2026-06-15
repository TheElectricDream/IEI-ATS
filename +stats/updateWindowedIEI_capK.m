function [mean_map, var_map, valid_map, iei] = updateWindowedIEI_capK(ev_x, ev_y, ev_t, imgSz, K, iei)
% Constant-sample-size variant: every pixel's IEI statistics are computed
% from EXACTLY the K most-recent events at that pixel (current window first,
% backfilled from the ring buffer only as needed). Dense pixels no longer get
% a larger sample than sparse ones, so mean/var/CV are directly comparable.
% A pixel is only marked valid once K events can actually be assembled.

    nPix = prod(imgSz);
    if ~isfield(iei,'buf') || isempty(iei.buf) || size(iei.buf,2) ~= K
        iei.buf = nan(nPix, K);
    end

    mean_map  = nan(imgSz);
    var_map   = nan(imgSz);
    valid_map = false(imgSz);
    if isempty(ev_t), return; end

    ev_x = ev_x(:);  ev_y = ev_y(:);  ev_t = ev_t(:);

    % ----------------------------------------------------------------------
    % 1. Index the window's events by active pixel.
    % ----------------------------------------------------------------------
    lin = sub2ind(imgSz, ev_x, ev_y);
    uid = unique(lin, 'stable');
    A   = numel(uid);
    row_of_pix      = zeros(nPix, 1);
    row_of_pix(uid) = (1:A)';
    win_row = row_of_pix(lin);
    counts  = accumarray(win_row, 1, [A 1]);

    % ----------------------------------------------------------------------
    % 2. Cap the window contribution to the newest K events per pixel.
    %    (Only matters for counts > K; for counts <= K every event is kept.)
    %    NOTE: the FULL window events are still used for the buffer refresh
    %    in step 6 -- the cap only shapes the statistics sample.
    % ----------------------------------------------------------------------
    Mw   = sortrows([win_row, ev_t], [1 2]);   % ascending time within pixel
    wr   = Mw(:,1);  wt = Mw(:,2);
    sW   = [true; wr(2:end) ~= wr(1:end-1)];
    gsW  = find(sW);
    grpW = cumsum(sW);
    nW   = accumarray(grpW, 1, [A 1]);
    posW = (1:numel(wr))' - gsW(grpW);         % 0-based position (oldest=0)
    rankW = nW(grpW) - posW;                    % 1 = newest
    keepW = rankW <= K;
    win_row_s = wr(keepW);                      % capped window stream
    ev_t_s    = wt(keepW);

    % ----------------------------------------------------------------------
    % 3. Backfill from the ring buffer: newest (K - min(counts,K)) per pixel.
    %    need == 0 wherever counts >= K, so dense pixels borrow nothing.
    % ----------------------------------------------------------------------
    need = K - min(counts, K);                  % == max(0, K - counts)

    oldbuf  = iei.buf(uid, :);
    [br, bc] = find(~isnan(oldbuf));
    br = br(:);  bc = bc(:);
    old_ts  = oldbuf(sub2ind([A K], br, bc));
    old_ts  = old_ts(:);

    rank_from_new = K - bc + 1;                  % 1 = newest
    borrow = rank_from_new <= need(br);
    bx_borrow = br(borrow);
    bt_borrow = old_ts(borrow);

    % ----------------------------------------------------------------------
    % 4. Statistics stream = capped window events + borrowed, grouped per
    %    pixel and sorted in time (borrowed are older, so they fall first).
    %    Every supplied pixel now contributes exactly K points => K-1 IEIs.
    % ----------------------------------------------------------------------
    R = [win_row_s;  bx_borrow];
    T = [ev_t_s;     bt_borrow];
    M = sortrows([R, T], [1 2]);
    R = M(:,1);  T = M(:,2);

    isStart = [true; R(2:end) ~= R(1:end-1)];
    gstart  = find(isStart);
    gend    = [gstart(2:end) - 1; numel(R)];
    n       = gend - gstart + 1;                % points used per pixel (<= K)
    nIEI    = n - 1;

    % IEI mean via telescoping identity: mean(diff) = (t_end - t_first)/nIEI
    t_first = T(gstart);
    t_last  = T(gend);
    mu_d        = (t_last - t_first) ./ max(nIEI, 1);
    mu_d(n < 2) = 0;

    % IEI sample variance: (sum d^2 - nIEI*mu^2)/(nIEI-1), within-pixel only.
    d      = diff(T);
    within = true(numel(T) - 1, 1);
    within(gend(1:end-1)) = false;
    glab   = cumsum(isStart);  glab = glab(1:end-1);
    sum_d2 = accumarray(glab(within), d(within).^2, [A 1]);

    var_d  = (sum_d2 - nIEI .* mu_d.^2) ./ max(nIEI - 1, 1);
    var_d  = max(var_d, 0);

    % ----------------------------------------------------------------------
    % 5. Emit only where the full K-event sample was assembled.
    % ----------------------------------------------------------------------
    valid = n >= K;                             % constant sample size reached
    mu_d(~valid)  = NaN;                         % keep 0-mean out of CV = std/mean
    var_d(~valid) = NaN;

    mean_map(uid)  = mu_d;
    var_map(uid)   = var_d;
    valid_map(uid) = valid;

    % ----------------------------------------------------------------------
    % 6. Refresh the ring buffer with the newest K REAL timestamps per pixel
    %    (full window events + retained buffer; the stats cap does not apply).
    % ----------------------------------------------------------------------
    Rb = [br;       win_row];
    Tb = [old_ts;   ev_t];
    Mb = sortrows([Rb, Tb], [1 2]);
    Rb = Mb(:,1);  Tb = Mb(:,2);

    sB    = [true; Rb(2:end) ~= Rb(1:end-1)];
    gsB   = find(sB);
    grpB  = cumsum(sB);
    nB    = accumarray(grpB, 1, [A 1]);
    posB  = (1:numel(Rb))' - gsB(grpB);
    rankB = nB(grpB) - posB;

    keep   = rankB <= K;
    col    = K - rankB(keep) + 1;
    newbuf = nan(A, K);
    newbuf(sub2ind([A K], grpB(keep), col)) = Tb(keep);
    iei.buf(uid, :) = newbuf;
end