function iei = updateWindowedIEI(iei, t_mean_diff, t_std_diff, ...
    t_min, t_2last, t_max, counts)
% UPDATEWINDOWEDIEI  Per-window IEI statistics with minimal backfill.
%
%   The .mean/.var maps are REBUILT FROM ZERO every call. A pixel gets
%   a nonzero statistic only if it formed one THIS window. Nothing is
%   held or accumulated across windows.
%
%     counts >= 3 : pure within-window statistics.
%     counts == 2 : not enough for a variance (one IEI). Backfill the
%                   single most recent retained event -> 3 points.
%     counts == 1 : backfill the two most recent retained events
%                   -> 3 points.
%     counts == 0, or insufficient history to reach 3 points:
%                   .mean/.var stay zero this window.
%
%   The only persistent state is the two-timestamp history (tp1, tp2),
%   which always holds the two most recent REAL events per pixel.
%   Backfilled events are used transiently to form this window's
%   statistic and are never re-inserted into any event stream.
%
%   State (all [imgSz]):
%     iei.mean / iei.var - THIS window's statistics (zero elsewhere)
%     iei.valid          - mask of pixels that formed a stat THIS window
%     iei.tp1 / iei.tp2  - newest / 2nd-newest real event timestamp
%                          (init NaN)
%
%   Inputs: outputs of stats.computeNeighborhoodStats plus t_2last
%   (2nd-latest timestamp this window, valid where counts >= 2) and
%   per-pixel counts.

new_mean = zeros(size(counts));      % fresh every window
new_var  = zeros(size(counts));
made     = false(size(counts));

% --- counts >= 3: pure within-window -------------------------------
dense = counts >= 3;
new_mean(dense) = t_mean_diff(dense);
new_var(dense)  = t_std_diff(dense).^2;
made(dense)     = true;

% --- counts == 2: backfill one retained event ----------------------
twoOK = (counts == 2) & ~isnan(iei.tp1);
d1 = t_min(twoOK) - iei.tp1(twoOK);
d2 = t_max(twoOK) - t_min(twoOK);
mB = (d1 + d2) ./ 2;
new_mean(twoOK) = mB;
new_var(twoOK)  = (d1 - mB).^2 + (d2 - mB).^2;   % n-1 = 1 (Bessel)
made(twoOK)     = true;

% --- counts == 1: backfill two retained events ---------------------
oneOK = (counts == 1) & ~isnan(iei.tp1) & ~isnan(iei.tp2);
e1 = iei.tp1(oneOK) - iei.tp2(oneOK);
e2 = t_max(oneOK)   - iei.tp1(oneOK);
mC = (e1 + e2) ./ 2;
new_mean(oneOK) = mC;
new_var(oneOK)  = (e1 - mC).^2 + (e2 - mC).^2;
made(oneOK)     = true;

iei.mean  = new_mean;
iei.var   = new_var;
iei.valid = made;                    % this window only, NOT sticky

% --- history update (AFTER use): two most recent REAL events -------
ge2 = counts >= 2;
iei.tp2(ge2) = t_2last(ge2);
iei.tp1(ge2) = t_max(ge2);

eq1 = counts == 1;
iei.tp2(eq1) = iei.tp1(eq1);         % shift old newest down
iei.tp1(eq1) = t_max(eq1);
end