%% Part 2 — While loop to a target precision + plots (DROP-IN)
% Reaches a user-chosen number of significant figures without using true pi.
% Plots:
% 1) Scatter (geometry)       2) Convergence of \hat{\pi}
% 3) CI half-width vs N       4) Required N (and time) vs s.f.

clear; clc; close all;
rng('shuffle');

s_list = [2 3];          % precision targets to achieve (add 4 if you want)
batch  = 1000;           % smaller batch => more steps -> nicer curves
z      = 1.96;           % ~95% CI
MaxN   = 5e7;            % safety cap to avoid runaway loops

% For report summary (per s)
result  = struct('s',[],'pi_hat',[],'pi_sf',[],'N',[],'time',[],'h',[],'relH',[]);
Results = [];            % array of structs

% Record detailed history ONLY for this target (so plots 2–3 have a curve)
hist_s     = 3;          % <-- make history for the harder target
hist_N     = [];         % running N
hist_pi    = [];         % running \hat{\pi}
hist_relH  = [];         % running relative half-width h/|\hat{\pi}|

for s = s_list
    rtol = 0.5 * 10^(1 - s);    % relative CI half-width threshold for s s.f.
    N = 0; H = 0; relHalfWidth = inf; h = NaN; pi_hat = NaN;
    tStart = tic;

    % --- Plot 1: geometry scatter (only for first item in s_list) ---
    doPlotScatter = (s == s_list(1));
    if doPlotScatter
        figure(1); clf; hold on; axis equal; box on
        xlim([0 1]); ylim([0 1]); xlabel('x'); ylabel('y');
        title(sprintf('Monte Carlo points (target: %d s.f., 95%% CI rule)', s));
        th = linspace(0, pi/2, 200); plot(cos(th), sin(th), 'k--','LineWidth',1);
    end

    while relHalfWidth > rtol && N < MaxN
        xb = rand(batch,1);
        yb = rand(batch,1);
        in = xb.^2 + yb.^2 <= 1;

        H = H + sum(in);
        N = N + batch;

        if doPlotScatter
            scatter(xb(in),  yb(in),  6, 'filled');   % inside
            scatter(xb(~in), yb(~in), 6);             % outside
            drawnow limitrate
        end

        p_hat  = H/N;
        pi_hat = 4*p_hat;
        se_pi  = 4*sqrt(p_hat*(1 - p_hat)/N);
        h      = z * se_pi;
        relHalfWidth = h / abs(pi_hat);

        % Record history for plots 2–3
        if s == hist_s
            hist_N(end+1)    = N;            %#ok<AGROW>
            hist_pi(end+1)   = pi_hat;       %#ok<AGROW>
            hist_relH(end+1) = relHalfWidth; %#ok<AGROW>
        end
    end

    elapsed = toc(tStart);

    % Round to s significant figures (robust to older MATLAB versions)
    try
        pi_sf = round(pi_hat, s, 'significant');
    catch
        pow = floor(log10(abs(pi_hat)));
        pi_sf = round(pi_hat / 10^(pow - (s-1))) * 10^(pow - (s-1));
    end

    if doPlotScatter
        text(0.05,0.95, sprintf('\\pi \\approx %.*g (N=%d)', s, pi_sf, N), ...
             'Units','normalized','VerticalAlignment','top', ...
             'FontWeight','bold','FontSize',12);
    end

    fprintf('s=%d: pi ≈ %.*g, N=%d, time=%.3fs, 95%% CI half-width=%.3g, relHalfWidth=%.3g\n', ...
            s, s, pi_sf, N, elapsed, h, relHalfWidth);

    % Store summary row
    r = result; r.s = s; r.pi_hat = pi_hat; r.pi_sf = pi_sf;
    r.N = N; r.time = elapsed; r.h = h; r.relH = relHalfWidth;
    Results = [Results; r]; %#ok<AGROW>
end

% ---- Plot 2: Convergence of \hat{\pi} (for hist_s only) ----
if ~isempty(hist_N)
    figure(2); clf
    semilogx(hist_N, hist_pi, 'o-'); hold on
    yline(pi,'--'); grid on
    xlabel('N (samples)','Interpreter','latex');
    ylabel('$\hat{\pi}$','Interpreter','latex');
    title(sprintf('Convergence of $\\hat{\\pi}$ (target %d s.f.)', hist_s), 'Interpreter','latex');
    legend({'$\hat{\pi}$','true $\pi$'},'Interpreter','latex','Location','best');
end

% ---- Plot 3: CI half-width (relative) vs N + stop threshold ----
if ~isempty(hist_N)
    figure(3); clf
    loglog(hist_N, hist_relH, 'o-'); hold on; grid on
    rtol_line = 0.5 * 10^(1 - hist_s);
    yline(rtol_line,'--','Interpreter','latex');
    xlabel('N (samples)','Interpreter','latex');
    ylabel('Relative CI half-width $h/|\hat{\pi}|$','Interpreter','latex');
    title('CI shrinkage vs N (95\% interval)','Interpreter','latex');
    legend({'history','$\text{stop threshold}$'},'Interpreter','latex','Location','southwest');
end

% ---- Plot 4: Required N (and time) vs precision level ----
S  = [Results.s]';      % s.f. levels
NN = [Results.N]';      % samples used
TT = [Results.time]';   % seconds

figure(4); clf
subplot(1,2,1)
bar(S, NN)
xlabel('Significant figures (s)'); ylabel('Samples required (N)');
title('Samples needed to reach target precision'); grid on

subplot(1,2,2)
bar(S, TT)
xlabel('Significant figures (s)'); ylabel('Runtime (s)');
title('Runtime to reach target precision'); grid on
