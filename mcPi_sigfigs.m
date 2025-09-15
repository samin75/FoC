function [pi_rounded, N, ci] = mcPi_sigfigs(s, varargin)
% MCPI_SIGFIGS  Monte Carlo estimate of π to s significant figures.
% Uses a 95% CI stop rule that does NOT use the true π.
%
% Usage:
%   [pi_rounded, N, ci] = mcPi_sigfigs(s, 'Batch', 2000, 'MaxN', 5e7, 'Plot', true)
%
% Inputs:
%   s       : integer >= 1, number of significant figures to guarantee
% Options (Name-Value pairs):
%   'Batch' : number of random points per update (default 5000)
%   'MaxN'  : maximum total samples (default 1e8)
%   'Plot'  : true/false, whether to plot the random points (default true)
%   'Figure': figure number to draw into (default new window)
%
% Outputs:
%   pi_rounded : π rounded to s significant figures
%   N          : total number of points used
%   ci         : 95% confidence interval around unrounded π estimate

    % --- Interactive prompt if s not provided ---
    if nargin < 1
        s = input('Enter desired significant figures (e.g., 2, 3, 4): ');
        if ~isscalar(s) || ~isfinite(s) || s ~= fix(s) || s < 1
            error('Please enter a positive integer (e.g., 2, 3, 4).');
        end
    end

    % --- Parse inputs ---
    p = inputParser;
    addRequired(p, 's', @(x) isscalar(x) && x==fix(x) && x>=1);
    addParameter(p, 'Batch', 5000, @(x) isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'MaxN', 1e8,   @(x) isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'Figure', [],  @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Plot', true,  @(x) islogical(x) || isnumeric(x));
    parse(p, s, varargin{:});

    B      = p.Results.Batch;
    MaxN   = p.Results.MaxN;
    figArg = p.Results.Figure;
    doPlot = logical(p.Results.Plot);

    % --- Setup ---
    rng('shuffle');
    z    = 1.96;                  % ~95% z-score
    rtol = 0.5 * 10^(1 - s);      % relative CI half-width threshold

    N = 0; H = 0;
    relHalfWidth = inf; h = NaN; pi_hat = NaN;

    % --- Plot scaffolding with throttling/buffering ---
    if doPlot
        if isempty(figArg)
            fig = figure('Name', sprintf('MC \\pi — %d s.f. (95%% CI)', s), 'NumberTitle','off');
        else
            fig = figure(figArg);
            set(fig, 'Name', sprintf('MC \\pi — %d s.f. (95%% CI)', s), 'NumberTitle','off');
        end
        hold on; axis equal; box on
        xlim([0 1]); ylim([0 1]);
        xlabel('x'); ylabel('y');
        title(sprintf('Monte Carlo \\pi (target: %d s.f., 95%% CI rule)', s));
        th = linspace(0, pi/2, 200);
        plot(cos(th), sin(th), 'k--', 'LineWidth', 1);  % quarter circle

        % Create ONE scatter for inside & ONE for outside; fixed colors
        hIn  = scatter(nan, nan, 6, 'g', 'filled');
        hOut = scatter(nan, nan, 6, 'r');

        % Buffers (keep only recent points on screen)
        maxPlotPts = 5e4;   % show at most ~50k points per class
        xIn = []; yIn = []; xOut = []; yOut = [];

        % Redraw throttling
        drawEvery = 20;     % update figure once every 20 batches
        iter = 0;
    end

    % --- While loop: sample until CI is tight enough ---
    tStart = tic;
    while relHalfWidth > rtol && N < MaxN
        xb = rand(B,1); yb = rand(B,1);
        in = xb.^2 + yb.^2 <= 1;

        H = H + sum(in);
        N = N + B;

        if doPlot
            % Append to buffers
            xIn  = [xIn;  xb(in)];  yIn  = [yIn;  yb(in)];
            xOut = [xOut; xb(~in)]; yOut = [yOut; yb(~in)];

            % Trim to last maxPlotPts
            if numel(xIn) > maxPlotPts
                xIn  = xIn(end-maxPlotPts+1:end);
                yIn  = yIn(end-maxPlotPts+1:end);
            end
            if numel(xOut) > maxPlotPts
                xOut = xOut(end-maxPlotPts+1:end);
                yOut = yOut(end-maxPlotPts+1:end);
            end

            % Throttled redraw
            iter = iter + 1;
            if mod(iter, drawEvery) == 0
                set(hIn,  'XData', xIn,  'YData', yIn);
                set(hOut, 'XData', xOut, 'YData', yOut);
                drawnow limitrate
            end
        end

        % Update stats
        p_hat  = H/N;
        pi_hat = 4*p_hat;
        se_pi  = 4*sqrt(p_hat*(1 - p_hat)/N);
        h      = z * se_pi;
        relHalfWidth = h / abs(pi_hat);
    end
    elapsed = toc(tStart);

    % Final refresh so last batch is visible
    if doPlot
        set(hIn,  'XData', xIn,  'YData', yIn);
        set(hOut, 'XData', xOut, 'YData', yOut);
        drawnow
    end

    % --- Final rounding to s significant figures ---
    try
        pi_rounded = round(pi_hat, s, 'significant');
    catch
        pow = floor(log10(abs(pi_hat)));
        pi_rounded = round(pi_hat / 10^(pow - (s-1))) * 10^(pow - (s-1));
    end

    % Guard in case loop never iterated (shouldn't happen with sane MaxN)
    if isnan(h)
        p_hat  = H/max(N,1);
        pi_hat = 4*p_hat;
        se_pi  = 4*sqrt(p_hat*(1 - p_hat)/max(N,1));
        h      = z * se_pi;
    end

    % 95% CI around unrounded estimate
    ci = [pi_hat - h, pi_hat + h];

    % --- Print + annotate ---
    fprintf('pi ≈ %.*g (s=%d s.f.), N=%d, time=%.3fs, 95%% CI: [%.6f, %.6f]\n', ...
            s, pi_rounded, s, N, elapsed, ci(1), ci(2));

    if doPlot
        text(0.05, 0.95, sprintf('\\pi \\approx %.*g (N=%d)', s, pi_rounded, N), ...
            'Units','normalized','VerticalAlignment','top', ...
            'FontWeight','bold','FontSize',12);
    end
end
