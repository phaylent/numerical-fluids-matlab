%% RUN_ADVECTION  Time and space schemes for linear advection
%
%   Solves u_t + a*u_x = 0 on a periodic domain with every pairing of two time
%   schemes (Trapezoidal, RK4) and two spatial schemes (CD2, CMPT4), each on
%   four combinations of time step and grid size. Reports the error against the
%   exact solution and saves one four-panel figure per scheme pairing.
%
%   Run from this folder. Output goes to ./results. Requires MATLAB R2020a or newer.

clear; clc; close all;

%% Parameters
a    = 1;
tEnd = 20;
xLim = [-10 10];            % periodic domain; the pulse makes exactly one lap

dtList = [0.20, 0.10, 0.20, 0.10];
NList  = [  40,   40,   80,   80];

timeSchemes  = {'Trapezoidal', 'RK4'};
spaceSchemes = {'CD2', 'CMPT4'};

outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Run all combinations
nRuns = numel(timeSchemes) * numel(spaceSchemes) * numel(dtList);
TimeScheme  = strings(nRuns, 1);
SpaceScheme = strings(nRuns, 1);
[dtCol, NCol, CFLCol, ErrL2, ErrMax] = deal(zeros(nRuns, 1));

xFine = linspace(xLim(1), xLim(2), 1000).';
uFine = exact_advection(xFine, tEnd, a, xLim);

row = 0;
for it = 1:numel(timeSchemes)
    for is = 1:numel(spaceSchemes)

        fig = figure('Color', 'w', 'Position', [100 100 900 750]);

        for p = 1:numel(dtList)
            dt = dtList(p);
            N  = NList(p);

            [x, u, uExact, CFL] = solve_advection(timeSchemes{it}, spaceSchemes{is}, ...
                                                  N, dt, tEnd, a, xLim);

            dx  = (xLim(2) - xLim(1)) / N;
            err = u - uExact;

            row = row + 1;
            TimeScheme(row)  = timeSchemes{it};
            SpaceScheme(row) = spaceSchemes{is};
            dtCol(row)  = dt;
            NCol(row)   = N;
            CFLCol(row) = CFL;
            ErrL2(row)  = sqrt(dx * sum(err.^2));     % discrete L2 norm
            ErrMax(row) = max(abs(err));

            subplot(2, 2, p);
            plot(xFine, uFine, '-k'); hold on;
            plot(x, u, 'o'); hold off;
            grid on; ylim([-0.4 1.1]);
            xlabel('x'); ylabel('u');
            title(sprintf('\\Delta t = %.2f,  N = %d,  CFL = %.2f', dt, N, CFL));
            legend('Exact', 'Numerical', 'Location', 'northwest');
        end

        sgtitle(sprintf('%s + %s', timeSchemes{it}, spaceSchemes{is}));
        exportgraphics(fig, fullfile(outDir, sprintf('advection_%s_%s.png', ...
            lower(timeSchemes{it}), lower(spaceSchemes{is}))), 'Resolution', 200);
    end
end

%% Error summary
Results = table(TimeScheme, SpaceScheme, dtCol, NCol, CFLCol, ErrL2, ErrMax, ...
    'VariableNames', {'TimeScheme','SpaceScheme','dt','N','CFL','ErrL2','ErrMax'});
disp(Results);
writetable(Results, fullfile(outDir, 'advection_errors.csv'));
