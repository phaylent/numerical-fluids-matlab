%% STUDY_CONVERGENCE  Observed order of accuracy of each scheme pairing
%
%   Refines the grid by factors of two while keeping CFL fixed, so dt shrinks
%   in proportion to dx. The error of the combined scheme then falls at the
%   rate of whichever of the spatial or temporal errors dominates:
%       observed order = log2( error(N) / error(2N) )
%
%   Run from this folder. Output goes to ./results. Requires MATLAB R2020a or newer.

clear; clc; close all;

a    = 1;
tEnd = 20;
xLim = [-10 10];
CFL  = 0.4;
NList = [40 80 160 320 640];

timeSchemes  = {'Trapezoidal', 'RK4'};
spaceSchemes = {'CD2', 'CMPT4'};
comboNames   = {'Trap+CD2', 'Trap+CMPT4', 'RK4+CD2', 'RK4+CMPT4'};

outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end

nN = numel(NList);
ErrL2 = zeros(nN, 4);

for n = 1:nN
    dx = (xLim(2) - xLim(1)) / NList(n);
    dt = CFL * dx / a;
    c = 0;
    for it = 1:2
        for is = 1:2
            c = c + 1;
            [~, u, uExact] = solve_advection(timeSchemes{it}, spaceSchemes{is}, ...
                                             NList(n), dt, tEnd, a, xLim);
            ErrL2(n, c) = sqrt(dx * sum((u - uExact).^2));
        end
    end
    fprintf('Finished N = %d\n', NList(n));
end

Order = log2(ErrL2(1:end-1, :) ./ ErrL2(2:end, :));

%% Print tables
fprintf('\nL2 error at t = %g, CFL = %.2f\n', tEnd, CFL);
fprintf('%6s', 'N'); fprintf('%14s', comboNames{:}); fprintf('\n');
for n = 1:nN
    fprintf('%6d', NList(n)); fprintf('%14.3e', ErrL2(n, :)); fprintf('\n');
end

fprintf('\nObserved order (between successive grids)\n');
fprintf('%10s', 'N pair'); fprintf('%14s', comboNames{:}); fprintf('\n');
for n = 1:nN-1
    fprintf('%10s', sprintf('%d->%d', NList(n), NList(n+1)));
    fprintf('%14.2f', Order(n, :)); fprintf('\n');
end

%% Plot
fig = figure('Color', 'w', 'Position', [100 100 650 500]);
loglog(NList, ErrL2, 'o-', 'LineWidth', 1.2); hold on;
loglog(NList, ErrL2(1,1) * (NList(1)./NList).^2, 'k--');
loglog(NList, ErrL2(1,4) * (NList(1)./NList).^4, 'k:');
grid on; xticks(NList);
xlabel('N (grid points)'); ylabel('L_2 error at t = 20');
title(sprintf('Convergence at fixed CFL = %.1f', CFL));
legend([comboNames, {'2nd-order slope', '4th-order slope'}], 'Location', 'southwest');

exportgraphics(fig, fullfile(outDir, 'advection_convergence.png'), 'Resolution', 200);
writematrix([NList(:), ErrL2], fullfile(outDir, 'advection_convergence.csv'));