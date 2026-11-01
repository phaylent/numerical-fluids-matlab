%% STUDY_GRID_CONVERGENCE  Convergence of the duct-flow solver against Fanno flow
%
%   Solves the forward problem on successively halved grids and measures the
%   largest error against the analytic Fanno solution for each variable.
%       observed order = log2( error(h) / error(h/2) )
%
%   Run from this folder. Output goes to ./results. Requires MATLAB R2020a or newer.

clear; clc; close all;

par.k       = 1.4;
par.LoverD  = 100;
par.f       = 0.01;
par.tol     = 1e-10;
par.maxIter = 50;
M1 = 0.5;

cells = [25 50 100 200 400];       % number of intervals; N = cells + 1 grid points
nC    = numel(cells);
err   = zeros(nC, 4);              % columns: M, P, rho, T

for n = 1:nC
    par.N = cells(n) + 1;
    sol = solve_duct(M1, par);
    [Me, Pe, rhoe, Te] = fanno_exact(sol.x, M1, par.k, par.f);
    err(n, :) = [max(abs(sol.M - Me)),     max(abs(sol.P - Pe)), ...
                 max(abs(sol.rho - rhoe)), max(abs(sol.T - Te))];
    fprintf('Finished %d intervals\n', cells(n));
end

order = log2(err(1:end-1, :) ./ err(2:end, :));

fprintf('\nMax error against Fanno solution (M1 = %.2f)\n', M1);
fprintf('%10s%14s%14s%14s%14s\n', 'intervals', 'M', 'P', 'rho', 'T');
for n = 1:nC
    fprintf('%10d%14.3e%14.3e%14.3e%14.3e\n', cells(n), err(n, :));
end

fprintf('\nObserved order\n');
fprintf('%10s%14s%14s%14s%14s\n', 'pair', 'M', 'P', 'rho', 'T');
for n = 1:nC-1
    fprintf('%10s%14.2f%14.2f%14.2f%14.2f\n', ...
            sprintf('%d->%d', cells(n), cells(n+1)), order(n, :));
end

fig = figure('Color', 'w', 'Position', [100 100 650 480]);
loglog(cells, err, 'o-', 'LineWidth', 1.2); hold on;
loglog(cells, err(1, 1) * (cells(1) ./ cells).^2, 'k--');
grid on; xticks(cells);
xlabel('Number of intervals'); ylabel('Max error against Fanno solution');
title('Grid convergence of the duct-flow solver');
legend('M', 'P', '\rho', 'T', '2nd-order slope', 'Location', 'southwest');

outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, 'duct_grid_convergence.png'), 'Resolution', 200);
writematrix([cells(:), err], fullfile(outDir, 'duct_grid_convergence.csv'));