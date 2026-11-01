%% RUN_DUCT_FLOW  Compressible flow with friction in a constant-area duct
%
%   Solves the steady quasi-1D equations (continuity, momentum with wall
%   friction, adiabatic energy, ideal-gas law) on N grid points with Newton's
%   method, and checks the result against the analytic Fanno-flow solution.
%
%   Task 1 (forward):    prescribe the inlet Mach number M1 and compute the flow.
%   Task 2 (backward):   prescribe the exit pressure ratio P_N and find M1.
%   Task 3 (parametric): repeat Task 2 for P_N = 0.60 ... 0.90.
%
%   Run from this folder. Output goes to ./results. Requires MATLAB R2020a or newer.

clear; clc; close all;

%% Parameters
par.N       = 200;       % grid points along the duct
par.k       = 1.4;       % specific heat ratio
par.LoverD  = 100;       % duct length in hydraulic diameters (x is in units of D_h)
par.f       = 0.01;      % Darcy friction factor
par.tol     = 1e-10;     % Newton tolerance: largest change in any unknown
par.maxIter = 50;
par.tolM1   = 1e-9;      % tolerance for the inlet-Mach iteration

M1Forward = 0.5;
PNTarget  = 0.6;
PNList    = 0.6:0.01:0.9;

outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Choking limit (analytic)
fL      = par.f * par.LoverD;
M1max   = fanno_choke_mach(par.k, fL);
PNchoke = 1 / ((1/M1max) * sqrt((par.k+1) / (2 + (par.k-1)*M1max^2)));
fprintf('Friction parameter f*L/D = %.4f\n', fL);
fprintf('Largest inlet Mach number without choking: %.4f\n', M1max);
fprintf('Lowest reachable exit pressure ratio (choked exit): %.4f\n', PNchoke);

%% Task 1: forward problem
fprintf('\n--- Task 1: forward problem, M1 = %.4f ---\n', M1Forward);
sol1 = solve_duct(M1Forward, par);
[Me, Pe, rhoe, Te] = fanno_exact(sol1.x, M1Forward, par.k, par.f);
ref1 = struct('x', sol1.x, 'M', Me, 'P', Pe, 'rho', rhoe, 'T', Te);

fprintf('Newton iterations: %d (last step %.2e), final residual %.2e\n', ...
        sol1.info.iterations, sol1.info.stepNorm, sol1.info.resNorm);
fprintf('Exit values       numerical     Fanno exact\n');
fprintf('  Mach number     %.6f      %.6f\n', sol1.M(end),   Me(end));
fprintf('  P/P_inlet       %.6f      %.6f\n', sol1.P(end),   Pe(end));
fprintf('  rho/rho_inlet   %.6f      %.6f\n', sol1.rho(end), rhoe(end));
fprintf('  T/T_inlet       %.6f      %.6f\n', sol1.T(end),   Te(end));
errs = [max(abs(sol1.M - Me)), max(abs(sol1.P - Pe)), ...
        max(abs(sol1.rho - rhoe)), max(abs(sol1.T - Te))];
fprintf('Max error over the duct: M %.2e, P %.2e, rho %.2e, T %.2e\n', errs);

fig = plot_duct_profiles(sol1, ref1, sprintf('Forward problem: M_1 = %.2f', M1Forward));
exportgraphics(fig, fullfile(outDir, 'duct_forward_profiles.png'), 'Resolution', 200);

fig = figure('Color', 'w', 'Position', [100 100 650 450]);
ix = 2:par.N;
semilogy(sol1.x(ix), abs(sol1.M(ix)   - Me(ix)),   ...
         sol1.x(ix), abs(sol1.P(ix)   - Pe(ix)),   ...
         sol1.x(ix), abs(sol1.rho(ix) - rhoe(ix)), ...
         sol1.x(ix), abs(sol1.T(ix)   - Te(ix)), 'LineWidth', 1.3);
grid on; xlabel('x/D_h'); ylabel('|numerical - exact|');
title('Error against the Fanno solution (forward problem)');
legend('M', 'P', '\rho', 'T', 'Location', 'northwest');
exportgraphics(fig, fullfile(outDir, 'duct_forward_error.png'), 'Resolution', 200);

%% Task 2: backward problem
fprintf('\n--- Task 2: backward problem, P_N = %.2f ---\n', PNTarget);
[M1b, sol2, nIter] = solve_inlet_mach(PNTarget, par, M1Forward, sol1.q);
M1bExact = fzero(@(m) fanno_exit_pressure_ratio(m, par.k, fL) - PNTarget, ...
                 [0.05, 0.999*M1max], optimset('TolX', 1e-14));
fprintf('Converged in %d iterations\n', nIter);
fprintf('Inlet Mach number: numerical %.6f, Fanno exact %.6f\n', M1b, M1bExact);
fprintf('Exit pressure ratio from the converged solution: %.8f (target %.2f)\n', ...
        sol2.P(end), PNTarget);

%% Task 3: parametric study
fprintf('\n--- Task 3: parametric study ---\n');
nP      = numel(PNList);
M1Num   = zeros(nP, 1);
M1Exact = zeros(nP, 1);
M1Guess = M1Forward;
qGuess  = sol1.q;

for j = 1:nP
    [M1Num(j), solj] = solve_inlet_mach(PNList(j), par, M1Guess, qGuess);
    M1Exact(j) = fzero(@(m) fanno_exit_pressure_ratio(m, par.k, fL) - PNList(j), ...
                       [0.05, 0.999*M1max], optimset('TolX', 1e-14));
    M1Guess = M1Num(j);       % continuation: start the next case from this one
    qGuess  = solj.q;
end

Parametric = table(PNList(:), M1Num, M1Exact, abs(M1Num - M1Exact), ...
    'VariableNames', {'PN', 'M1_numerical', 'M1_exact', 'AbsError'});
disp(Parametric);
writetable(Parametric, fullfile(outDir, 'duct_parametric_study.csv'));

fig = figure('Color', 'w', 'Position', [100 100 650 480]);
plot(PNList, M1Exact, 'k-', 'LineWidth', 1.2); hold on;
plot(PNList, M1Num, 'o', 'LineWidth', 1.2);
grid on; xlabel('P_{exit}/P_{inlet}'); ylabel('Inlet Mach number M_1');
title('Inlet Mach number required for a given exit pressure');
legend('Fanno exact', 'Numerical', 'Location', 'northeast');
exportgraphics(fig, fullfile(outDir, 'duct_parametric_study.png'), 'Resolution', 200);