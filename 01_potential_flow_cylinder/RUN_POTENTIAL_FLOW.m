%% RUN_POTENTIAL_FLOW  Potential flow over a circular cylinder (finite differences)
%
%   Solves the Laplace equation for the stream function of uniform flow past a
%   cylinder, recovers the velocity and pressure fields, and measures how the
%   solution depends on (a) mesh resolution and (b) outer-domain size.
%
%   Pressure coefficient (Bernoulli):  Cp = 1 - (Vr^2 + Vtheta^2) / Uinf^2
%   Exact surface result, infinite domain:  Cp = 1 - 4*sin(theta)^2
%
%   Run from this folder. Figures and a CSV summary are written to ./results.
%   Requires MATLAB R2020a or newer (exportgraphics).

clear; clc; close all;

%% Parameters
prm.Rmin = 1;        % cylinder radius (length scale)
prm.Uinf = 1;        % free-stream speed (velocity scale)
prm.C0   = 0;        % psi on the cylinder surface
prm.Cinf = 0;        % constant offset of psi at the outer boundary
alpha    = 0.5;      % radial stretching parameter

assert(prm.C0 == 0 && prm.Cinf == 0, ...
    'Exact-solution comparisons assume C0 = Cinf = 0 (no circulation).');

% Cases: [Nr, Ntheta, Rmax]
%   Rows 1-4: mesh refinement at fixed domain (Rmax = 10)
%   Rows 5-7: domain growth at fixed mesh (160 x 80)
caseTable = [ ...
     20,  20, 10;
     40,  40, 10;
     80,  80, 10;
    160, 160, 10;
    160,  80, 10;
    160,  80, 20;
    160,  80, 40];

plotCases = [1 2 6];      % cases that get the four-panel figure

outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Solve all cases
nCases = size(caseTable, 1);
[PsiErr, CpMin, CpErrDisc, CpErrTotal, FyTop, FxLeft, FxTotal, FyTotal] = ...
    deal(zeros(nCases, 1));

for k = 1:nCases
    Nr = caseTable(k,1);  Nth = caseTable(k,2);  Rmax = caseTable(k,3);
    prm.Rmax = Rmax;
    fprintf('Case %d/%d: Nr = %d, Ntheta = %d, Rmax = %g\n', k, nCases, Nr, Nth, Rmax);

    % Stream function and velocity
    [psi, r, theta] = solve_laplace_cylinder(Nr, Nth, prm.Rmin, Rmax, ...
                                             prm.Uinf, prm.C0, prm.Cinf, alpha);
    [Vr, Vtheta, Vx, Vy] = compute_velocity_field(psi, theta, prm.Rmin, Rmax, alpha);

    % Pressure coefficient from Bernoulli's equation
    Cp     = 1 - (Vr.^2 + Vtheta.^2) / prm.Uinf^2;
    CpSurf = Cp(:, 1);

    % Reference solutions
    [psiExact, A] = exact_stream_function(r, theta, prm.Rmin, Rmax, prm.Uinf);
    CpExactInf = 1 - 4*sin(theta).^2;                         % Rmax -> infinity
    CpExactFin = 1 - 4*(A/prm.Uinf)^2 * sin(theta).^2;        % this finite domain

    % Errors and forces
    PsiErr(k)      = max(abs(psi(:) - psiExact(:)));
    CpMin(k)       = min(CpSurf);
    CpErrDisc(k)   = max(abs(CpSurf - CpExactFin));   % discretization error only
    CpErrTotal(k)  = max(abs(CpSurf - CpExactInf));   % discretization + domain truncation
    loads          = compute_surface_loads(CpSurf, theta, prm.Rmin, prm.Uinf);
    FyTop(k)       = loads.Fy_top;
    FxLeft(k)      = loads.Fx_left;
    FxTotal(k)     = loads.Fx_total;
    FyTotal(k)     = loads.Fy_total;

    % Figures
    if ismember(k, plotCases)
        sol = struct('r', r, 'theta', theta, 'psi', psi, 'Vx', Vx, 'Vy', Vy, 'Cp', Cp);
        caseTitle = sprintf('N_r = %d, N_\\theta = %d, R_\\infty = %g', Nr, Nth, Rmax);
        fig = plot_case_summary(sol, prm, caseTitle);
        exportgraphics(fig, fullfile(outDir, ...
            sprintf('case_Nr%d_Nth%d_Rmax%d.png', Nr, Nth, Rmax)), 'Resolution', 200);
    end
end

%% Summary table
Results = table(caseTable(:,1), caseTable(:,2), caseTable(:,3), PsiErr, CpMin, ...
                CpErrDisc, CpErrTotal, FyTop, FxLeft, FxTotal, FyTotal, ...
                'VariableNames', {'Nr','Ntheta','Rmax','PsiErr','CpMin', ...
                'CpErrDisc','CpErrTotal','FyTop','FxLeft','FxTotal','FyTotal'});
disp(Results);
writetable(Results, fullfile(outDir, 'results_summary.csv'));

fprintf('\nReferences (infinite domain, rho = 1):\n');
fprintf('  Cp,min = -3\n');
fprintf('  Fy (upper half) = %.4f,  Fx (left half) = %.4f\n', ...
        (5/3)*prm.Uinf^2*prm.Rmin, -(1/3)*prm.Uinf^2*prm.Rmin);
fprintf('  Full-cylinder Fx, Fy = 0 (d''Alembert''s paradox)\n');

%% Convergence figure
idxA = caseTable(:,3) == 10 & caseTable(:,1) == caseTable(:,2);   % mesh refinement
idxB = caseTable(:,1) == 160 & caseTable(:,2) == 80;              % domain growth

fig = figure('Color', 'w', 'Position', [100 100 950 400]);

subplot(1,2,1);
NrA = caseTable(idxA, 1);  errA = CpErrDisc(idxA);
loglog(NrA, errA, 'o-', 'LineWidth', 1.2); hold on;
loglog(NrA, errA(1) * (NrA(1) ./ NrA).^2, 'k--');
grid on; xlabel('N_r = N_\theta');
ylabel('Max surface C_p error vs. finite-domain exact');
title('Mesh refinement (R_\infty = 10)');
legend('Numerical', 'Second-order slope', 'Location', 'southwest');

subplot(1,2,2);
RmaxB    = caseTable(idxB, 3);
CpMinFin = 1 - 4*(RmaxB.^2 ./ (RmaxB.^2 - prm.Rmin^2)).^2;
plot(RmaxB, CpMin(idxB), 'o-', 'LineWidth', 1.2); hold on;
plot(RmaxB, CpMinFin, 's--', 'LineWidth', 1.2);
yline(-3, 'k:');
grid on; xlabel('R_\infty'); ylabel('Minimum surface C_p');
title('Domain size (160 x 80 mesh)');
legend('Numerical', 'Exact, finite domain', 'Exact, infinite domain', 'Location', 'east');

exportgraphics(fig, fullfile(outDir, 'convergence_studies.png'), 'Resolution', 200);