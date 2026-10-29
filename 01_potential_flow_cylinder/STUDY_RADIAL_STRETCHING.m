%% STUDY_RADIAL_STRETCHING  Choosing the radial stretching parameter alpha
%
%   Stronger stretching reaches farther with fewer points, but makes the cell
%   size change faster from one cell to the next, which can degrade accuracy.
%   This script sweeps alpha on an 80 x 80 mesh (Rmax = 8), records the
%   largest cell-to-cell growth in the radial spacing and the maximum error of
%   psi against the exact solution, and picks the most accurate alpha that
%   keeps cell growth under a chosen limit.

clear; clc; close all;

Rmin = 1;  Rmax = 8;  Uinf = 1;  C0 = 0;  Cinf = 0;
Nr = 80;   Ntheta = 80;

alphaList    = 0.5:0.5:5.0;
stretchLimit = 0.08;       % design choice: max 8% change between adjacent cells

nA         = numel(alphaList);
maxStretch = zeros(nA, 1);
psiErr     = zeros(nA, 1);

for k = 1:nA
    [psi, r, theta] = solve_laplace_cylinder(Nr, Ntheta, Rmin, Rmax, Uinf, C0, Cinf, alphaList(k));

    dr = diff(r);
    maxStretch(k) = max(abs(diff(dr) ./ dr(1:end-1)));

    psiExact  = exact_stream_function(r, theta, Rmin, Rmax, Uinf);
    psiErr(k) = max(abs(psi(:) - psiExact(:)));

    fprintf('alpha = %.2f   max stretch = %.4f   max error = %.6e\n', ...
            alphaList(k), maxStretch(k), psiErr(k));
end

candidates = find(maxStretch <= stretchLimit);
if isempty(candidates)
    error('No alpha satisfies the stretch limit of %.2f.', stretchLimit);
end
[~, m] = min(psiErr(candidates));
bestAlpha = alphaList(candidates(m));
fprintf('\nBest alpha within the stretch limit = %.2f\n', bestAlpha);

fig = figure('Color', 'w', 'Position', [100 100 900 380]);
subplot(1,2,1);
plot(alphaList, maxStretch, 'o-'); hold on; yline(stretchLimit, 'k--');
grid on; xlabel('\alpha'); ylabel('Max cell-to-cell growth'); title('Mesh smoothness');
subplot(1,2,2);
semilogy(alphaList, psiErr, 'o-');
grid on; xlabel('\alpha'); ylabel('Max |\psi - \psi_{exact}|'); title('Accuracy');

if ~exist('results', 'dir'), mkdir('results'); end
exportgraphics(fig, fullfile('results', 'stretching_study.png'), 'Resolution', 200);