%% PLOT_SCHEME_SPECTRA  Dispersion and dissipation of each scheme pairing
%
%   Each Fourier mode exp(i*k*x) is advected at exactly speed a. A numerical
%   scheme changes the speed (dispersion) and the amplitude (dissipation) of
%   each mode. For a mode with phi = k*dx, the spatial operator has the symbol
%       CD2    :  i * sin(phi)
%       CMPT4  :  i * 1.5*sin(phi) / (1 + 0.5*cos(phi))
%   and each time scheme turns this into a per-step amplification factor G.
%   The quantities plotted are
%       phase speed ratio  =  -angle(G) / (CFL*phi)       (1 = exact)
%       amplitude left     =  |G|^nSteps                  (1 = no damping)
%
%   Settings match the first panel of run_advection (N = 40, dt = 0.2).
%   Run from this folder. Output goes to ./results. Requires MATLAB R2020a or newer.

clear; clc; close all;

a = 1;  tEnd = 20;  xLim = [-10 10];
N0 = 40;  dt0 = 0.2;

dx     = (xLim(2) - xLim(1)) / N0;
CFL    = a * dt0 / dx;
nSteps = round(tEnd / dt0);

phi = linspace(0.01, pi, 400);      % phi = k*dx; phi = pi is the Nyquist wavenumber

symbol.CD2   = sin(phi);
symbol.CMPT4 = 1.5 * sin(phi) ./ (1 + 0.5*cos(phi));

timeSchemes  = {'Trapezoidal', 'RK4'};
spaceSchemes = {'CD2', 'CMPT4'};

colors = lines(2);
styles = {'-', '--'};

fig = figure('Color', 'w', 'Position', [100 100 950 400]);
for it = 1:2
    for is = 1:2
        z = -1i * CFL * symbol.(spaceSchemes{is});

        switch timeSchemes{it}
            case 'Trapezoidal'
                G = (1 + z/2) ./ (1 - z/2);
            case 'RK4'
                G = 1 + z + z.^2/2 + z.^3/6 + z.^4/24;
        end

        label = sprintf('%s + %s', timeSchemes{it}, spaceSchemes{is});

        subplot(1, 2, 1); hold on;
        plot(phi, -angle(G) ./ (CFL*phi), styles{it}, 'Color', colors(is,:), ...
             'LineWidth', 1.4, 'DisplayName', label);

        subplot(1, 2, 2); hold on;
        plot(phi, abs(G).^nSteps, styles{it}, 'Color', colors(is,:), ...
             'LineWidth', 1.4, 'DisplayName', label);
    end
end

subplot(1, 2, 1);
yline(1, 'k:', 'HandleVisibility', 'off');
grid on; xlabel('\phi = k\Delta x'); ylabel('Numerical / exact phase speed');
title(sprintf('Dispersion (N = %d, CFL = %.2f)', N0, CFL));
legend('Location', 'southwest'); ylim([0 1.1]);

subplot(1, 2, 2);
grid on; xlabel('\phi = k\Delta x'); ylabel(sprintf('Amplitude remaining after %d steps', nSteps));
title('Dissipation'); legend('Location', 'southwest'); ylim([0 1.05]);

outDir = 'results';
if ~exist(outDir, 'dir'), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, 'scheme_spectra.png'), 'Resolution', 200);