function fig = plot_case_summary(sol, prm, caseTitle)
%PLOT_CASE_SUMMARY  Four-panel summary of one solved case.
%
%   sol: struct with fields r, theta, psi, Vx, Vy, Cp  (see run_potential_flow)
%   prm: struct with fields Rmin, Rmax, Uinf
%
%   Panels: stream function contours, velocity vectors, Cp field, and the
%   surface Cp compared with both the infinite-domain and finite-domain
%   exact solutions.

    [Ntheta, Nr] = size(sol.psi);
    lim = [-3 3];

    Xg = cos(sol.theta) * sol.r.';
    Yg = sin(sol.theta) * sol.r.';

    % Repeat the first theta row to close the periodic seam in contour/pcolor
    Xc   = [Xg;       Xg(1,:)];
    Yc   = [Yg;       Yg(1,:)];
    psiC = [sol.psi;  sol.psi(1,:)];
    CpC  = [sol.Cp;   sol.Cp(1,:)];

    fig = figure('Color', 'w', 'Position', [100 100 900 800]);

    subplot(2,2,1);
    contour(Xc, Yc, psiC, linspace(-3, 3, 40));
    axis equal; grid on; xlim(lim); ylim(lim);
    xlabel('x'); ylabel('y'); title('Stream function \psi');

    subplot(2,2,2);
    iq = 1:max(1, round(Ntheta/40)):Ntheta;      % thin out vectors on fine meshes
    jq = 1:max(1, round(Nr/25)):Nr;
    quiver(Xg(iq,jq), Yg(iq,jq), sol.Vx(iq,jq), sol.Vy(iq,jq), 0.8);
    axis equal; grid on; xlim(lim); ylim(lim);
    xlabel('x'); ylabel('y'); title('Velocity field (V_x, V_y)');

    subplot(2,2,3);
    pcolor(Xc, Yc, CpC); shading interp; colorbar;
    caxis([-3 1]);
    axis equal; grid on; xlim(lim); ylim(lim);
    xlabel('x'); ylabel('y'); title('Pressure coefficient C_p');

    subplot(2,2,4);
    Anorm  = prm.Rmax^2 / (prm.Rmax^2 - prm.Rmin^2);   % A/Uinf for the finite domain
    thFine = linspace(0, 2*pi, 400);
    plot(sol.theta, sol.Cp(:,1), 'o'); hold on;
    plot(thFine, 1 - 4*sin(thFine).^2, '-');
    plot(thFine, 1 - 4*Anorm^2*sin(thFine).^2, '--');
    grid on; xlabel('\theta (rad)'); ylabel('C_p');
    title('C_p on cylinder surface (r = R_{min})');
    legend('Numerical', ...
           'Exact, infinite domain', ...
           sprintf('Exact, R_\\infty = %g', prm.Rmax), ...
           'Location', 'south');

    sgtitle(caseTitle);
end