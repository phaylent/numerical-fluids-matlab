function fig = plot_duct_profiles(sol, ref, titleStr)
%PLOT_DUCT_PROFILES  Four-panel plot of M, P, rho, T along the duct.
%
%   sol: numerical solution (fields x, M, P, rho, T)
%   ref: analytic solution with the same fields, or [] to omit
%   P, rho, and T are normalized by their inlet values.

    fields  = {'M', 'P', 'rho', 'T'};
    names   = {'Mach number', 'Pressure', 'Density', 'Temperature'};
    ylabels = {'M', 'P/P_{inlet}', '\rho/\rho_{inlet}', 'T/T_{inlet}'};

    fig = figure('Color', 'w', 'Position', [100 100 900 700]);
    for n = 1:4
        subplot(2, 2, n);
        plot(sol.x, sol.(fields{n}), '-', 'LineWidth', 1.6); hold on;
        if ~isempty(ref)
            plot(ref.x, ref.(fields{n}), 'k--', 'LineWidth', 1.0);
        end
        grid on; xlabel('x/D_h'); ylabel(ylabels{n}); title(names{n});
        if n == 1 && ~isempty(ref)
            legend('Numerical', 'Fanno exact', 'Location', 'northwest');
        end
    end
    sgtitle(titleStr);
end