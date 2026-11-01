function PN = fanno_exit_pressure_ratio(M1, k, fL)
%FANNO_EXIT_PRESSURE_RATIO  Analytic exit/inlet pressure ratio of a Fanno duct.
%
%   fL is the total friction parameter f*L/D (Darcy f). Returns NaN if the
%   flow would choke before the exit.

    target = fanno_parameter(M1, k) - fL;
    if target <= 0
        PN = NaN;
        return
    end

    M2    = fzero(@(m) fanno_parameter(m, k) - target, [M1, 1], optimset('TolX', 1e-14));
    pStar = @(m) (1 ./ m) .* sqrt((k+1) ./ (2 + (k-1)*m.^2));
    PN    = pStar(M2) / pStar(M1);
end