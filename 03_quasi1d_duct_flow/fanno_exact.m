function [M, P, rho, T] = fanno_exact(x, M1, k, f)
%FANNO_EXACT  Analytic subsonic Fanno-flow profiles along the duct.
%
%   Inputs:  x   positions in units of hydraulic diameter (column or row)
%            M1  inlet Mach number
%            k   specific heat ratio
%            f   Darcy friction factor
%
%   At each x the Mach number satisfies
%       fanno_parameter(M) = fanno_parameter(M1) - f*x
%   Pressure, density, and temperature then follow from the Fanno ratios
%   (relative to the sonic reference state), re-normalized by their inlet values.
%   Outputs are column vectors normalized by inlet values (M is the Mach number).

    x   = x(:);
    FL1 = fanno_parameter(M1, k);
    opts = optimset('TolX', 1e-14);

    M = zeros(size(x));
    for j = 1:numel(x)
        target = FL1 - f * x(j);
        if target <= 0
            error('fanno_exact:choked', ...
                'Duct is choked: f*x/D = %.4f exceeds the limit %.4f for M1 = %.4f.', ...
                f * x(j), FL1, M1);
        end
        if x(j) == 0
            M(j) = M1;
        else
            M(j) = fzero(@(m) fanno_parameter(m, k) - target, [M1, 1], opts);
        end
    end

    pStar = @(m) (1 ./ m) .* sqrt((k+1) ./ (2 + (k-1)*m.^2));     % P / P*
    tStar = @(m) (k+1) ./ (2 + (k-1)*m.^2);                       % T / T*
    rStar = @(m) (1 ./ m) .* sqrt((2 + (k-1)*m.^2) / (k+1));      % rho / rho*

    P   = pStar(M) / pStar(M1);
    T   = tStar(M) / tStar(M1);
    rho = rStar(M) / rStar(M1);
end