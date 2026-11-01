function F = fanno_parameter(M, k)
%FANNO_PARAMETER  Friction parameter f*L*/D for adiabatic flow with friction.
%
%   For a subsonic Mach number M, F is the dimensionless duct length (with the
%   Darcy friction factor f) over which friction would bring the flow to M = 1:
%       F = (1 - M^2)/(k*M^2) + (k+1)/(2k) * ln( (k+1)*M^2 / (2 + (k-1)*M^2) )
%   F decreases monotonically from infinity at M = 0 to 0 at M = 1.

    F = (1 - M.^2) ./ (k*M.^2) ...
        + (k+1)/(2*k) * log((k+1)*M.^2 ./ (2 + (k-1)*M.^2));
end