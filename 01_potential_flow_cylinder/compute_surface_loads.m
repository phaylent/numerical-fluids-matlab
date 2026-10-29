function loads = compute_surface_loads(CpSurface, theta, Rmin, Uinf)
%COMPUTE_SURFACE_LOADS  Pressure forces on the cylinder (per unit span, rho = 1).
%
%   The pressure force on the body is F = -integral( p * n dA ) with outward
%   normal n = (cos(theta), sin(theta)). With the free-stream pressure removed,
%   p - p_inf = 0.5*rho*Uinf^2*Cp.
%
%   Returned fields:
%       Fy_top    y-force from the upper half (0 <= theta <= pi)
%       Fx_left   x-force from the left half  (pi/2 <= theta <= 3*pi/2)
%       Fx_total  x-force over the full cylinder
%       Fy_total  y-force over the full cylinder
%
%   For steady inviscid flow the full-cylinder forces vanish (d'Alembert's
%   paradox), so Fx_total and Fy_total are a check on the numerical solution.
%   Ntheta must be a multiple of 4 so the quarter-turn angles are grid nodes.

    theta  = theta(:);
    pExcess = 0.5 * Uinf^2 * CpSurface(:);
    dtheta  = 2*pi / numel(theta);

    integrandX = -pExcess .* cos(theta) * Rmin;
    integrandY = -pExcess .* sin(theta) * Rmin;

    tol = 1e-9;
    topMask  = theta >= -tol       & theta <= pi + tol;
    leftMask = theta >= pi/2 - tol & theta <= 3*pi/2 + tol;

    loads.Fy_top   = trapz(theta(topMask),  integrandY(topMask));
    loads.Fx_left  = trapz(theta(leftMask), integrandX(leftMask));
    loads.Fx_total = dtheta * sum(integrandX);    % periodic trapezoid rule
    loads.Fy_total = dtheta * sum(integrandY);
end