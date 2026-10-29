function [Vr, Vtheta, Vx, Vy] = compute_velocity_field(psi, theta, Rmin, Rmax, alpha)
%COMPUTE_VELOCITY_FIELD  Velocity components from the stream function.
%
%   In polar coordinates:
%       Vr     =  (1/r) * dpsi/dtheta
%       Vtheta = -dpsi/dr = -(dpsi/ds) / (dr/ds)
%
%   theta-derivative: second-order central difference, periodic.
%   s-derivative:     second-order central difference in the interior and
%                     second-order one-sided differences at the two radial
%                     boundaries (there is no node beyond them).
%
%   Outputs are (Ntheta x Nr). Vx, Vy are the Cartesian components.

    [Ntheta, Nr] = size(psi);
    dtheta = 2*pi / Ntheta;
    [~, ds, r, drds] = radial_mesh(Nr, Rmin, Rmax, alpha);

    % Radial velocity
    dpsi_dtheta = (circshift(psi, -1, 1) - circshift(psi, 1, 1)) / (2*dtheta);
    Vr = dpsi_dtheta ./ r.';

    % Azimuthal velocity
    dpsi_ds = zeros(Ntheta, Nr);
    dpsi_ds(:, 2:Nr-1) = (psi(:, 3:Nr) - psi(:, 1:Nr-2)) / (2*ds);
    dpsi_ds(:, 1)  = (-3*psi(:, 1)  + 4*psi(:, 2)    - psi(:, 3))    / (2*ds);
    dpsi_ds(:, Nr) = ( 3*psi(:, Nr) - 4*psi(:, Nr-1) + psi(:, Nr-2)) / (2*ds);
    Vtheta = -dpsi_ds ./ drds.';

    % Polar to Cartesian
    Vx = Vr .* cos(theta) - Vtheta .* sin(theta);
    Vy = Vr .* sin(theta) + Vtheta .* cos(theta);
end