function [psi, r, theta, D] = solve_laplace_cylinder(Nr, Ntheta, Rmin, Rmax, Uinf, C0, Cinf, alpha)
%SOLVE_LAPLACE_CYLINDER  Stream function for uniform flow past a circular cylinder.
%
%   Solves the Laplace equation for the stream function psi on the annulus
%   Rmin <= r <= Rmax using second-order finite differences and a sparse
%   direct solve.
%
%   Governing equation, transformed to the uniform coordinates (s, theta):
%       (1/r_s^2) psi_ss + [1/(r*r_s) - r_ss/r_s^3] psi_s + (1/r^2) psi_thth = 0
%   where r_s and r_ss are the first and second derivatives of r(s).
%
%   Boundary conditions:
%       r = Rmin :  psi = C0                              (cylinder surface)
%       r = Rmax :  psi = Uinf*Rmax*sin(theta) + Cinf     (free-stream flow)
%       theta    :  periodic
%
%   Inputs:  Nr, Ntheta  grid points in r and theta
%            Rmin, Rmax  inner and outer radius
%            Uinf        free-stream speed
%            C0, Cinf    constant offsets of psi at the inner/outer boundary
%            alpha       radial stretching parameter (see radial_mesh)
%
%   Outputs: psi    (Ntheta x Nr) stream function
%            r      (Nr x 1) radial coordinates
%            theta  (Ntheta x 1) azimuthal coordinates, 0 <= theta < 2*pi
%            D      sparse system matrix (use spy(D) to inspect its structure)

    [~, ds, r, drds, d2rds2] = radial_mesh(Nr, Rmin, Rmax, alpha);

    theta  = 2*pi * (0:Ntheta-1).' / Ntheta;
    dtheta = 2*pi / Ntheta;

    % Unknowns are ordered with theta varying fastest: psi(i,j) -> (j-1)*Ntheta + i
    Ntotal = Nr * Ntheta;
    node   = @(i, j) (j-1)*Ntheta + i;

    % Preallocate sparse triplets (at most 5 entries per row)
    rows = zeros(5*Ntotal, 1);
    cols = zeros(5*Ntotal, 1);
    vals = zeros(5*Ntotal, 1);
    f    = zeros(Ntotal, 1);
    nnzCount = 0;

    for j = 1:Nr
        for i = 1:Ntheta
            row = node(i, j);

            if j == 1
                % Cylinder surface: Dirichlet condition
                nnzCount = nnzCount + 1;
                rows(nnzCount) = row;  cols(nnzCount) = row;  vals(nnzCount) = 1;
                f(row) = C0;

            elseif j == Nr
                % Far field: Dirichlet condition for uniform flow
                nnzCount = nnzCount + 1;
                rows(nnzCount) = row;  cols(nnzCount) = row;  vals(nnzCount) = 1;
                f(row) = Uinf * Rmax * sin(theta(i)) + Cinf;

            else
                % Interior node: five-point stencil from central differences
                metric1 = 1 / drds(j)^2;
                metric2 = 1 / (r(j)*drds(j)) - d2rds2(j) / drds(j)^3;

                cLower  = metric1/ds^2 - metric2/(2*ds);     % psi(i, j-1)
                cUpper  = metric1/ds^2 + metric2/(2*ds);     % psi(i, j+1)
                cAzim   = 1 / (r(j)^2 * dtheta^2);           % psi(i-1, j) and psi(i+1, j)
                cCenter = -2*metric1/ds^2 - 2*cAzim;         % psi(i, j)

                % Periodic wrap in theta
                iPrev = mod(i-2, Ntheta) + 1;
                iNext = mod(i,   Ntheta) + 1;

                k = nnzCount + (1:5);
                rows(k) = row;
                cols(k) = [node(i,j-1), node(iPrev,j), node(i,j), node(iNext,j), node(i,j+1)];
                vals(k) = [cLower, cAzim, cCenter, cAzim, cUpper];
                nnzCount = nnzCount + 5;
            end
        end
    end

    D = sparse(rows(1:nnzCount), cols(1:nnzCount), vals(1:nnzCount), Ntotal, Ntotal);
    psi = reshape(D \ f, Ntheta, Nr);
end