function [psiExact, A, B] = exact_stream_function(r, theta, Rmin, Rmax, Uinf)
%EXACT_STREAM_FUNCTION  Analytical solution on the finite annulus Rmin <= r <= Rmax.
%
%   Separation of variables gives psi = (A*r + B/r) * sin(theta). Imposing
%   psi = 0 at r = Rmin and psi = Uinf*Rmax*sin(theta) at r = Rmax gives
%       A = Uinf*Rmax^2 / (Rmax^2 - Rmin^2),   B = -A*Rmin^2.
%
%   This is the exact solution of the problem the numerical solver actually
%   solves. It approaches the infinite-domain result psi = Uinf*(r - Rmin^2/r)*sin(theta)
%   only as Rmax -> infinity. Valid for C0 = Cinf = 0.

    A = Uinf * Rmax^2 / (Rmax^2 - Rmin^2);
    B = -A * Rmin^2;
    psiExact = (A*r.' + B./r.') .* sin(theta);     % (Ntheta x Nr)
end
