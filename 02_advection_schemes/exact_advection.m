function u = exact_advection(x, t, a, xLim)
%EXACT_ADVECTION  Exact solution of u_t + a*u_x = 0 on a periodic domain.
%
%   The initial profile u0(x) = exp(-x^2/2) is translated by a*t. On a
%   periodic domain, whatever leaves one boundary re-enters at the other, so
%   the shifted coordinate is wrapped back into the domain with mod().
%
%   Inputs:  x     query points
%            t     time
%            a     wave speed
%            xLim  [xMin xMax] periodic domain

    L  = xLim(2) - xLim(1);
    xi = mod(x - a*t - xLim(1), L) + xLim(1);
    u  = exp(-xi.^2 / 2);
end
