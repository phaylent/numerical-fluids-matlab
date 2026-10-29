function [s, ds, r, drds, d2rds2] = radial_mesh(Nr, Rmin, Rmax, alpha)
%RADIAL_MESH  Exponentially stretched radial grid and its metric terms.
%
%   Maps a uniform computational coordinate s in [0,1] to the physical radius
%       r(s) = Rmin + (Rmax - Rmin) * (exp(alpha*s) - 1) / (exp(alpha) - 1)
%   so grid points cluster near the cylinder (s = 0) and spread out toward the
%   far field (s = 1). Larger alpha gives stronger clustering.
%
%   Outputs (all Nr x 1): s, r, dr/ds, d2r/ds2. ds is the uniform spacing in s.

    assert(alpha > 0, 'alpha must be positive.');

    s  = linspace(0, 1, Nr).';
    ds = s(2) - s(1);

    scale  = (Rmax - Rmin) / (exp(alpha) - 1);
    r      = Rmin + scale * (exp(alpha*s) - 1);
    drds   = scale * alpha   * exp(alpha*s);
    d2rds2 = scale * alpha^2 * exp(alpha*s);
end