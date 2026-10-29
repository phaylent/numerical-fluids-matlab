function D = spatial_operator(scheme, N)
%SPATIAL_OPERATOR  Periodic first-derivative matrix, scaled by the grid spacing.
%
%   Returns D such that  dx * du/dx ~= D*u  on a uniform periodic grid of N
%   nodes. Multiplying by CFL = a*dt/dx later gives the semi-discrete system
%   du/dt = -(a/dx) * D * u.
%
%   scheme:
%     'CD2'    second-order central difference
%                  (u_{i+1} - u_{i-1}) / (2 dx)
%     'CMPT4'  fourth-order tridiagonal compact scheme
%                  (1/4) u'_{i-1} + u'_i + (1/4) u'_{i+1} = (3/4)(u_{i+1} - u_{i-1}) / dx
%              The derivative is defined implicitly, so applying it requires
%              solving a tridiagonal (periodic) system: D = A \ B.

    I  = eye(N);
    Ip = circshift(I,  1, 2);     % (Ip*u)_i = u_{i+1}, periodic
    Im = circshift(I, -1, 2);     % (Im*u)_i = u_{i-1}, periodic

    switch upper(scheme)
        case 'CD2'
            D = 0.5 * (Ip - Im);

        case 'CMPT4'
            A = I + 0.25 * (Ip + Im);
            B = 0.75 * (Ip - Im);
            D = A \ B;

        otherwise
            error('Unknown spatial scheme "%s". Use CD2 or CMPT4.', scheme);
    end
end