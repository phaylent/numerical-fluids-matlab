function G = build_propagator(scheme, D, CFL)
%BUILD_PROPAGATOR  One-step update matrix for du/dt = -(a/dx) D u.
%
%   Because the problem is linear with constant coefficients, each time scheme
%   reduces to u^{n+1} = G*u^n with a fixed matrix G.
%
%   scheme:
%     'Trapezoidal'  (I + CFL/2 D) u^{n+1} = (I - CFL/2 D) u^n
%                    implicit, second-order in time
%     'RK4'          classical fourth-order Runge-Kutta. For a linear system
%                    this equals the degree-4 Taylor polynomial of exp(-CFL*D):
%                    G = I + Z + Z^2/2 + Z^3/6 + Z^4/24,  Z = -CFL*D

    N = size(D, 1);
    I = eye(N);

    switch lower(scheme)
        case 'trapezoidal'
            G = (I + 0.5*CFL*D) \ (I - 0.5*CFL*D);

        case 'rk4'
            Z = -CFL * D;
            G = I;
            term = I;
            for m = 1:4
                term = term * Z / m;
                G = G + term;
            end

        otherwise
            error('Unknown time scheme "%s". Use Trapezoidal or RK4.', scheme);
    end
end