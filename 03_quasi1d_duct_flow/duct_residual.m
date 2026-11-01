function R = duct_residual(q, M1, par)
%DUCT_RESIDUAL  Residuals of the discretized duct-flow equations.
%
%   Steady, adiabatic, constant-area flow of an ideal gas with wall friction.
%   P, rho, and T are normalized by their inlet values, and velocity by
%   sqrt(R*T_inlet), so the gas law reads P = rho*T and the speed is
%   u = sqrt(k) * M * sqrt(T). Length is measured in hydraulic diameters.
%
%   Equations (one set per node 2..N, or per interval 1..N-1 for momentum):
%     continuity :  rho*M*sqrt(T) = M1               (rho*u constant)
%     momentum   :  P_i - P_{i+1} = G*(u_{i+1} - u_i) + (f*dx/4)*G*(u_i + u_{i+1})
%                   with G = rho*u = sqrt(k)*M1. Using rho*u = G in the
%                   convective and friction terms makes the interval equation
%                   algebraic; friction is evaluated with the trapezoidal rule.
%     energy     :  T*(1 + (k-1)/2*M^2) = 1 + (k-1)/2*M1^2      (T0 constant)
%     state      :  P = rho*T
%
%   par fields used: N, k, LoverD, f (Darcy friction factor).

    k  = par.k;
    N  = par.N;
    dx = par.LoverD / (N - 1);
    beta = par.f * dx / 4;

    [M, P, rho, T] = unpack_duct_state(q, M1);
    V  = M .* sqrt(T);

    iL = 1:N-1;            % left node of each interval
    iR = 2:N;              % right node of each interval

    R_mass   = rho(iR) .* M(iR) .* sqrt(T(iR)) - M1;
    R_mom    = k*M1*(V(iL) - V(iR)) + (P(iL) - P(iR)) - beta*k*M1*(V(iL) + V(iR));
    R_energy = T(iR) .* (1 + 0.5*(k-1)*M(iR).^2) - (1 + 0.5*(k-1)*M1^2);
    R_state  = P(iR) - rho(iR) .* T(iR);

    R = [R_mass; R_mom; R_energy; R_state];
end