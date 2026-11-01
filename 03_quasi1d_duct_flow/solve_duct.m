function sol = solve_duct(M1, par, q0)
%SOLVE_DUCT  Solve the duct-flow problem for a given inlet Mach number.
%
%   sol = solve_duct(M1, par)        uses the linear-ramp initial guess
%   sol = solve_duct(M1, par, q0)    starts Newton from q0 (e.g. a nearby solution)
%
%   Returned fields: x (in units of D_h), M, P, rho, T (N x 1, P/rho/T
%   normalized by inlet values), q (unknown vector), M1, info (Newton report).

    if nargin < 3 || isempty(q0)
        q0 = duct_initial_guess(M1, par);
    end

    resFun = @(q) duct_residual(q, M1, par);
    [q, info] = newton_solve(resFun, q0, par.tol, par.maxIter);

    [M, P, rho, T] = unpack_duct_state(q, M1);

    sol.x    = linspace(0, par.LoverD, par.N).';
    sol.M    = M;
    sol.P    = P;
    sol.rho  = rho;
    sol.T    = T;
    sol.q    = q;
    sol.M1   = M1;
    sol.info = info;
end