function [M1, sol, nIter] = solve_inlet_mach(PN, par, M1Guess, q0)
%SOLVE_INLET_MACH  Find the inlet Mach number that gives a prescribed exit pressure.
%
%   Newton iteration on M1 for r(M1) = P_N(M1) - PN = 0, where P_N is the
%   exit pressure divided by the inlet pressure and dP_N/dM1 comes from a
%   finite difference between two flow solves. Each solve starts from the
%   previous one, so only a few inner Newton steps are needed.
%
%   Inputs:  PN       target exit pressure ratio
%            par      parameter struct (see run_duct_flow)
%            M1Guess  starting inlet Mach number
%            q0       (optional) starting unknown vector
%
%   Outputs: M1 (converged), sol (flow solution at M1), nIter (outer iterations)

    if nargin < 4, q0 = []; end

    M1max = fanno_choke_mach(par.k, par.f * par.LoverD);
    dM    = 1e-6;
    M1    = M1Guess;
    q     = q0;
    nIter = 0;
    converged = false;

    while ~converged && nIter < par.maxIter
        if M1 <= 0 || M1 >= M1max
            error('solve_inlet_mach:outOfRange', ...
                'Iterate M1 = %.4f left the unchoked range (0, %.4f).', M1, M1max);
        end

        solA = solve_duct(M1,      par, q);
        solB = solve_duct(M1 + dM, par, solA.q);

        r    = solA.P(end) - PN;
        drdM = (solB.P(end) - solA.P(end)) / dM;
        step = -r / drdM;

        M1    = M1 + step;
        q     = solA.q;
        nIter = nIter + 1;
        converged = abs(step) < par.tolM1;
    end

    if ~converged
        warning('solve_inlet_mach:noConvergence', ...
            'Inlet-Mach iteration did not converge for PN = %.3f.', PN);
    end

    sol = solve_duct(M1, par, q);
end