function [q, info] = newton_solve(resFun, q0, tol, maxIter)
%NEWTON_SOLVE  Newton's method with a finite-difference Jacobian.
%
%   Solves resFun(q) = 0. Iteration stops when the largest change in any
%   unknown falls below tol. info reports the iteration count, final step
%   size, final residual, and whether the tolerance was met.

    q = q0;
    stepNorm = inf;
    iter = 0;

    while stepNorm > tol && iter < maxIter
        R  = resFun(q);
        J  = fd_jacobian(resFun, q, R);
        dq = -(J \ R);
        q  = q + dq;
        stepNorm = max(abs(dq));
        iter = iter + 1;
    end

    info.iterations = iter;
    info.stepNorm   = stepNorm;
    info.resNorm    = max(abs(resFun(q)));
    info.converged  = stepNorm <= tol;

    if ~info.converged
        warning('newton_solve:noConvergence', ...
            'Newton iteration did not converge in %d iterations (step %.2e).', ...
            maxIter, stepNorm);
    end
end

function J = fd_jacobian(resFun, q, R0)
%FD_JACOBIAN  Forward-difference Jacobian, one column per unknown.
    n = numel(q);
    J = zeros(n, n);
    for j = 1:n
        h  = 1e-7 * max(1, abs(q(j)));
        qp = q;
        qp(j) = qp(j) + h;
        J(:, j) = (resFun(qp) - R0) / h;
    end
end