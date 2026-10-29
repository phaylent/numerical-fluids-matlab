function [x, u, uExact, CFL] = solve_advection(timeScheme, spaceScheme, N, dt, tEnd, a, xLim)
%SOLVE_ADVECTION  Advance u_t + a*u_x = 0 on a periodic domain to time tEnd.
%
%   Inputs:  timeScheme   'Trapezoidal' or 'RK4'
%            spaceScheme  'CD2' or 'CMPT4'
%            N            number of grid nodes (the node at xMax is the same
%                         point as the node at xMin, so it is not repeated)
%            dt, tEnd     time step and final time (tEnd must be a multiple of dt)
%            a            wave speed
%            xLim         [xMin xMax]
%
%   Outputs: x       grid nodes (N x 1)
%            u       numerical solution at tEnd
%            uExact  exact solution at the same nodes
%            CFL     a*dt/dx

    dx  = (xLim(2) - xLim(1)) / N;
    x   = xLim(1) + dx * (0:N-1).';
    CFL = a * dt / dx;

    nSteps = round(tEnd / dt);
    if abs(nSteps*dt - tEnd) > 1e-9 * tEnd
        error('tEnd must be an integer multiple of dt.');
    end

    u = exact_advection(x, 0, a, xLim);          % initial condition

    D = spatial_operator(spaceScheme, N);
    G = build_propagator(timeScheme, D, CFL);

    for n = 1:nSteps
        u = G * u;
    end

    uExact = exact_advection(x, tEnd, a, xLim);
end