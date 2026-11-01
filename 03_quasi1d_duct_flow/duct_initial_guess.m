function q0 = duct_initial_guess(M1, par)
%DUCT_INITIAL_GUESS  Linear-ramp starting guess for the Newton solver.
%
%   For subsonic Fanno flow, M rises along the duct while P, rho, and T fall.
%   The exit Mach number is guessed as min(1.5*M1, 0.95), the exit
%   temperature follows from conservation of stagnation temperature, and the
%   exit density from conservation of mass flux. All variables ramp linearly
%   from the inlet to those exit values.

    k = par.k;
    N = par.N;

    T0      = 1 + 0.5*(k-1)*M1^2;
    Mexit   = min(1.5*M1, 0.95);
    Texit   = T0 / (1 + 0.5*(k-1)*Mexit^2);
    rhoExit = M1 / (Mexit * sqrt(Texit));
    Pexit   = rhoExit * Texit;

    frac = (1:N-1).' / (N-1);
    q0 = [M1 + frac*(Mexit   - M1);
          1  + frac*(Pexit   - 1);
          1  + frac*(rhoExit - 1);
          1  + frac*(Texit   - 1)];
end