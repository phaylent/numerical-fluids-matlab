function M1max = fanno_choke_mach(k, fL)
%FANNO_CHOKE_MACH  Largest inlet Mach number a duct can pass without choking.
%
%   For total friction parameter fL = f*L/D, the exit reaches M = 1 exactly
%   when fanno_parameter(M1) = fL. Any larger inlet Mach number is not possible.

    M1max = fzero(@(m) fanno_parameter(m, k) - fL, [0.05, 1 - 1e-9], ...
                  optimset('TolX', 1e-14));
end