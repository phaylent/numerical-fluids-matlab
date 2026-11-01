function [M, P, rho, T] = unpack_duct_state(q, M1)
%UNPACK_DUCT_STATE  Split the unknown vector into full profiles.
%
%   The unknowns are stored variable by variable for nodes 2..N:
%       q = [M(2:N); P(2:N); rho(2:N); T(2:N)]
%   The inlet node is fixed (M1, 1, 1, 1), since P, rho, and T are normalized
%   by their inlet values, and is prepended here. Outputs are N x 1 columns.

    n   = numel(q) / 4;          % n = N - 1
    M   = [M1; q(1:n)];
    P   = [1;  q(n+1:2*n)];
    rho = [1;  q(2*n+1:3*n)];
    T   = [1;  q(3*n+1:4*n)];
end