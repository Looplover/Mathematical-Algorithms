function [x] = doolittle(A,b)
[L,U] = LU_decomp(A);
b_new = forwardsub(L,b);
x = backwardsub(U,b_new);