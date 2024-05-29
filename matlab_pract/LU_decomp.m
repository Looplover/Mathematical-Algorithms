function [L,U] = LU_decomp(A)
U = A;
L = zeros(size(A));
n = length(A);
for k = 1:n
    L(k,k) = 1;
    for j = k+1:n
        L(j,k) = U(j,k)/U(k,k);
        U(j,:) = U(j,:) - (U(j,k)/U(k,k))*U(k,:);
    end
end


