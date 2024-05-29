function [x] = cholesky(A,b)
n = length(A);
L = zeros(size(A));
for v = 1:n
    temp = 0;
    for u = 1:v-1
        temp = temp + L(v,u)*L(v,u);
    end
    L(v,v) = sqrt(A(v,v) - temp);
    for t = v+1:n
        temp = 0;
        for u = 1:v-1
            temp = temp + L(t,u)*L(v,u);
        end
        L(t,v) = (1/L(v,v))*(A(t,v)-temp);
    end
end
disp(L)
b_new = forwardsub(L,b);
x = backwardsub(L.',b_new);