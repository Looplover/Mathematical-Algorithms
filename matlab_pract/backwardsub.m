function [x] = backwardsub(A,b)
n = length(b);
x = zeros(1,n);
x(n) = b(n)/A(n,n);
for i = n-1:-1:1
    x(i) = b(i);
    for j = n:-1:i+1
        x(i) = x(i) - A(i,j)*x(j);
    end
    x(i) = x(i)/A(i,i);
end
end