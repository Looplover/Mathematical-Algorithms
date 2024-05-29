function [x] = forwardsub(A,b)
n = length(b);
x = zeros(1,n);
x(1) = b(1)/A(1,1);
for i = 2:n
    x(i) = b(i);
    for j = 1:i-1
        x(i) = x(i) - A(i,j)*x(j);
    end
    x(i) = x(i)/A(i,i);
end
end
