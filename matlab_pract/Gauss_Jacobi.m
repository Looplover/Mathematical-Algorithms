function[x,num_iter,r] = Gauss_Jacobi(A,b,x,tol)
d = diag(A);
r = (b.') - A*(x.');
num_iter = 0;
while norm(r) > tol
num_iter = num_iter + 1;
r = (b.') - A*(x.');
x = x + (r./d).';
end

