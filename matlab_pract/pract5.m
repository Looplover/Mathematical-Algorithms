disp("Q1")
A = [1 1 -1; 1 2 -2;-2 1 1];
[L,U] = LU_decomp(A)
disp("L")
disp(L)
disp("U")
disp(U)

disp("Q2")
disp(doolittle(A,[1 1 1]))

disp("Q3")
disp("ans using LU")
A = [16 4 4 -4; 4 10 4 2; 4 4 6 -2;-4 2 -2 4];
b = [32 26 20 -6];
disp(doolittle(A,b))

disp("Q5")
tol = 10^(-4);
x = [1/2 1/2];
A = [10 1;1 10];
b = [11 11];
[x,num_iter,r] = Gauss_Jacobi(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))

disp("Q6")
tol = 10^(-4);
x = [0 0 0];
A = [4 1 -1;2 7 1; 1 -3 12];
b = [3 19 31];
[x,num_iter,r] = Gauss_Jacobi(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))

disp("Q7")
tol = 10^(-4);
x = [0 0 0 0 0];
A = [5 -2 3 0 6;-3 9 1 -2 7.4;2 -1 -7 1 6.7;4 3 -5 7 9;2 3.5 6.1 -4 -8.1];
b = [-1 2 3 0.5 3.1];
[x,num_iter,r] = Gauss_Jacobi(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))

disp("Q8")
tol = 10^(-4);
x = [0 0 0];
A = [1 2 3;2 -1 2;3 1 -2];
b = [5 1 -1];
[x,num_iter,r] = Gauss_Jacobi(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))

disp("Q10: doing 5  and 6 with Siedel")

disp("Q5")
tol = 10^(-4);
x = [1/2 1/2];
A = [10 1;1 10];
b = [11 11];
[x,num_iter,r] = Gauss_Siedel(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))

disp("Q6")
tol = 10^(-4);
x = [0 0 0];
A = [4 1 -1;2 7 1; 1 -3 12];
b = [3 19 31];
[x,num_iter,r] = Gauss_Siedel(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))

disp("Q11")
tol = 10^(-4);
x = [0 0 0 0];
A = [2 8 3 1;0 2 -1 4;7 -2 1 2; -1 0 5 2];
b = [-2 4 3 5];
[x,num_iter,r] = Gauss_Siedel(A,b,x,tol);
disp(x)
disp("num_iter")
disp(num_iter)
disp("final error:")
disp(norm(r))


