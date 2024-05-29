disp("Q1");
A1 = diag([1 2 3]);
b1 = [1 1 1];
disp(diagsolve(A1,b1));

disp("Q2");
A2 = [1 0 0;1 1 0;3 0.5 1];
b2 = [1 2 1];
disp(forwardsub(A2,b2));

disp("Q3");
A3 = [1 -1 3; 0 2 -3; 0 0 -6.5];
b3 = [1 7 6.5];
disp(backwardsub(A3,b3));

disp("Q5");
A5 = [4 1 -1;5 1 2; 6 1 1];
b5 = [-2 4 6];
disp(gausselim(A5,b5));

disp("Q6");
A6 = [10 -7 0;-3 2.099 6;5 -1 5];
n = length(A6);
b6 = zeros(1,n);
[x,Ass] = gausselimsimp(A6,b6);
disp(Ass);
answer = 1;
for j = 1:n
    answer = answer*Ass(j,j);
end
disp(answer);

