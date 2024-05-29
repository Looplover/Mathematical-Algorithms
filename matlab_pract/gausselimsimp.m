function [det_A,x] = gausselimsimp(A,b)
augm = [A b.'];
n = length(b);
det_A = 1;
for k = 1:n
    for j = k+1:n
        augm(j,:) = augm(j,:) - (augm(j,k)/augm(k,k))*augm(k,:);
    end
    det_A = det_A*augm(k,k);
end
x = backwardsub(augm(:,1:n),augm(:,n+1));
end

