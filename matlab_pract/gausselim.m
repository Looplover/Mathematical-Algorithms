function [x] = gausselim(A,b)
augm = [A b.'];
n = length(b);
for k = 1:n
     maxind = k;
     maxval = augm(k,k);
     for h = k+1:n
         if abs(augm(h,k)) > abs(maxval)
             maxval = augm(h,k);
             maxind = h;
         end
    end
    temp = augm(k,:);
    augm(k,:) = augm(maxind,:);
    augm(maxind,:) = temp;
    for j = k+1:n
        augm(j,:) = augm(j,:) - (augm(j,k)/augm(k,k))*augm(k,:);
    end
end
disp(augm)
x = backwardsub(augm(:,1:n),augm(:,n+1));
end

