function [sigmad,sigmadh,sigmadv] = poolw(in,w,r,c,p) 

%
% [sigmad,sigmadh,sigmadv] = poolw(in,w,r,c,p) 
%
% w - weights
% r - number of rows
% c - number of columns
% p - number of variables

% Allan Aasbjerg Nielsen
% aa@space.dtu.dk

r1 = r-1;
c1 = c-1;
w = reshape(w, r, c);
w = w(1:r1,1:c1);
for i = 1:p
    A = reshape(in(:,i), r, c);
    Dh = A(1:r1, 1:c1) - A(1:r1, 2:c);
    H(:,i) = reshape(Dh, r1*c1, 1);
    Dv = A(1:r1, 1:c1) - A(2:r, 1:c1);
    V(:,i) = reshape(Dv, r1*c1, 1);
end

sigmadh = covw(H,w);
sigmadv = covw(V,w);

% pool
sigmad = c1/(r1+c1)*sigmadh + r1/(r1+c1)*sigmadv;

return

% simple pool
sigmad = 0.5*(sigmadh+sigmadv);
