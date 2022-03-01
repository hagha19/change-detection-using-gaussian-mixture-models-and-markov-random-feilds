function [X, R] = imstack2vectors(S, MASK)
[M, N, n] = size(S);
if nargin == 1% Number of function input arguments
   MASK = true(M, N);%true(m, n) or true([m, n]) is an m-by-n matrix of logical ones.
else
   MASK = MASK ~= 0;
end

[P, Q] = find(MASK);
%ind = find(X) locates all nonzero elements of array X, and returns the linear indices of those elements in vector ind.

R = [P, Q];

Q = M*N;
X = reshape(S, Q, n); %B = reshape(A,m,n) returns the m-by-n matrix B whose elements are taken column-wise from A. An error results if A does not have m*n elements.
%[M, N, n] = size(X)
MASK = reshape(MASK, Q, 1);
%[M, N, n] = size(MASK)
X = X(MASK, :);
end