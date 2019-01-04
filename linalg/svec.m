function [x] = svec(X)
% svec(X) returns the symmetric vectorization 
%         x = [X(1,1) sqrt(2)*X(1,2) ... sqrt(2)*X(1,n) X(2,2) ... sqrt(2)*X(2,n) X(n,n)]'
%         for symmetrix matrices X

assert(all(all(X'==X)),'X must be symmetric')

X = diag(diag(X)) + tril(X,-1) * sqrt(2);
x = X(tril(true(size(X))));

end
