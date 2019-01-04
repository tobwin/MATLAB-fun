function x = vech(X)
% vech(X) returns the vectorization x = [X(1,1) X(1,2) ... X(1,n) X(2,2) X(2,1) ... X(n,n)]' 
%         for symmetric matrices X.

assert(all(all(X==X')), 'X must be symmetric.')

x = X(triu(true(size(X))));

end
