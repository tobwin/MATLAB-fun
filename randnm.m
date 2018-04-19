function X = randnm(M, V, U, vargargin)
%randnm Random matrix
%   randnm(M, V, U) returns a random normal matrix ~NM(X|M,V,U)
%   randnm(M, V, U, 'chol') accepts lower-triangular U and V and samples
%   from ~NM(X|M,VV',UU')
%
% Tobias Winner, 19.04.2018
% email: winner.tobias@gmail.com

n = size(V, 1);
p = size(U, 1);
if all(all(M==0))
    M = zeros(n, p);
end
if nargin > 3
    if ismember('chol',varargin)
        vecX = kron(V, U) * randn(n * p, 1);
    end
else
    vecX = chol(kron(V, U),'lower') * randn(n * p, 1);
end
X = reshape(vecX, n, p) + M;
end
