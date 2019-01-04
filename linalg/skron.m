function M = skron(A,B)
% skron(A,B) returns the symmetric Kronecker product for two n x n matrices A and B.
%
% Credit goes to David Goodmanson whose original code can be found at 
%   https://nl.mathworks.com/matlabcentral/answers/361452-symmetric-kronecker-product-in-matlab

n = size(A,1);
U = eye(n^2);
a = reshape(1:n^2,n,n);
b = a';
U = U + U(b(:),:);
c = tril(a);
c = c(:);
c(c==0) = [];
U = U(c,:);
U(U==1) = 1/sqrt(2);
U(U==2) = 1;
M = (1/2)*U*(kron(A,B)+kron(B,A))*U';

end
