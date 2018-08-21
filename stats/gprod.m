function [N,c] = gprod(N1,N2,varargin)
% gprod: product of multivariate normal pdfs
%
%     [N,c] = gprod(N1,N2,N3,...Nk) with inputs N1 = MVN(mu1,Sigma1), ..., Nk = MVN(muk,Sigmak)
%             returns N = MVN(mu,Sigma) together with the constant c such that
%             c * MVN(x|mu,Sigma) = MVN(x|mu1,Sigma1) * ... * MVN(x|muk,Sigmak).
%
%
% Tobias Winner, 20.08.2018
% email: winner.tobias@gmail.com


assert(N1.dim==N2.dim, 'Dimensions must match')
d = N1.dim;

S = inv(N1.Lambda + N2.Lambda);
m = S*(N1.eta + N2.eta);
N = MVN(m,S);

c = (2*pi)^(-d/2) * det(S)^(1/2) * det(N1.Sigma)^(-1/2) * det(N2.Sigma)^(-1/2) * ...
    exp( -(N1.mu'*N1.eta + N2.mu'*N2.eta - m'*inv(S)*m)/2 );

if nargin > 2
    [N,c(2)] = gprod(N,varargin{:});
    c = prod(c);
end

end
