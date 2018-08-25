classdef MVN < matlab.mixin.Copyable
    % MVN: Multivariate Normal distribution
    %
    %     N = MVN(mu,Sigma) defines a Gaussian distribution with mean mu and
    %     covariance Sigma.
    %
    %
    %     ---------- Parameters ----------
    %
    %     N.mu and N.Sigma return the moment parameters mu and Sigma, respectively.
    %     N.c returns the normalization constant for this parameterization.
    %
    %     N.eta and N.Lambda return the canonical parameters eta and Lambda, respectively.
    %     N.a returns the normalization constant for this parameterization.
    %
    %
    %     ---------- Sampling and Evaluation ----------
    %
    %     N.sample(N) draws N samples x(:,1), x(:,2), ..., x(:,N) from
    %     MVN(mu,Sigma).
    %
    %     N.eval(x) evaluates the pdf at x(:,1), x(:,2), ..., x(:,N).
    %
    %
    % Tobias Winner, 20.08.2018
    % email: winner.tobias@gmail.com
    
    properties ( SetAccess=private )
        dim
    end
    
    properties
        mean
        cov
    end
    
    properties ( Dependent, Hidden )
        mu, Sigma
    end
    
    properties ( Hidden, SetAccess=private )
        c
        eta, Lambda, a
        V, D
    end
    
    methods
        
        function self = MVN(mu, Sigma, varargin)
            assert(iscolumn(mu), 'mu must be a column vector or scalar')
            self.dim = length(mu);
            
            assert(all((size(Sigma)==[self.dim self.dim])), 'Input dimensions must match')
            assert(rank(Sigma) == self.dim, 'Sigma must have full rank')
            
            [V,D] = eig(Sigma);
            assert(all(diag(D) > 0), 'Sigma must be positive')
            
            self.cov = Sigma;
            self.mean = mu;
            self.c = 1/(sqrt(det(2*pi*Sigma)));
            self.V = V;
            self.D = D;
        end
        
        function mu = get.mu(self)
            mu = self.mean;
        end
        
        function Sigma = get.Sigma(self)
            Sigma = self.cov;
        end
        
        function Lambda = get.Lambda(self)
            if isempty(self.Lambda)
                self.Lambda = inv(self.cov);
            end
            Lambda = self.Lambda;
        end
        
        function eta = get.eta(self)
            if isempty(self.eta)
                self.eta = self.Lambda*self.mu;
            end
            eta = self.eta;
        end
        
        function a = get.a(self)
            if isempty(self.a)
                self.a = -(self.dim * log(2*pi) - log(det(self.Lambda)) + self.eta'*inv(self.Lambda)*self.eta )/2;
            end
            a = self.a;
        end
        
        function V = get.V(self)
            if isempty(self.V)
                [V,D] = eig(self.cov);
                self.V = V;
                self.D = D;
            else
                V = self.V;
            end
        end
        
        function D = get.D(self)
            if isempty(self.D)
                [V,D] = eig(self.cov);
                self.V = V;
                self.D = D;
            else
                D = self.D;
            end
        end
        
        function self = set.mean(self, mu)
            assert(iscolumn(mu), 'mu must be a column vector or scalar')
            assert(length(mu)==self.dim, 'Dimensions must match')
            
            self.mean = mu;
            
            self.eta = [];
            self.a = [];
        end
        
        function self = set.cov(self, Sigma)
            assert(all((size(Sigma)==[self.dim self.dim])), 'Input dimensions must match')
            assert(rank(Sigma) == self.dim, 'Sigma must have full rank')
            
            [V,D] = eig(Sigma);
            assert(all(diag(D) > 0), 'Sigma must be positive')
            
            self.cov = Sigma;
            self.c = 1/(sqrt(det(2*pi*Sigma)));
            
            self.Lambda = [];
            self.eta = [];
            self.a = [];
            self.V = V;
            self.D = D;
        end
        
        function m = max(self)
            m = self.c;
        end
        
        function H = entropy(self)
            H = -log(self.c);
        end
        
        function p = eval(self,x)
            p = mvnpdf(x',self.mean',self.cov)';
        end
        
        function x = sample(self, N)
            x = mvnrnd(self.mean',self.cov,N)';
        end
        
    end
    
end
