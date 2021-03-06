function [R] = rot(theta,k)
%rot Rotation matrix
%   rot(theta) returns a 2D rotation matrix
%   rot(theta,k) returns a 3D rotation matrix, rotating about the k-th axis
%   rot(theta,v) returns a 3D rotation matrix, rotating about the vector v
%
% Tobias Winner, 20.10.2017
% email: winner.tobias@gmail.com

if nargin == 1
    R = [cos(theta) -sin(theta);sin(theta) cos(theta)];
elseif nargin == 2
    if length(k) == 1
        R = eye(3);
        idx = find(~ismember(1:3,k));
        R(idx(1),idx(1)) = cos(theta);
        R(idx(2),idx(2)) = cos(theta);
        R(idx(1),idx(2)) = -sin(theta);
        R(idx(2),idx(1)) = sin(theta);
    elseif length(k) == 3
        u = k/norm(k);
        c = cos(theta); c1 = 1-c ; s = sin(theta);
        R = [c+u(1)^2*c1        u(1)*u(2)*c1-u(3)*s     u(1)*u(3)*c1+u(2)*s
            u(2)*u(1)*c1+u(3)*s c+u(2)^2*c1             u(2)*u(3)*c1-u(1)*s
            u(3)*u(1)*c1-u(2)*s u(3)*u(2)*c1+u(1)*s     c+u(3)^2*c1];
    end
end
