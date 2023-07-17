function [X_new, Y_new, XY_new, R_new, S_new] = orthogonal_projection_deflation(X, Y, u, v, R, S)

% Input sizes: 
% X - nxp, Y - nxq, u - px1, v - qx1
% R and S keep track of previous vectors for orthogonalization

p = size(X, 2);
q = size(Y, 2);

if numel(R) == 0
    R = zeros(p, 1);
    S = zeros(q, 1);
end

r_new = (eye(p) - R*R')*u;
r_new = r_new / norm(r_new);

s_new = (eye(q) - S*S')*v;
s_new = s_new / norm(s_new);

X_new = (eye(p) - r_new*r_new')*X';
X_new = X_new';

Y_new = (eye(q) - s_new*s_new')*Y';
Y_new = Y_new';

XY_new = X_new'*Y_new;

R_new = [R, r_new];
S_new = [S, s_new];
