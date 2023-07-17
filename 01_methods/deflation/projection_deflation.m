function [X_new, Y_new, XY_new] = projection_deflation(X, Y, u, v)

% Input sizes: 
% X - nxp, Y - nxq, u - px1, v - qx1

p = size(X, 2);
q = size(Y, 2);

X_new = (eye(p) - u*u')*X';
X_new = X_new';

Y_new = (eye(q) - v*v')*Y';
Y_new = Y_new';

XY_new = X_new'*Y_new;
