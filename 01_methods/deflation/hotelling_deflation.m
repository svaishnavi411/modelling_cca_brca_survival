function [XY_new] = hotelling_deflation(XY, u, v)

% Input sizes: 
% X - nxp, Y - nxq, u - px1, v - qx1
%  
%  X_new = X;
%  Y_new = Y;

XY_new = XY - (u*u')*XY*(v*v');
