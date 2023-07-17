function [u, v, ecorr] = cca(paras)

%  Cxx = X'* X;
%  Cyy = Y'* Y;
%  Cxy = X'* Y;

[U, S, V]  = svd(paras.Cxy);

u = U(:, 1);
v = V(:, 1);
%  u = Cxx^(-0.5)*U(:, 1);
%  v = Cyy^(-0.5)*V(:, 1);

u = u / norm(u);
v = v / norm(v);

% estimated correlation coefficient
ecorr = S(1,1);
