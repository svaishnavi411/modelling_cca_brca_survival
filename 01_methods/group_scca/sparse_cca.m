function [u, v] = sparse_cca(X, Y, paras)
    
[u, v, ~] = SparseCCAGroup(X, Y, paras.lambda1, paras.lambda2, [], [], [], [], []);
    
u = u / norm(u);
v = v / norm(v);
end

