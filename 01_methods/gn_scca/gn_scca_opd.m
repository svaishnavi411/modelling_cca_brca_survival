function [U, V] = gn_scca_opd(X, Y, paras)

    paras.hd = false;
    p = size(X,2);
    q = size(Y,2);

    U = zeros(p, paras.dim);
    V = zeros(q, paras.dim);
    
    R = [];
    S = [];
    
    for k = 1:paras.dim
        [u, v] = gn_scca(X, Y, paras);
        
        if any(any(isnan(u))) || any(any(isnan(v))) || any(any(isinf(u))) || any(any(isinf(v)))
            U(:, k) = u;
            V(:, k) = v;
            break
        end
        
        u = u/norm(u);
        v = v/norm(v);
        
        U(:, k) = u;
        V(:, k) = v;
        
        [X, Y, ~, R, S] = orthogonal_projection_deflation(X, Y, u, v, R, S);
    end
        
end
