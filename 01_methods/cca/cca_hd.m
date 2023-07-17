function [U, V] = cca_hd(X, Y, paras)

    paras.hd = true;
            
    paras.Cxy = X'* Y;

    for k = 1:paras.dim

        [u, v, ecorr] = cca(paras);
    
        u = u/norm(u);
        v = v/norm(v);

        U(:, k) = u;
        V(:, k) = v;

        [XY_new] = hotelling_deflation(paras.Cxy, u, v);
        paras.Cxy = XY_new;
    end
end
