function [U, V]=SparseCCAGroupMultiple(X, Y, lambda_x, lambda_y, T_x, w_x, T_y, w_y, K, deflation)
% Group structured sparse CCA  
% Input:
%   X, Y: data
%   lambda_x, lambda_y: regularization parameter for X and Y (c_1 or theta
%   in Eq.(3) in SSCCA paper)
%   T_x, w_x, group structure and weights for X side (if they are empty,
%   then simply use L1-constraint)
%   T_y, w_y  group structure and weights for Y side (if they are empty,
%   then simply use L1-constraint)
%   para: parameter
%   For simplicity, we assume the parameter tau is 1
% Output
%   coefficients u for X side
%   coefficients v for Y side
%   cor: correlation

%% Set up parameters

    disp(K)
    
    if exist('para', 'var') && isfield(para, 'verbose')
        verb=para.verbose;  % verbose for printing
    else
        verb=true;
    end
    
    if exist('para', 'var') && isfield(para, 'tol')
        tol=para.tol;  % tol for algorithm
    else
        tol=1e-4;
    end
    
    if exist('para', 'var') && isfield(para, 'maxiter')
        maxiter= para.maxiter;
    else
        maxiter=50;
    end
    
    if exist('para', 'var') && isfield(para, 'threshold')
        threshold= para.threshold;
    else
        threshold=1e-5;
    end
    
    if exist('para', 'var') && isfield(para, 'inner_verbose')
        inner_verb=para.inner_verbose;  % verbose for printing
    elsep
        inner_verb=false;
    end
    
    if exist('para', 'var') && isfield(para, 'inner_tol')
        inner_tol=para.inner_tol;  % tol for algorithm
    else
        inner_tol=1e-6;
    end
    
    if exist('para', 'var') && isfield(para, 'inner_maxiter')
        inner_maxiter= para.inner_maxiter;
    else
        inner_maxiter=1000;
    end
    
%     inner_para=struct('verbose', inner_verb, 'tol', inner_tol, 'maxiter', inner_maxiter);

    
    if isempty(T_x) || isempty(w_x)
        group_x_tag=false;  % group on x side
    else
        group_x_tag=true;
    end
    
    if isempty(T_y) || isempty(w_y)
        group_y_tag=false;  % group on y side
    else
        group_y_tag=true;
    end   
%     
%     if (verb && group_x_tag && group_y_tag)
%         disp('Sparse CCA with group structured penalties on X and Y');
%     elseif (verb && group_x_tag)
%         disp('Sparse CCA with group structured penalty on X and L1 on Y');
%     elseif (verb && group_y_tag)
%         disp('Sparse CCA with group structured penalty on Y and L1 on X');
%     else
%         disp('Sparse CCA with L1 Penalty');
%     end
    
    %% Alternating Procedure
    
    p = size(X, 2);
    q = size(Y, 2);
    
    U = zeros(p, K);
    V = zeros(q, K);
%     F = [];
    
    XY = X'*Y;
%     XY_orig = XY;
    YX = XY';
      
    R = [];
    S = [];
    
    for k = 1:K
        % initialization
        col_X=size(X,2);
        col_Y=size(Y,2);    
    
        u_old=projl2(randn(col_X,1));
        v_old=projl2(randn(col_Y,1));    
        
        iter=0;
        while (iter<maxiter)

            if (group_x_tag)
                u=group_exgap(XY*v_old, lambda_x, T_x, w_x, inner_para);
            else
                u=maxL1L2(XY*v_old, lambda_x);
            end
            u_gap=max(abs(u-u_old));
            u_old=u;

            if (group_y_tag)
                v=group_exgap(YX*u_old, lambda_y, T_y, w_y, inner_para);                
            else
                v=maxL1L2(YX*u_old, lambda_y);
            end

            v_gap=max(abs(v-v_old));
            v_old=v;
            
            if (verb)
                fprintf('Iter=%d, u_gap=%f, v_gap=%f\n', iter, u_gap, v_gap);
            end

            if (u_gap<tol && v_gap<tol)
                break
            end
            
            iter=iter+1;
          
        end       
        u(abs(u)<threshold)=0;
        v(abs(v)<threshold)=0;
        
        if (any(u))
            u=u/sqrt(sum(u.^2));        
        end
        if (any(v))
            v=v/sqrt(sum(v.^2));
        end
        
        u = u./norm(u);
        v = v./norm(v);
        
        U(:, k) = u;
        V(:, k) = v;
        
        if deflation == "hd"
            [XY] = hotelling_deflation(XY, u, v);
        
        elseif deflation == "pd"
            [X, Y, XY] = projection_deflation(X, Y, u, v);
        
        elseif deflation == "opd"
            [X, Y, XY, R, S] = orthogonal_projection_deflation(X, Y, u, v, R, S);
        end 
        
%         factor = sum(sum((XY .*(u*v')))) / sum(sum(((u*v') .*(u*v'))));
%         XY = (XY - factor*(u*v'));
%         XY = XY/norm(XY, "fro");
%         disp(norm(XY, "fro"));
%         norm_uv =  norm(u*v', "fro");
%         factor = sum(sum((XY .*(u*v')))) /norm_uv;  %  sum(sum(((u*v') .*(u*v'))));
%         XY = XY - factor*(u*v') /norm_uv;  %/norm(XY, "fro");
%         F = [F, factor];
    
    end
    
% disp(norm(XY_orig - U*diag(F)*V', "fro"));
% disp(norm(XY_orig, "fro"))
end
