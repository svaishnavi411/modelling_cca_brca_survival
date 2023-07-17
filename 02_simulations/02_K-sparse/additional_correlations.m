function [corrs] = additional_correlations(U, V, X, Y)

p = size(U, 1);
q = size(V, 1);
dim = size(U, 2);

XY = X*Y';

corrs = zeros(1, dim);
R = [];
S = [];

for d = 1:dim
    if numel(R) == 0
        R = zeros(p, 1);
        S = zeros(q, 1);
    end
    
    r_new = (eye(p) - R*R')*U(:, d);
    s_new = (eye(q) - S*S')*V(:, d);
    corrs(1, d) = (r_new'*(XY)*s_new)/(norm(r_new'*X)*(norm(s_new'*Y)));  
        
    r_new = r_new / norm(r_new);
    s_new = s_new / norm(s_new);
    
    R = [R, r_new];
    S = [S, s_new];

end

end
