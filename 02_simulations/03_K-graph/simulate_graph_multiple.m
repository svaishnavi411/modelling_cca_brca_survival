function [X_train, Y_train, X_valid, Y_valid, X_test, Y_test, L1, L2, G1, G2, U, V, w_train, w_valid, w_test]  = simulate_graph_multiple(n, p, q, d, sigma_factor, valid_factor, test_factor)

%% Step 1: Simulating data

num_valid = ceil(valid_factor * n);
num_test = ceil(test_factor * n);
num_train = n - num_test - num_valid;


L1 = [];
G1 = [];

% Constructing a fully connected component here. 

prob_edge = 1.25*log(p)/p;

while 1
  A1 = rand(p, p) < prob_edge;
  A1 = triu(A1,1);
  A1 = A1 + A1';

  G1 = graph(A1);
  L1 = laplacian(G1);

  if max(G1.conncomp) == 1
      break;
  end
end

[E1, ~] = eig(full(L1));

U = [];

for k=1:d
    u = (rand(1, k*5)* E1(:, 2:(k*5+1))')';
    if k > 1
        u = (eye(p) - U*U')*u;
    end
    u = u ./ norm(u);
    U = [U, u];
end
disp('U generated')

L2 = [];
G2 = [];
prob_edge = 1.25*log(q)/q;

while 1
    A2 = rand(q, q) < prob_edge;
    A2 = triu(A2,1);
    A2 = A2 + A2';
    G2 = graph(A2);
    L2 = laplacian(G2);
    
    if max(G2.conncomp) == 1
        break;
    end
end

[E2, ~] = eig(full(L2));

V = [];

for k=1:d
    v = (rand(1, k*5)* E2(:, 2:(k*5+1))')';
    if k > 1
        v = (eye(q) - V*V')*v;
    end
    v = v ./ norm(v);
    V = [V, v];
end

disp('V generated')

%% Noise independent of u and v
w_train = randn(num_train, d);
w_valid = randn(num_valid, d);
w_test = randn(num_test, d);

X_train = U * w_train' + sigma_factor * randn(p, num_train);
Y_train = V * w_train' + sigma_factor * randn(q, num_train);

X_valid = U * w_valid' + sigma_factor * randn(p, num_valid);
Y_valid = V * w_valid' + sigma_factor * randn(q, num_valid);

X_test = U * w_test' + sigma_factor * randn(p, num_test);
Y_test = V * w_test' + sigma_factor * randn(q, num_test);
            
% Run normalization on data
X_valid = getNormalization(X_valid', X_train')';
Y_valid = getNormalization(Y_valid', Y_train')';

% Run normalization on data
X_test = getNormalization(X_test', X_train')';
Y_test = getNormalization(Y_test', Y_train')';

 % Run normalization on data
X_train = getNormalization(X_train', X_train')';
Y_train = getNormalization(Y_train', Y_train')';
          
end

