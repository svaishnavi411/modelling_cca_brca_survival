function [X_train, Y_train, X_valid, Y_valid, X_test, Y_test, U, V, w_train, w_valid, w_test] = simulate_basic_multiple(n, p, q, d,...
    sigma_factor, valid_ratio, test_ratio)

%%% X_train is p x num_train, X_valid is p x num_valid, X_test is p x num_test
%%% Y_train is q x num_train, Y_valid is q x num_valid, Y_test is q x num_test

%% Step 1: Simulating data

num_valid = ceil(valid_ratio * n);
num_test = ceil(test_ratio * n);
num_train = n - num_test - num_valid;

U = randn(p, d);
disp('U generated')
V = randn(q, d);
disp('V generated')

w_train = randn(num_train, d);
w_valid = randn(num_valid, d);
w_test = randn(num_test, d);

X_train = U * w_train' + sigma_factor * randn(p, num_train);
Y_train = V * w_train' + sigma_factor * randn(q, num_train);

X_valid = U * w_valid' + sigma_factor * randn(p, num_valid);
Y_valid = V * w_valid' + sigma_factor * randn(q, num_valid);

X_test = U * w_test' + sigma_factor * randn(p, num_test);
Y_test = V * w_test' + sigma_factor * randn(q, num_test);

end

