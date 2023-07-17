addpath('../../01_methods/data_preprocessing')

rng(0)

n = 100;
p = 200;
q = 200;
d = 5;
 
for sigma = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
 for itr = 1:50
     tic

     [X_train, Y_train, X_valid, Y_valid, X_test, Y_test, U, V, w_train, w_valid, w_test] = simulate_basic_single_cca_simulation( ...
                             n, p, q, d, sigma, 0.1, 0.3);

     save(strcat('../data/basic/', ...
                 num2str(n), '_', num2str(p), '_', num2str(q), '_', num2str(d), ...
                 '_', num2str(sigma*100), '_', num2str(itr), '.mat'));
     
     
 end
end
 