clear all;
close all;

addpath('../../01_methods/data_preprocessing')

rng(42)

n = 100;
p = 200;
q = 200;
d = 5;

for sparsity = [0.25, 0.5, 0.75]
    for sigma = [0.1, 0.25]
        for itr = 1:10
            tic

            %% Generate the simulated data
            [X_train, Y_train, X_valid, Y_valid, X_test, Y_test, U, V, w_train, w_valid, w_test] = simulate_sparse_multiple( ...
                                    n, p, q, d, sparsity, sigma, 0.15, 0.25);

            %% Save the simulation data
            save(strcat('../data/sparse/', num2str(sparsity*100), '_', ...
                        num2str(n), '_', num2str(p), '_', num2str(q), '_', num2str(d), ...
                        '_', num2str(sigma*100), '_', num2str(itr), '.mat'));
        end
    end
end
