clear all;
close all;

addpath('../../01_methods/cca')
addpath('../../01_methods/data_preprocessing')
addpath('../../01_methods/deflation')

num_experiments = 10000;
Table_all = table('Size', [num_experiments, 10], ...
        'VariableTypes', {'int16', 'int16', ... 
                          'double', 'int16', ...
                          'double', 'double', 'double', 'string', 'string', 'double'}, ...
        'VariableNames', {'n', 'p', ...
                          'noise_sigma', 'itr', ...
                          'dim', 'add_var', 'sum_var', 'deflation', 'method', 'elapsed_time'});


save_name = '../results/simulation_results_basic';

%% First generate data for a specific setting
  
table_iter = 1;

n = 100;
p = 200;
q = 200;
d = 5;

for sigma = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    for itr = 1:50
        tic

        %% Load the simulation data
        load(strcat('../data/basic/', ...
                    num2str(n), '_', num2str(p), '_', num2str(q), '_', num2str(d), ...
                    '_', num2str(sigma*100), '_', num2str(itr), '.mat'));
                    
        X_valid = getNormalization(X_valid', X_train')';
        Y_valid = getNormalization(Y_valid', Y_train')';

        X_test = getNormalization(X_test', X_train')';
        Y_test = getNormalization(Y_test', Y_train')';

        X_train = getNormalization(X_train', X_train')';
        Y_train = getNormalization(Y_train', Y_train')';
        
        X_train = X_train';
        Y_train = Y_train';
        X_valid = X_valid';
        Y_valid = Y_valid';
        X_test = X_test';
        Y_test = Y_test';

        %% CCA-HD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        cca_params = struct("dim", d);

        add_paths = ["/methods/cca/"];
        settings = struct('add_path', add_paths, ...
            'need_laplacian', false, ...
            'params', cca_params, ...
            'function_call', @cca_hd, ...
            'itr', 0);

        [U_CCA_HD, V_CCA_HD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
        correlations_curr = additional_correlations(U_CCA_HD, V_CCA_HD, X_test', Y_test');
        sum_correlations = cumsum(correlations_curr);
        
        for dim = 1:d
            temp_struct = [];
            temp_struct.n = n;
            temp_struct.p = p;
            temp_struct.noise_sigma = sigma;
            temp_struct.itr = itr;
            temp_struct.dim = dim;
            temp_struct.add_var = correlations_curr(dim);
            temp_struct.sum_var = sum_correlations(dim);
            temp_struct.deflation = "HD";
            temp_struct.method = 'CCA';
            temp_struct.time = elapsed_time;
            disp(temp_struct)
            Table_all(table_iter, :) = struct2table(temp_struct);
            table_iter = table_iter + 1;
        end
        
        disp('CCA')
        
        %% Writing
        writetable(Table_all(1:table_iter-1, :),  save_name)
        save(strcat('../results/basic/', ...
                    num2str(n), '_', num2str(p), '_', num2str(q), '_', num2str(d), ...
                    '_', num2str(sigma*100), '_', num2str(itr), '.mat'))
        
        toc
    end
end
