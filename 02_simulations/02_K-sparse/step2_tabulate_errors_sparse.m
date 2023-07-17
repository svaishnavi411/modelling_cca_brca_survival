clear all;
close all;

rng(42)

addpath('../../01_methods/cca')
addpath('../../01_methods/group_scca')
addpath('../../01_methods/gn_scca')
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


%% First generate data for a specific setting

save_name = '../results/simulation_results_sparse.txt';

table_iter = 1;

n = 100;
p = 200;
q = 200;
d = 5;

for sparsity = [0.25] 
    for sigma = [0.1] 
        for itr = 1:10
            tic

            %% Load the simulation data
            load(strcat('../data/sparse/', num2str(sparsity*100), '_', ...
                        num2str(n), '_', num2str(p), '_', num2str(q), '_', num2str(d), ...
                        '_', num2str(sigma*100), '_', num2str(itr), '.mat'));
                
            X_train = X_train';
            Y_train = Y_train';
            X_valid = X_valid';
            Y_valid = Y_valid';
            X_test = X_test';
            Y_test = Y_test';

            %%% Run the different deflation methods
            
            %% Sparse CCA-HD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            sparse_cca_params = struct("lambda1", 0.1:sqrt(p)/10:sqrt(p), ...
                "lambda2", 0.1:sqrt(q)/10:sqrt(q), "dim", d);

            add_paths = ["/methods/group_scca/"];
            settings = struct('add_path', add_paths, ...
                'need_laplacian', false, ...
                'params', sparse_cca_params, ...
                'function_call', @sparse_cca_hd, ...
                'itr', 0);

            [U_SCCA_HD, V_SCCA_HD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
            correlations = additional_correlations(U_SCCA_HD, V_SCCA_HD, X_test', Y_test');
            sum_correlations = cumsum(correlations);
            
            for dim = 1:d
                temp_struct = [];
                temp_struct.n = n;
                temp_struct.p = p;
                temp_struct.noise_sigma = sigma;
                temp_struct.itr = itr;
                temp_struct.dim = dim;
                temp_struct.add_var = correlations(dim);
                temp_struct.sum_var = sum_correlations(dim);
                temp_struct.deflation = "HD";
                temp_struct.method = 'SCCA';
                temp_struct.time = elapsed_time;
                Table_all(table_iter, :) = struct2table(temp_struct);
                disp(temp_struct)
                table_iter = table_iter + 1;
            end
            
            disp('Sparse CCA')

            %% GN-SCCA-HD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            params_struct_1 = struct("alpha1", [10^(-1), 1, 10], ...
                "alpha2", [10^(-1), 1, 10], ...
                "lambda1", [10^(-1), 1, 10], ... % 0.2:0.2:0.6 ... %1:2:10, ...
                "lambda2", [10^(-1), 1, 10], ... % 0.2:1:2:10, ...
                "beta1", [1, 10, 100], ...
                "beta2", [1, 10, 100], ...
                "dim", d);

            add_paths = ["/methods/gn_scca/", "/methods/data_processing/"];
            settings = struct('add_path', add_paths, ...
                'need_laplacian', false, ...
                'params', params_struct_1, ...
                'function_call', @gn_scca_hd, ...
                'run_mode', "orig", ...
                'itr', 0);

            [U_GNSCCA_HD, V_GNSCCA_HD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
            correlations = additional_correlations(U_GNSCCA_HD, V_GNSCCA_HD, X_test', Y_test');
            sum_correlations = cumsum(correlations);
            
            for dim = 1:d
                temp_struct = [];
                temp_struct.n = n;
                temp_struct.p = p;
                temp_struct.noise_sigma = sigma;
                temp_struct.itr = itr;
                temp_struct.dim = dim;
                temp_struct.add_var = correlations(dim);
                temp_struct.sum_var = sum_correlations(dim);
                temp_struct.deflation = "HD";
                temp_struct.method = 'GN-SCCA';
                temp_struct.time = elapsed_time;
                Table_all(table_iter, :) = struct2table(temp_struct);
                disp(temp_struct)
                table_iter = table_iter + 1;
            end

            disp('GN-SCCA')

            %% Sparse CCA-PD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            add_paths = ["/methods/group_scca/"];
            settings = struct('add_path', add_paths, ...
                'need_laplacian', false, ...
                'params', sparse_cca_params, ...
                'function_call', @sparse_cca_pd, ...
                'itr', 0);

            [U_SCCA_PD, V_SCCA_PD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
            correlations = additional_correlations(U_SCCA_PD, V_SCCA_PD, X_test', Y_test');
            sum_correlations = cumsum(correlations);
            
            for dim = 1:d
                temp_struct = [];
                temp_struct.n = n;
                temp_struct.p = p;
                temp_struct.noise_sigma = sigma;
                temp_struct.itr = itr;
                temp_struct.dim = dim;
                temp_struct.add_var = correlations(dim);
                temp_struct.sum_var = sum_correlations(dim);
                temp_struct.deflation = "PD";
                temp_struct.method = 'SCCA';
                temp_struct.time = elapsed_time;
                Table_all(table_iter, :) = struct2table(temp_struct);
                disp(temp_struct)
                table_iter = table_iter + 1;
            end

            disp('SCCA')

            %% GN-SCCA-PD

            add_paths = ["/methods/gn_scca/", "/methods/data_processing/"];
            settings = struct('add_path', add_paths, ...
                'need_laplacian', false, ...
                'params', params_struct_1, ...
                'function_call', @gn_scca_pd, ...
                'run_mode', "orig", ...
                'itr', 0);

            [U_GNSCCA_PD, V_GNSCCA_PD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
            correlations = additional_correlations(U_GNSCCA_PD, V_GNSCCA_PD, X_test', Y_test');
            sum_correlations = cumsum(correlations);
            
            for dim = 1:d
                temp_struct = [];
                temp_struct.n = n;
                temp_struct.p = p;
                temp_struct.noise_sigma = sigma;
                temp_struct.itr = itr;
                temp_struct.dim = dim;
                temp_struct.add_var = correlations(dim);
                temp_struct.sum_var = sum_correlations(dim);
                temp_struct.deflation = "PD";
                temp_struct.method = 'GN-SCCA';
                temp_struct.time = elapsed_time;
                Table_all(table_iter, :) = struct2table(temp_struct);
                disp(temp_struct)
                table_iter = table_iter + 1;
            end

            disp('GN-SCCA')

            
            %% Sparse CCA-OPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            add_paths = ["/methods/group_scca/"];
            settings = struct('add_path', add_paths, ...
                'need_laplacian', false, ...
                'params', sparse_cca_params, ...
                'function_call', @sparse_cca_opd, ...
                'itr', 0);

            [U_SCCA_OPD, V_SCCA_OPD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
            correlations = additional_correlations(U_SCCA_OPD, V_SCCA_OPD, X_test', Y_test');
            sum_correlations = cumsum(correlations);
            
            for dim = 1:d
                temp_struct = [];
                temp_struct.n = n;
                temp_struct.p = p;
                temp_struct.noise_sigma = sigma;
                temp_struct.itr = itr;
                temp_struct.dim = dim;
                temp_struct.add_var = correlations(dim);
                temp_struct.sum_var = sum_correlations(dim);
                temp_struct.deflation = "OPD";
                temp_struct.method = 'SCCA';
                temp_struct.time = elapsed_time;
                Table_all(table_iter, :) = struct2table(temp_struct);
                disp(temp_struct)
                table_iter = table_iter + 1;
            end

            disp('SCCA')
            
            %% GN-SCCA-OPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            add_paths = ["/methods/gn_scca/", "/methods/data_processing/"];
            settings = struct('add_path', add_paths, ...
                'need_laplacian', false, ...
                'params', params_struct_1, ...
                'function_call', @gn_scca_opd, ...
                'run_mode', "orig", ...
                'itr', 0);

            [U_GNSCCA_OPD, V_GNSCCA_OPD, elapsed_time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings);
            correlations = additional_correlations(U_GNSCCA_OPD, V_GNSCCA_OPD, X_test', Y_test');
            sum_correlations = cumsum(correlations);
            
            for dim = 1:d
                temp_struct = [];
                temp_struct.n = n;
                temp_struct.p = p;
                temp_struct.noise_sigma = sigma;
                temp_struct.itr = itr;
                temp_struct.dim = dim;
                temp_struct.add_var = correlations(dim);
                temp_struct.sum_var = sum_correlations(dim);
                temp_struct.deflation = "OPD";
                temp_struct.method = 'GN-SCCA';
                temp_struct.time = elapsed_time;
                Table_all(table_iter, :) = struct2table(temp_struct);
                disp(temp_struct)
                table_iter = table_iter + 1;
            end

            disp('GN-SCCA')
            
            %% Writing
            writetable(Table_all(1:table_iter-1, :),  save_name)
            save(strcat('../results/sparse/', num2str(sparsity*100), '_', ...
                        num2str(n), '_', num2str(p), '_', num2str(q), '_', num2str(d), ...
                        '_', num2str(sigma*100), '_', num2str(itr), '.mat'))
            
            toc
        end
    end
end
