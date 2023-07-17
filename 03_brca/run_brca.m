function run_brca(num_genes, data_loc, settings, L_graph)
% Code to run the specified method on BRCA dataset with the specific number
% of genes.

disp('Num genes:')
disp(num_genes)

cd ../
for idx = 1:numel(settings.add_path)
    addpath(genpath([pwd settings.add_path(idx)]));
end
cd run

for fold_num = ["0", "1", "2", "3", "4"]
    % load data sets with current num_genes
    X_train = table2array(readtable(strcat(data_loc, 'train/', num2str(num_genes), "/", fold_num, '/genes.csv')));
    Y_train = table2array(readtable(strcat(data_loc, 'train/', num2str(num_genes), "/", fold_num,'/image.csv')));

    % NOTE THE CHANGE BETWEEN TEST AND VALID!!!
    
    X_valid = table2array(readtable(strcat(data_loc, 'test/', num2str(num_genes), "/", fold_num, '/genes.csv')));
    Y_valid = table2array(readtable(strcat(data_loc, 'test/', num2str(num_genes), "/", fold_num,'/image.csv')));

    X_test = table2array(readtable(strcat(data_loc, 'valid/', num2str(num_genes), "/", fold_num, '/genes.csv')));
    Y_test = table2array(readtable(strcat(data_loc, 'valid/', num2str(num_genes), "/", fold_num,'/image.csv')));

    % Run normalization on data
    X_valid = getNormalization(X_valid, X_train);
    Y_valid = getNormalization(Y_valid, Y_train);
    
    % Run normalization on data
    X_test = getNormalization(X_test, X_train);
    Y_test = getNormalization(Y_test, Y_train);

     % Run normalization on data
    X_train = getNormalization(X_train, X_train);
    Y_train = getNormalization(Y_train, Y_train);
    
    [N, p] = size(X_train)
    [N1, q] = size(Y_train)
        
    % Get number of experiments to run
    number_of_experiments = 1;

    fn = fieldnames(settings.params);
    for k=1:numel(fn)
        number_of_experiments = number_of_experiments * numel(settings.params.(fn{k}));
    end

    disp('Number of experiments to run: ')
    disp(number_of_experiments)

    idx = 1;
    best_correlation = -1;
    best_struct = struct();

    fn = fieldnames(settings.params);

if numel(fn) == 6

    Param1 =  settings.params.(fn{1});
    Param2 =  settings.params.(fn{2});
    Param3 =  settings.params.(fn{3});
    Param4 =  settings.params.(fn{4});
    Param5 =  settings.params.(fn{5});
    Param6 =  settings.params.(fn{6});
    for param1 = Param1
        for param2 = Param2
            for param3 = Param3
                for param4 = Param4
                    for param5 = Param5
                        for param6 = Param6
                            if settings.need_laplacian == true
                                curr_struct = struct(fn{1}, param1, ...
                                                     fn{2}, param2, ...
                                                     fn{3}, param3, ...
                                                     fn{4}, param4, ...
                                                     fn{5}, param5, ...
                                                     fn{6}, param6, ...
                                                     'L1', L_graph, ...
                                                     'mode', settings.run_mode);
                            else
                                curr_struct = struct(fn{1}, param1, ...
                                                     fn{2}, param2, ...
                                                     fn{3}, param3, ...
                                                     fn{4}, param4, ...
                                                     fn{5}, param5, ...
                                                     fn{6}, param6, ...
                                                     'mode', settings.run_mode);
                            end
                            curr_struct.track = false;
                            [u, v] = settings.function_call(X_train, Y_train, curr_struct);
                            if any(any(isnan(u))) || any(any(isnan(v))) || any(any(isinf(u))) || any(any(isinf(v)))
                                continue
                            end
                            
                            x_variate = u'*X_valid';
                            y_variate = v'*Y_valid';
                            curr_correlation = x_variate*y_variate'/(sqrt(x_variate*x_variate')*sqrt(y_variate*y_variate'));

                            if curr_correlation > best_correlation
                                best_correlation = curr_correlation;
                                best_struct = curr_struct;
                            end
                            idx = idx + 1;
                        end
                    end
                end
            end
        end
    end

    if numel(fieldnames(best_struct)) ~= 0
        
        best_struct.track = true;
        best_struct.prefix = strcat(num2str(settings.itr), '_');
        [u, v] = settings.function_call(X_train, Y_train, best_struct);
    end
 
elseif numel(fn) == 2

    Param1 =  settings.params.(fn{1});
    Param2 =  settings.params.(fn{2});
    for param1 = Param1
        for param2 = Param2
            curr_struct = struct(fn{1}, param1, ...
                                fn{2}, param2);

            curr_struct.track = false;

            [u, v] = settings.function_call(X_train, Y_train, curr_struct);

            x_variate = u'*X_valid';
            y_variate = v'*Y_valid';
            curr_correlation = x_variate*y_variate'/(sqrt(x_variate*x_variate')*sqrt(y_variate*y_variate'));

            if curr_correlation > best_correlation
                best_correlation = curr_correlation;
                best_struct = curr_struct;
            end
            idx = idx + 1;
        end
    end

    best_struct.track = true;
    best_struct.prefix = strcat(num2str(settings.itr), '_');
    [u, v] = settings.function_call(X_train, Y_train, best_struct);
end   

fileID = fopen(strcat(settings.save_loc, ...
                    settings.save_name, num2str(num_genes), '_', fold_num, '.txt'),'w');
fprintf(fileID,'%6s \n','correlations');

x_variate = u'*X_train';
y_variate = v'*Y_train';
train_correlation = x_variate*y_variate'/(sqrt(x_variate*x_variate')*sqrt(y_variate*y_variate'));

disp('On training data')
disp(train_correlation)
fprintf(fileID,'%6.5f \n',train_correlation);

x_variate = u'*X_valid';
y_variate = v'*Y_valid';
valid_correlation = x_variate*y_variate'/(sqrt(x_variate*x_variate')*sqrt(y_variate*y_variate'));

disp('On validation data')
disp(valid_correlation)
fprintf(fileID,'%6.5f \n', valid_correlation);

x_variate = u'*X_test';
y_variate = v'*Y_test';
test_correlation = x_variate*y_variate'/(sqrt(x_variate*x_variate')*sqrt(y_variate*y_variate'));

fprintf(fileID,'%6.5f \n', test_correlation);
fclose(fileID);

disp('On testing data')
disp(test_correlation)


save(strcat(settings.save_loc, ...
    settings.save_name, num2str(num_genes),"_", fold_num, '_best.mat'), ...
    'train_correlation', 'valid_correlation', 'test_correlation', 'u', 'v', 'best_struct')

end
