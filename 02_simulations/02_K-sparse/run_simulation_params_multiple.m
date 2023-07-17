function [u, v, time] = run_simulation_params_multiple(X_train, Y_train, X_valid, Y_valid, settings)
% Code to run simulations

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

if numel(fn) == 7

    Param1 =  settings.params.(fn{1});
    Param2 =  settings.params.(fn{2});
    Param3 =  settings.params.(fn{3});
    Param4 =  settings.params.(fn{4});
    Param5 =  settings.params.(fn{5});
    Param6 =  settings.params.(fn{6});
    Param7 =  settings.params.(fn{7});
    for param1 = Param1
        for param2 = Param2
            for param3 = Param3
                for param4 = Param4
                    for param5 = Param5
                        for param6 = Param6
                            for param7 = Param7
                                if settings.need_laplacian == true
                                    curr_struct = struct(fn{1}, param1, ...
                                                         fn{2}, param2, ...
                                                         fn{3}, param3, ...
                                                         fn{4}, param4, ...
                                                         fn{5}, param5, ...
                                                         fn{6}, param6, ...
                                                         fn{7}, param7, ...
                                                         'L1', L_graph_u, ...
                                                         'L2', L_graph_v, ...
                                                         'mode', settings.run_mode);
                                else
                                    curr_struct = struct(fn{1}, param1, ...
                                                         fn{2}, param2, ...
                                                         fn{3}, param3, ...
                                                         fn{4}, param4, ...
                                                         fn{5}, param5, ...
                                                         fn{6}, param6, ...
                                                         fn{7}, param7, ...
                                                         'mode', settings.run_mode);
                                end

                                curr_struct.track = false;
                                [u, v] = settings.function_call(X_train, Y_train, curr_struct);
                                if any(any(isnan(u))) || any(any(isnan(v))) || any(any(isinf(u))) || any(any(isinf(v)))
                                    continue
                                end
                                
                                curr_correlation = sum(additional_correlations(u, v, X_valid', Y_valid'));
                                
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
    end

    if numel(fieldnames(best_struct)) ~= 0
        
        best_struct.track = true;
        best_struct.prefix = strcat(num2str(settings.itr), '_');
        tStart = tic;
        [u, v] = settings.function_call(X_train, Y_train, best_struct);
        time = toc(tStart);
    else
        time = inf;
    end
 
elseif numel(fn) == 3

    Param1 =  settings.params.(fn{1});
    Param2 =  settings.params.(fn{2});
    Param3 =  settings.params.(fn{3});
    for param1 = Param1
        for param2 = Param2
            for param3 = Param3
                curr_struct = struct(fn{1}, param1, ...
                                    fn{2}, param2, ...
                                    fn{3}, param3); 

                curr_struct.track = false;
                        
                [u, v] = settings.function_call(X_train, Y_train, curr_struct);
                curr_correlation = sum(additional_correlations(u, v, X_valid', Y_valid'));

                if curr_correlation > best_correlation
                    best_correlation = curr_correlation;
                    best_struct = curr_struct;
                end
                idx = idx + 1;
            end
        end
    end

    best_struct.track = true;
    best_struct.prefix = strcat(num2str(settings.itr), '_');
    tStart = tic;
    [u, v] = settings.function_call(X_train, Y_train, best_struct);
    time = toc(tStart);
    
elseif numel(fn) == 1
    curr_struct = struct(fn{1}, settings.params.(fn{1}));
    tStart = tic;
    [u, v] = settings.function_call(X_train, Y_train, curr_struct);
    time = toc(tStart);
end
