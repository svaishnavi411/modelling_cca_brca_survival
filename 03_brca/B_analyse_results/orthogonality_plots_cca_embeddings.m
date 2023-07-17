clear all;
%  close all;
addpath('/Users/vs5/Documents/repos/cca_brca_prediction/final/3_cca/methods/cca')
addpath('/Users/vs5/Documents/repos/cca_brca_prediction/final/3_cca/methods/group_scca')
addpath('/Users/vs5/Documents/repos/cca_brca_prediction/final/3_cca/methods/data_preprocessing')
addpath('/Users/vs5/Documents/repos/cca_brca_prediction/final/3_cca/methods/gn_scca')
addpath('/Users/vs5/Documents/repos/cca_brca_prediction/final/3_cca/methods/graph_cca_proximal')
addpath('/Users/vs5/Documents/repos/cca_brca_prediction/final/3_cca/methods/deflation')

num_experiments = 10000;
T = table('Size', [num_experiments, 10], ...
        'VariableTypes', {'int16', 'int16', ... 
                          'double', 'int16', ...
                          'double', 'double', 'double', 'string', 'string', 'double'}, ...
        'VariableNames', {'n', 'p', ...
                          'noise_sigma', 'itr', ...
                          'dim', 'add_var', 'sum_var', 'cca_scheme', 'method', 'elapsed_time'});


cca_scheme_names = ["cca", "sparse", "gnscca", "gcca"];
print_names = ["CCA", "SCCA", "GN-SCCA", "GN-SCCA-PG"];
deflation_names = ["hd", "pd", "opd"];
print_deflation = ["HD", "PD", "OPD"];
data_loc = '../data/splits_journal/';
suffix = '_log_histogram';
num_genes = 1000;
fold_num = 3;

fig = figure();

for cca_scheme = 2:3

    %% Load the simulation data

    X_test = table2array(readtable(strcat(data_loc, 'test', suffix, '/', num2str(num_genes), "/", num2str(fold_num), "/", 'genes.csv')));
    Y_test = table2array(readtable(strcat(data_loc, 'test', suffix, '/', num2str(num_genes), "/", num2str(fold_num),"/", 'image.csv')));

        load((strcat('../results/', ...
                     cca_scheme_names(cca_scheme), '_', deflation_names(1), '_1000_0_log_histogram_best.mat')))

        values_cca = zeros(50, 50, 4, 3, 2);

        Cxy = X_test'*Y_test;
        for j=1:50
            Cxy = hotelling_deflation(Cxy, u(:, j), v(:, j));
            values_cca(:, j, cca_scheme, 1, 1) = vecnorm(Cxy'*u);
            values_cca(:, j, cca_scheme, 1, 2) = vecnorm(Cxy *v);

        end
        values_cca(:, :, cca_scheme, 1, 1) = triu(values_cca(:, :, cca_scheme, 1, 1)); 
        values_cca(:, :, cca_scheme, 1, 2) = triu(values_cca(:, :, cca_scheme, 1, 2));

        load((strcat('../results/', ...
                     cca_scheme_names(cca_scheme), '_', deflation_names(2), '_1000_0_log_histogram_best.mat')))

        Cxy = X_test'*Y_test;
        X_PD = X_test;
        Y_PD = Y_test;
        for j=1:50

            [XY_PD, Y_PD, Cxy] = projection_deflation(X_PD, Y_PD, u(:, j), v(:, j));
            values_cca(:, j, cca_scheme, 2, 1) = vecnorm(Cxy'*u);
            values_cca(:, j, cca_scheme, 2, 2) = vecnorm(Cxy *v);
        end

        values_cca(:, :, cca_scheme, 2, 1) = triu(values_cca(:, :, cca_scheme, 2, 1)); 
        values_cca(:, :, cca_scheme, 2, 2) = triu(values_cca(:, :, cca_scheme, 2, 2));

        % load((strcat('../results/', ...
        %              cca_scheme_names(cca_scheme), '_', deflation_names(3), '_1000_0_log_histogram_best.mat')))

        % Cxy = X_test'*Y_test;
        % X_OPD = X_test;
        % Y_OPD = Y_test;

        % R = [];
        % S = [];

        % for j=1:50

        %     [X_OPD, Y_OPD, Cxy, R, S] = orthogonal_projection_deflation(X_OPD, Y_OPD, u(:, j), v(:, j), R, S);
        %     values_cca(:, j, cca_scheme, 3, 1) = vecnorm(Cxy'*u);
        %     values_cca(:, j, cca_scheme, 3, 2) = vecnorm(Cxy *v);
        % end

        values_cca(:, :, cca_scheme, 3, 1) = triu(values_cca(:, :, cca_scheme, 3, 1)); 
        values_cca(:, :, cca_scheme, 3, 2) = triu(values_cca(:, :, cca_scheme, 3, 2)); 

    max_val =  [0, 0];
    values_cca = mean(values_cca, 5);
    new_max = max(values_cca, [], 'all');

    values_cca(values_cca == 0) = NaN;

    for deflation=1:3
        ax = subplot(3, 3, (deflation-1)*3 + (cca_scheme-1));
        s = pcolor(values_cca(:, :, cca_scheme, deflation));
        set(s, 'EdgeColor', 'none');
        ax.CLim = [1 new_max];
        colormap(ax,jet(10));
        ax.XLabel.String = print_names(cca_scheme);
        ax.YLabel.String = print_deflation(deflation);
        set(gca, 'FontSize', 15);
    end
end


cb = colorbar;
cb.Position = [0.925, 0.1, 0.025, 0.825];
cb.Ticks = linspace(0, 2.5*10^8, 6) ; 

saveas(fig, strcat('../plots/', cca_scheme_names(cca_scheme), '_journal.png'))

% close;
