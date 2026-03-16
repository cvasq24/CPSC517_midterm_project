 matrices = {'msc01050', 'bcsstk38', 'bcsstk36', 's3dkt3m2', 'apache2', 'audikw_1'}; 
%  matrices = {'s3dkt3m2'}; 
colNames = {'Matrix', 'N_Rows', 'NNZ_Original', ...
            'NNZ_Natural', 'Time_Ord_Natural', 'Time_Chol_Natural', ...
            'NNZ_RCM', 'Time_Ord_RCM', 'Time_Chol_RCM', ...
            'NNZ_Exact_MD', 'Time_Ord_Exact_MD', 'Time_Chol_Exact_MD', ...
            'NNZ_Approx_AMD', 'Time_Ord_Approx_AMD', 'Time_Chol_Approx_AMD', ...
            'NNZ_MATLAB_AMD', 'Time_Ord_MATLAB_AMD', 'Time_Chol_MATLAB_AMD'};

for m = 1:length(matrices)
    load(['data/', matrices{m}, '.mat']);
    A = Problem.A;
    n = size(A,1);
    nnz_A = nnz(A);
    fprintf('\n=== %s (n=%d, nnz=%d) ===\n', matrices{m}, n, nnz_A);
    
    % try with different orderings
    % natural, RCM, our exact MD, our AMD, and MATLAB's AMD
    tic; F_nat = chol(A, 'lower'); t_chol_nat = toc;
    nnz_nat = nnz(F_nat);
    
    tic; P_rcm = symrcm(A); t_ord_rcm = toc;
    tic; F_rcm = chol(A(P_rcm,P_rcm), 'lower'); t_chol_rcm = toc;
    nnz_rcm = nnz(F_rcm);
    
    tic; P_exa = min_degree(A); t_ord_exa = toc;
    tic; F_exa = chol(A(P_exa,P_exa), 'lower'); t_chol_exa = toc;
    nnz_exa = nnz(F_exa);

    tic; P_app = ap_min_degree(A); t_ord_app = toc;
    tic; F_app = chol(A(P_app,P_app), 'lower'); t_chol_app = toc;
    nnz_app = nnz(F_app);
    
    tic; P_mat = symamd(A); t_ord_mat = toc;
    tic; F_mat = chol(A(P_mat,P_mat), 'lower'); t_chol_mat = toc;
    nnz_mat = nnz(F_mat);
    
    % save results
    rowData = {matrices{m}, n, nnz_A, ...
        nnz_nat, 0, t_chol_nat, ...
        nnz_rcm, t_ord_rcm, t_chol_rcm, ...
        nnz_exa, t_ord_exa, t_chol_exa, ...
        nnz_app, t_ord_app, t_chol_app, ...
        nnz_mat, t_ord_mat, t_chol_mat};
    
    resultsTable = cell2table(rowData, 'VariableNames', colNames);
    csv_filename = sprintf('results/%s_benchmark_results.csv', matrices{m});
    writetable(resultsTable, csv_filename);
    fprintf('Results saved to %s\n', csv_filename);

    % save figures
    fig = figure('Name', matrices{m}, 'Position', [100, 100, 1200, 800]);    
    
    subplot(2,3,1); spy(F_nat); title(sprintf('Natural\nnz = %d', nnz_nat));
    subplot(2,3,2); spy(F_rcm); title(sprintf('RCM\nnz = %d', nnz_rcm));
    
    % Handle the Exact MD Subplot conditionally
    subplot(2,3,3); spy(F_exa); title(sprintf('Exact MD\nnz = %d', nnz_exa));
    subplot(2,3,4); spy(F_app); title(sprintf('Approx AMD\nnz = %d', nnz_app));
    subplot(2,3,5); spy(F_mat); title(sprintf('MATLAB AMD\nnz = %d', nnz_mat));
    subplot(2,3,6); spy(A); title(sprintf('Original A\nnz = %d', nnz_A));
    
    sgtitle(sprintf('%s (n=%d)', matrices{m}, n), 'Interpreter', 'none', 'FontSize', 16);
    
    img_filename = sprintf('results/%s_spy_plots.png', matrices{m});
    exportgraphics(fig, img_filename, 'Resolution', 300);
    close(fig); 
    fprintf('Saved spy plots to %s\n', img_filename);
end
disp('All evaluations complete.');