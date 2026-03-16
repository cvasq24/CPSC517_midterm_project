% ablations.m
dataset = 's3dkt3m2';
scriptPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(scriptPath, '..', 'data', [dataset, '.mat']);
load(dataPath);

A = Problem.A;

testNames = {'No_PostOrdering', 'No_Dense', 'No_Absorption'};
nTests = length(testNames);

times = zeros(nTests, 1);
nnz_L = zeros(nTests, 1); 

fprintf('--- Running Ablations on %s ---\n', dataset);

fig_sub = figure('Position', [100, 100, 1000, 800], 'Visible', 'off');
for t = 1:nTests
    fprintf('Running %s... ', testNames{t});
    
    % Call the 4 distinct files you created
    tic;
    if t == 1
        P = ap_min_degree_nopost(A);
    elseif t == 2
        P = ap_min_degree_nodense(A);
    elseif t == 3
        P = ap_min_degree_noabsorb(A);
    end
    times(t) = toc;
    
    % Quick Cholesky to get the fill-in (NNZ) for Slide 9
    F = chol(A(P,P), 'lower');
    nnz_L(t) = nnz(F);

    subplot(2,2,t+1);
    spy(F);
    title(strrep(testNames{t}, '_', ' '));
    
    fprintf('Done! Time: %.4f s | NNZ: %d\n', times(t), nnz_L(t));
end
saveas(fig_sub, fullfile(scriptPath, 'ablation_comparison_spy.png'));
close(fig_sub);

% Print table to the console and save it to a CSV
ablationTable = table(testNames', times, nnz_L, ...
    'VariableNames', {'Configuration', 'Ordering_Time_Sec', 'Factor_NNZ'});
disp(ablationTable);
writetable(ablationTable, 'ablation_times.csv');
disp('Saved to /ablation_times.csv');