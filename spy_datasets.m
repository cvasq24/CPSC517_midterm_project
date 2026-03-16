matrices = {'msc01050', 'bcsstk38', 'bcsstk36', 's3dkt3m2', 'apache2', 'audikw_1'}; 

fig = figure('Name', 'Sparsity Patterns', 'Position', [100 100 1400 500]);
t = tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'compact');
for m = 1:length(matrices)
    load(['data/', matrices{m}, '.mat']);
    A = Problem.A;

    subplot(2,3,m);
    spy(A,1);
    title(sprintf('%s  (n=%d)', matrices{m}, size(A,1)), ...
      'Interpreter', 'none', 'FontSize', 9);
end
exportgraphics(fig, 'results/spy_plots.png', 'Resolution', 500);