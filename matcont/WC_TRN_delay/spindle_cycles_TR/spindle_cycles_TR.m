% Wilson-Cowan spindle limit cycles in the (T, R) plane, coloured by
% physiological proximity. Cut at the rho turning point.
% On-figure title and colour-bar wording live in the LaTeX caption.

clear; clc; close all;

warning('off','all');
data1 = load('spindle_cycles_part1.mat');
data2 = load('spindle_cycles_part2.mat');
warning('on','all');

% biological target windows
target_R = [10, 15];
target_T = [4, 6];
n_vars   = 8;

all_x      = [data1.x, data2.x];
num_cycles = size(all_x, 2);
n_mesh     = (size(all_x, 1) - 2) / n_vars;
if mod(size(all_x, 1) - 2, n_vars) ~= 0
    error('Data size is not compatible with n_vars = %d.', n_vars);
end


rho_vals = all_x(end, :);
drho     = diff(rho_vals);
idx_turn = find(drho < 0, 1, 'first');
if isempty(idx_turn)
    idx_turn = num_cycles;
end
valid_idx = 1:idx_turn;

% proximity score per cycle: distance of (T_peak, R_peak) from the target windows
scores = zeros(1, num_cycles);
for i = 1:num_cycles
    cycle_data = reshape(all_x(1:end-2, i), n_vars, n_mesh);
    T_peak = max(cycle_data(1, :));
    R_peak = max(cycle_data(2, :));
    dT = max([0, target_T(1) - T_peak, T_peak - target_T(2)]);
    dR = max([0, target_R(1) - R_peak, R_peak - target_R(2)]);
    scores(i) = sqrt(dT^2 + dR^2);
end
if max(scores) == 0
    norm_scores = ones(1, num_cycles);
else
    norm_scores = 1 - (scores / max(scores));
end

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 9, 6.5]);
hold on;
cmap = flipud(parula(256));

for i = valid_idx
    cycle_data = reshape(all_x(1:end-2, i), n_vars, n_mesh);
    T_trace = cycle_data(1, :);
    R_trace = cycle_data(2, :);
    c_idx = max(1, min(256, round(norm_scores(i) * 255) + 1));
    plot(T_trace, R_trace, 'Color', [cmap(c_idx, :), 0.55], 'LineWidth', 1.3);
end

xlabel('Relay Population Rate \it{T} \rm{[Hz]}', 'FontSize', 15);
ylabel('TRN Population Rate \it{R} \rm{[Hz]}', 'FontSize', 15);
set(gca, 'TickDir', 'out', 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
xlim([0 6]);
ylim([0 20.5]);
grid off;
box on;
daspect([1 2.8 1]);

colormap(cmap);
cb = colorbar('Ticks', [0, 1], ...
'TickLabels', {'Escape State', 'Biological Spindle'}, ...
'TickLabelInterpreter', 'latex');
hold off;