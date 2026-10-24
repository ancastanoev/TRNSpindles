% Spindle-proximity diagrams in the (rho, R) plane for three connectivity
% settings (a) wRR = 0.15, (b) wRR = 0.25, (c) wTR = 0.6.
%
% Each diagram shows the TRN limit-cycle amplitude band coloured by proximity
% to the biological spindle target window, as the mean delay rho varies.
%
% Output: spindle_proximity_panel.png

clear;
clc;
close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

output_file = fullfile(script_dir, 'spindle_proximity_panel.png');

files = { ...
    fullfile(script_dir, 'rho_LC_wRT_0_3_wPP_1_1_wRR_0_15_wTR_0_5.mat'), ...
    fullfile(script_dir, 'rho_LC_wRT_0_3_wPP_1_1_wRR_0_25_wTR_0_5.mat'), ...
    fullfile(script_dir, 'rho_LC_wRT_0_25_wPP_1_1_wRR_0_2_wTR_0_6.mat') ...
};

tile_titles = { ...
    '(a)\ \ $w_{RR}=0.15$', ...
    '(b)\ \ $w_{RR}=0.25$', ...
    '(c)\ \ $w_{TR}=0.6,\ w_{RR}=0.2$' ...
};

target_rho      = [0.30, 0.60];
target_R_top    = [10, 15];
target_R_bottom = [0, 4];

n_vars = 8;

cmap = flipud(parula(256));

panel_data = struct([]);
global_R_min = inf;
global_R_max = -inf;

for k = 1:numel(files)

    fname = files{k};

    panel_data(k).exists = false;
    panel_data(k).fname  = fname;

    if ~isfile(fname)
        warning('File not found: %s. Tile %d will be marked missing.', fname, k);
        continue;
    end

    fprintf('Loading %s\n', fname);

    previous_warning_state = warning;
    warning('off', 'all');
    data = load(fname);
    warning(previous_warning_state);

    if ~isfield(data, 'x')
        warning('Variable x is missing from file: %s', fname);
        continue;
    end

    all_x = data.x;
    num_cycles = size(all_x, 2);

    mesh_value = (size(all_x, 1) - 2) / n_vars;
    if mesh_value ~= round(mesh_value)
        error(['Invalid dimensions in file:\n%s\n' ...
               '(size(x,1)-2)/n_vars is not an integer.'], fname);
    end

    n_mesh = round(mesh_value);

    rho_vals  = zeros(1, num_cycles);
    R_tops    = zeros(1, num_cycles);
    R_bottoms = zeros(1, num_cycles);
    scores    = zeros(1, num_cycles);
    invalid   = false(1, num_cycles);

    for i = 1:num_cycles
        cyc = reshape(all_x(1:end-2, i), n_vars, n_mesh);

        T_trace = cyc(1, :);
        R_trace = cyc(2, :);
        rho_val = all_x(end, i);

        R_top    = max(R_trace);
        R_bottom = min(R_trace);

        rho_vals(i)  = rho_val;
        R_tops(i)    = R_top;
        R_bottoms(i) = R_bottom;

        if any(~isfinite(T_trace)) || any(~isfinite(R_trace)) || ~isfinite(rho_val)
            invalid(i) = true;
        elseif any(T_trace < 0) || any(R_trace < 0)
            invalid(i) = true;
        end

        % Distance from the biological spindle target window
        drho = max([0, target_rho(1) - rho_val, rho_val - target_rho(2)]);
        dTop = max([0, target_R_top(1) - R_top, R_top - target_R_top(2)]);
        dBot = max([0, target_R_bottom(1) - R_bottom, R_bottom - target_R_bottom(2)]);

        scores(i) = sqrt(drho^2 + dTop^2 + dBot^2);
    end

    max_score = max(scores);

    if max_score == 0
        norm_scores = ones(1, num_cycles);
    else
        norm_scores = 1 - scores / max_score;   
        % 1 = biological spindle, 0 = escape/farthest
    end

    panel_data(k).exists       = true;
    panel_data(k).rho_vals     = rho_vals;
    panel_data(k).R_tops       = R_tops;
    panel_data(k).R_bottoms    = R_bottoms;
    panel_data(k).scores       = scores;
    panel_data(k).norm_scores  = norm_scores;
    panel_data(k).invalid      = invalid;

    valid_rho     = rho_vals(isfinite(rho_vals));
    valid_bottoms = R_bottoms(isfinite(R_bottoms));
    valid_tops    = R_tops(isfinite(R_tops));

    if ~isempty(valid_rho)
        panel_data(k).xlim = [min(valid_rho) - 0.1, max(valid_rho) + 0.1];
    else
        panel_data(k).xlim = [0 1];
    end

    if ~isempty(valid_bottoms)
        global_R_min = min(global_R_min, min(valid_bottoms));
    end

    if ~isempty(valid_tops)
        global_R_max = max(global_R_max, max(valid_tops));
    end
end

if ~isfinite(global_R_min) || ~isfinite(global_R_max)
    error('No valid data found in the selected panels.');
end

shared_ylim = [global_R_min - 1.0, global_R_max + 1.0];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1 1 13.5 4.6], ...
    'PaperPositionMode', 'auto');

panel_gap    = 0.018;
left_margin  = 0.065;
top_margin   = 0.10;
bot_margin   = 0.15;
cb_gap       = 0.020;
cb_width     = 0.018;
right_margin = 0.105;

x_cb        = 1 - right_margin - cb_width;
panels_area = (x_cb - cb_gap) - left_margin;
panel_w     = (panels_area - 2 * panel_gap) / 3;
plot_h      = 1 - top_margin - bot_margin;

panel_pos = @(k) [left_margin + (k-1) * (panel_w + panel_gap), ...
                  bot_margin, panel_w, plot_h];

colormap(fig, cmap);

ax_handles  = gobjects(1, numel(files));
colorbar_ax = [];

for k = 1:numel(files)

    ax = axes('Parent', fig, 'Position', panel_pos(k));
    ax_handles(k) = ax;

    hold(ax, 'on');

    set(ax, ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'TickDir', 'out', ...
        'LineWidth', 1.3, ...
        'FontSize', 11, ...
        'FontName', 'Times New Roman');

    caxis(ax, [0 1]);

    if ~panel_data(k).exists
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'missing .mat', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 11);
        continue;
    end

    colorbar_ax = ax;

    rho_vals     = panel_data(k).rho_vals;
    R_tops       = panel_data(k).R_tops;
    R_bottoms    = panel_data(k).R_bottoms;
    norm_scores  = panel_data(k).norm_scores;
    invalid      = panel_data(k).invalid;

    plot(ax, rho_vals, R_tops,    'k-', 'LineWidth', 1.2);
    plot(ax, rho_vals, R_bottoms, 'k-', 'LineWidth', 1.2);

    blocks = contiguous_blocks(~invalid);

    for b = 1:size(blocks, 1)

        idx = blocks(b, 1):blocks(b, 2);

        if numel(idx) < 2
            continue;
        end

        X = [rho_vals(idx); rho_vals(idx)];
        Y = [R_bottoms(idx); R_tops(idx)];
        C = [norm_scores(idx); norm_scores(idx)];

        surf(ax, X, Y, zeros(size(X)), C, ...
            'EdgeColor', 'none', ...
            'FaceColor', 'interp', ...
            'FaceAlpha', 0.95);
    end

    view(ax, 2);

    if k == 2
        xlabel(ax, 'Control Parameter, $\rho$', ...
            'Interpreter', 'latex', ...
            'FontSize', 13);
    end

    if k == 1
        ylabel(ax, 'TRN Rate $R$ [Hz]', ...
            'Interpreter', 'latex', ...
            'FontSize', 13);
    end

    title(ax, tile_titles{k}, ...
        'Interpreter', 'latex', ...
        'FontSize', 13);

    xlim(ax, panel_data(k).xlim);
    ylim(ax, shared_ylim);

    if k > 1
        ax.YTickLabel = [];
    end

    grid(ax, 'off');
end

if isempty(colorbar_ax) || ~isgraphics(colorbar_ax, 'axes')
    error('No valid diagram was loaded, so the shared colorbar could not be created.');
end

cb = colorbar(colorbar_ax);

for k = 1:numel(files)
    if isgraphics(ax_handles(k), 'axes')
        ax_handles(k).Position = panel_pos(k);
    end
end

cb.Ticks = [0 1];
cb.TickLabels = {'Escape State', 'Biological Spindle'};
cb.FontSize = 10;
cb.FontName = 'Times New Roman';
cb.Position = [x_cb, bot_margin, cb_width, plot_h];

exportgraphics(fig, output_file, 'Resolution', 300);

fprintf('Saved %s\n', output_file);

function blocks = contiguous_blocks(mask)
    mask = mask(:)';
    d = diff([false, mask, false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    blocks = [starts(:), ends(:)];
end