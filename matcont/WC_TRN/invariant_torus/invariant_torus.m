% Invariant torus rendered as a thin tube in the (P, T, R) phase space,
% coloured by evolution time. Trajectory loaded from invariant_torus.mat.

clear; clc; close all;
load('invariant_torus.mat');
P = y(:, 1); T = y(:, 2); R = y(:, 3);
time_vec = t;

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 9, 7.5]);
hold on;

% colour map: indigo -> blue -> cyan -> green -> gold -> yellow, low t to high t
anchors = [ ...
    0.20, 0.05, 0.45;
    0.05, 0.35, 0.85;
    0.05, 0.78, 0.88;
    0.40, 0.85, 0.45;
    0.95, 0.62, 0.15;
    0.98, 0.92, 0.30];
n_colors   = 256;
pos        = linspace(0, 1, size(anchors, 1));
interp_pts = linspace(0, 1, n_colors)';
aurora_map = [interp1(pos, anchors(:,1), interp_pts), ...
              interp1(pos, anchors(:,2), interp_pts), ...
              interp1(pos, anchors(:,3), interp_pts)];

radius  = 0.011;
n_sides = 12;
[xt, yt, zt] = tube_generate(P, T, R, radius, n_sides);

% colour each tube vertex by its time value, then render with transparency
c_data = repmat(time_vec', n_sides + 1, 1);
h = surf(xt, yt, zt, c_data, ...
'EdgeColor', 'none', ...
'FaceColor', 'interp', ...
'FaceAlpha', 0.80);

colormap(aurora_map);
lighting phong;
camlight('headlight');
material([0.7 0.3 0.1 5 0.2]);

view(135, 30);
grid off;
ax = gca;
ax.LineWidth = 1.2;
ax.Box = 'on';
ax.FontName = 'Times New Roman';
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';

xlabel('$P$ component', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$T$ component', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('$R$ component', 'Interpreter', 'latex', 'FontSize', 22);

cb = colorbar('eastoutside');
cb.Label.String = 'Evolution Time ($t$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 18;
cb.TickLabelInterpreter = 'latex';

axis tight;
daspect([1 1 0.38]);
hold off;

% Build a tubular surface of radius r (n facets) around the 3-D curve (px,py,pz)
% using a parallel-transport-style local frame along the curve.
function [x, y, z] = tube_generate(px, py, pz, r, n)
    n_pts = length(px); theta = linspace(0, 2*pi, n + 1);
    x = zeros(n + 1, n_pts); y = x; z = x;
    for i = 1:n_pts
        if i < n_pts, v = [px(i+1)-px(i), py(i+1)-py(i), pz(i+1)-pz(i)];
        else, v = [px(i)-px(i-1), py(i)-py(i-1), pz(i)-pz(i-1)]; end
        v = v / (norm(v) + eps);
        if abs(v(1)) < 0.9, n1 = [1, 0, 0]; else, n1 = [0, 1, 0]; end
        u = cross(v, n1); u = u / norm(u); w = cross(v, u);
        for j = 1:n+1
            x(j,i) = px(i) + r * (cos(theta(j))*u(1) + sin(theta(j))*w(1));
            y(j,i) = py(i) + r * (cos(theta(j))*u(2) + sin(theta(j))*w(2));
            z(j,i) = pz(i) + r * (cos(theta(j))*u(3) + sin(theta(j))*w(3));
        end
    end
end