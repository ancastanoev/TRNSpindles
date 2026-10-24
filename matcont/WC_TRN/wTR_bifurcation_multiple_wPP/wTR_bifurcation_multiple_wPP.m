% Bifurcation diagram of the WC equilibria under varying reticular inhibition
% w_TR, for three cortical self-couplings w_PP = 1.1 (blue), 1.0 (orange),
% 0.9 (brown). Two panels show the T and R components (j = 3, 4) of the
% equilibrium branch with the limit-cycle envelope.
% Fixed: w_TP = 1.2, w_RT = 0.4, w_RR = 0.2.

figure('Position',[100 50 1500 550]);

for panel = 1:2
    if panel == 1
        j = 3;
        subplot(1,2,1);
        axisLimits = [0.5 2 -0.5 3];
        ylab = '$T^*$ component of equilibrium';
    else
        j = 4;
        subplot(1,2,2);
        axisLimits = [0.5 1.2 0 25];
        ylab = '$R^*$ component of equilibrium';
    end

    axis(axisLimits);
    ax = gca; ax.FontSize = 12;
    xlabel('Connectivity $w_{TR}$','FontSize',20,'Interpreter','latex');
    ylabel(ylab,'FontSize',20,'Interpreter','latex');
    hold on;
    box on;

    % wPP = 1.1 (blue)
    load('wTR_cycles_wPP_1_1.mat')
    L = 34;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices);
    curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, 'c', 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color','b','LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color','b', 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    load('wTR_equilibrium_wPP_1_1.mat');
    plot(x(5,1:42),x(j,1:42),'Color','b','Linewidth',1);
    plot(x(5,42:end),x(j,42:end),':','Color','b','Linewidth',1);
    plot(x(5,42),x(j,42),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor','b');

    % wPP = 1.0 (orange) 
    load('wTR_cycles_wPP_1_0.mat')
    L = 44;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices);
    curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, [0.93,0.69,0.13], 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color',[1.00,0.41,0.16],'LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color',[1.00,0.41,0.16], 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    load('wTR_equilibrium_wPP_1_0.mat');
    plot(x(5,1:48),x(j,1:48),'Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,48:end),x(j,48:end),':','Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,48),x(j,48),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[1.00,0.41,0.16]);

    % wPP = 0.9 (brown) 
    load('wTR_cycles_wPP_0_9.mat')
    L = 63;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices);
    curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, [1.00,0.07,0.65], 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color',[0.64,0.08,0.18],'LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color',[0.64,0.08,0.18], 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    load('wTR_equilibrium_wPP_0_9.mat');
    plot(x(5,1:52),x(j,1:52),'Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,52:end),x(j,52:end),':','Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,52),x(j,52),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[0.64,0.08,0.18]);

end