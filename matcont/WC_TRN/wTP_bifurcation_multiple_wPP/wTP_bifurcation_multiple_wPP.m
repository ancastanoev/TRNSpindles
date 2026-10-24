% Bifurcation diagram of the WC equilibria under varying
% corticothalamic drive w_TP. Two panels show the T and R components (j = 3, 4)
% of the equilibrium branch with the limit-cycle envelope, for three cortical
% self-couplings w_PP = 1.0 (brown), 0.75 (blue), 0.5 (orange).
% Fixed: w_TR = 0.8, w_RT = 0.4, w_RR = 0.2.

figure('Position',[100 50 1500 550]);

for panel = 1:2
    if panel == 1
        j = 3;
        subplot(1,2,1);
        axisLimits = [0 6 -0.5 6];
        ylab = '$T^*$ component of equilibrium';
    else
        j = 4;
        subplot(1,2,2);
        axisLimits = [0 6 2 24];
        ylab = '$R^*$ component of equilibrium';
    end

    axis(axisLimits);
    ax = gca; ax.FontSize = 12;
    xlabel('Connectivity $w_{TP}$','FontSize',20,'Interpreter','latex');
    ylabel(ylab,'FontSize',20,'Interpreter','latex');
    hold on;
    box on;

    % wPP = 1 (brown)
    load('wTP_cycles_wPP_1_00_wRR_0_2(a).mat')
    L = 41;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, [1.00,0.07,0.65], 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color',[0.64,0.08,0.18],'LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color',[0.64,0.08,0.18], 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    L = 71;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:72
        m(k) = min(x(j:4:end-4,k+40));
        M(k) = max(x(j:4:end-4,k+40));
    end
    y_original = x(646,41:112)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y,'w', 'EdgeColor', 'none','FaceAlpha',0.2);
    plot(y_sorted, curve1_sorted, 'Color',[0.64,0.08,0.18],'LineWidth', 0.5);
    plot(y_sorted, curve2_sorted, 'Color',[0.64,0.08,0.18], 'LineWidth', 0.5);
    p1 = x(end,112); q1 = m(L); p2 = x(end,112); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',0.5,'Color','k');

    load('wTP_cycles_wPP_1_00_wRR_0_2(b).mat')
    L = 39;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, [1.00,0.07,0.65], 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color',[0.64,0.08,0.18],'LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color',[0.64,0.08,0.18], 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    load('wTP_bif_wPP_1_00_wRR_0_2(a).mat');
    plot(x(5,1:34),x(j,1:34),'Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,34:end),x(j,34:end),':','Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,34),x(j,34),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[0.64,0.08,0.18]);
    load('wTP_bif_wPP_1_00_wRR_0_2(b).mat');
    plot(x(5,1:28),x(j,1:28),'Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,28:end),x(j,28:end),':','Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,28),x(j,28),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[0.64,0.08,0.18]);

    % wPP = 0.75 (blue)
    load('wTP_cycles_wPP_0_75_wRR_0_2(a).mat')
    L = 94;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, 'c', 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color','b','LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color','b', 'LineWidth', 1);

    load('wTP_cycles_wPP_0_75_wRR_0_2(b).mat')
    L = 188;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, 'c', 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color','b','LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color','b', 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    L = 298;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:298
        m(k) = min(x(j:4:end-4,k+187));
        M(k) = max(x(j:4:end-4,k+187));
    end
    y_original = x(646,188:485)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y,'w', 'EdgeColor', 'none','FaceAlpha',0.2);
    plot(y_sorted, curve1_sorted,'Color','b','LineWidth', 0.5);
    plot(y_sorted, curve2_sorted, 'Color','b', 'LineWidth', 0.5);
    p1 = x(end,485); q1 = m(L); p2 = x(end,485); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',0.5,'Color','k');

    load('wTP_bif_wPP_0_75_wRR_0_2(a)');
    plot(x(5,1:18),x(j,1:18),'Color','b','Linewidth',1);
    plot(x(5,18:end),x(j,18:end),':','Color','b','Linewidth',1);
    plot(x(5,18),x(j,18),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor','b');
    load('wTP_bif_wPP_0_75_wRR_0_2(b)');
    plot(x(5,1:63),x(j,1:63),'Color','b','Linewidth',1);
    plot(x(5,63:end),x(j,63:end),':','Color','b','Linewidth',1);
    plot(x(5,63),x(j,63),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor','b');

    % wPP = 0.50 (orange)
    load('wTP_cycles_wPP_0_50_wRR_0_2.mat')
    L = 94;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:L
        m(k) = min(x(j:4:end-4,k));
        M(k) = max(x(j:4:end-4,k));
    end
    y_original = x(646,1:L)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y, [0.93,0.69,0.13], 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color',[1.00,0.41,0.16],'LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color',[1.00,0.41,0.16], 'LineWidth', 1);

    load('wTP_bif_wPP_0_50_wRR_0_2(a).mat');
    plot(x(5,:),x(j,:),'Color',[1.00,0.41,0.16],'Linewidth',1);
    load('wTP_bif_wPP_0_50_wRR_0_2(b).mat');
    plot(x(5,1:58),x(j,1:58),'Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,58:end),x(j,58:end),':','Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,58),x(j,58),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[1.00,0.41,0.16]);

end