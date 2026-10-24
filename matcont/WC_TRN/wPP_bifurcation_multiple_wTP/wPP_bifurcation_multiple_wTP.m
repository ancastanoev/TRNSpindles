% Bifurcation diagram of the WC equilibria under varying cortical
% self-coupling w_PP. Each of the four panels plots one component (P, I, T, R;
% j = 1..4) of the equilibrium branch together with the limit-cycle envelope,
% for three corticothalamic drives w_TP = 1.0 (orange), 1.25 baseline (blue),
% 1.5 (brown). Fixed: w_TR = 0.8, w_RT = 0.4, w_RR = 0.2.

figure('Position',[100 50 1400 900]);

for j = 1:4
    subplot(2,2,j);

    switch j
        case 1
            axisLimits = [0.7 1.4 2 15];
            ylab = '$P^*$ component of equilibrium';
        case 2
            axisLimits = [0.7 1.4 2 22];
            ylab = '$I^*$ component of equilibrium';
        case 3
            axisLimits = [0.7 1.4 -1 7];
            ylab = '$T^*$ component of equilibrium';
        case 4
            axisLimits = [0.7 1.4 4 24];
            ylab = '$R^*$ component of equilibrium';
    end

    axis(axisLimits);
    ax = gca; ax.FontSize = 12;
    xlabel('Connectivity $w_{PP}$','FontSize',20,'Interpreter','latex');
    ylabel(ylab,'FontSize',20,'Interpreter','latex');
    hold on;
    box on;

    % wTP = 1 (orange)
    load('wPP_cycles_wTP_1_00_wRR_0_2.mat')
    L = 23;
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
    fill(fill_x, fill_y,[0.93,0.69,0.13], 'EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color',[1.00,0.41,0.16],'LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color',[1.00,0.41,0.16], 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    L = 31;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:31
        m(k) = min(x(j:4:end-4,k+22));
        M(k) = max(x(j:4:end-4,k+22));
    end
    y_original = x(646,23:53)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y,'w', 'EdgeColor', 'none','FaceAlpha',0.2);
    plot(y_sorted, curve1_sorted,'Color',[1.00,0.41,0.16],'LineWidth', 0.5);
    plot(y_sorted, curve2_sorted, 'Color',[1.00,0.41,0.16], 'LineWidth', 0.5);
    p1 = x(end,53); q1 = m(L); p2 = x(end,53); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',0.5,'Color','k');

    load('wPP_bif_wTP_1_00_wRR_0_2.mat');
    plot(x(5,1:73),x(j,1:73),'Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,73:246),x(j,73:246),':','Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,246:333),x(j,246:333),'--','Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,333:end),x(j,333:end),'Color',[1.00,0.41,0.16],'Linewidth',1);
    plot(x(5,73),x(j,73),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[1.00,0.41,0.16]);
    plot(x(5,246),x(j,246),'o','Markersize',6,'Markeredgecolor','k','Markerfacecolor',[1.00,0.41,0.16]);
    plot(x(5,333),x(j,333),'o','Markersize',6,'Markeredgecolor','k','Markerfacecolor',[1.00,0.41,0.16]);

    %wTP = 1.25 (blue)
    load('wPP_cycles_wTP_1_25_wRR_0_2.mat')
    L = 23;
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
    fill(fill_x, fill_y,'c','EdgeColor', 'none','FaceAlpha',0.6);
    plot(y_sorted, curve1_sorted, 'Color','b','LineWidth', 1);
    plot(y_sorted, curve2_sorted, 'Color','b', 'LineWidth', 1);
    p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

    L = 31;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:31
        m(k) = min(x(j:4:end-4,k+22));
        M(k) = max(x(j:4:end-4,k+22));
    end
    y_original = x(646,23:53)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y,'w', 'EdgeColor', 'none','FaceAlpha',0.2);
    plot(y_sorted, curve1_sorted, 'Color',[0.64,0.08,0.18],'LineWidth', 0.5);
    plot(y_sorted, curve2_sorted, 'Color',[0.64,0.08,0.18], 'LineWidth', 0.5);
    p1 = x(end,53); q1 = m(L); p2 = x(end,53); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',0.5,'Color','k');

    load('wPP_bif_wTP_1_25_wRR_0_2.mat');
    plot(x(5,1:76),x(j,1:76),'Color','b','Linewidth',1);
    plot(x(5,76:238),x(j,76:238),':','Color','b','Linewidth',1);
    plot(x(5,238:374),x(j,238:374),'--','Color','b','Linewidth',1);
    plot(x(5,374:end),x(j,374:end),'Color','b','Linewidth',1);
    plot(x(5,76),x(j,76),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor','b');
    plot(x(5,238),x(j,238),'o','Markersize',6,'Markeredgecolor','k','Markerfacecolor','b');
    plot(x(5,374),x(j,374),'o','Markersize',6,'Markeredgecolor','k','Markerfacecolor','b');

    % wTP = 1.50 (brown)
    load('wPP_cycles1_wTP_1_50_wRR_0_2.mat')
    L = 18;
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

    L = 31;
    m = zeros(L,1); M = zeros(L,1);
    for k=1:31
        m(k) = min(x(j:4:end-4,k+17));
        M(k) = max(x(j:4:end-4,k+17));
    end
    y_original = x(646,18:48)';
    [y_sorted, sort_indices] = sort(y_original);
    curve1_sorted = M(sort_indices); curve2_sorted = m(sort_indices);
    fill_x = [y_sorted; flipud(y_sorted)];
    fill_y = [curve1_sorted; flipud(curve2_sorted)];
    fill(fill_x, fill_y,'w', 'EdgeColor', 'none','FaceAlpha',0.2);
    plot(y_sorted, curve1_sorted, 'Color',[0.64,0.08,0.18],'LineWidth', 0.5);
    plot(y_sorted, curve2_sorted, 'Color',[0.64,0.08,0.18], 'LineWidth', 0.5);
    p1 = x(end,48); q1 = m(L); p2 = x(end,48); q2 = M(L);
    line([p1,p2],[q1,q2],'Linewidth',0.5,'Color','k');

    load('wPP_bif_wTP_1_50_wRR_0_2.mat');
    plot(x(5,1:73),x(j,1:73),'Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,73:217),x(j,73:217),':','Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,217:483),x(j,217:483),'--','Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,483:end),x(j,483:end),'Color',[0.64,0.08,0.18],'Linewidth',1);
    plot(x(5,73),x(j,73),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',[0.64,0.08,0.18]);
    plot(x(5,217),x(j,217),'o','Markersize',6,'Markeredgecolor','k','Markerfacecolor',[0.64,0.08,0.18]);
    plot(x(5,483),x(j,483),'o','Markersize',6,'Markeredgecolor','k','Markerfacecolor',[0.64,0.08,0.18]);

end