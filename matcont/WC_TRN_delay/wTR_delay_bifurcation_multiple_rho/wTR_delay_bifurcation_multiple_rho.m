% Delay model: wTR bifurcation with limit cycles, for rho = 0.5, 1, 3, 5.
% Four panels:
%   top row    -> wRT = 0.4 :  T (j=3, left) , R (j=4, right)
%   bottom row -> wRT = 0.2 :  T (j=3, left) , R (j=4, right)
% wPP = 1.1, wRR = 0.2 fixed.
%
% Each panel reproduces the original single-panel scripts
% (wTR_for_wRT_0_4_multiple_rho / wTR_for_wRT_0_2_multiple_rho), changing
% only the component index j and the per-panel axis limits.

figure('Position',[100 50 1400 900]);


panels = { ...
    1, '0_4', 3, [0.15 1 -1 10], '$T^*$ component of equilibrium'; ...
    2, '0_4', 4, [0.15 1 -1 25], '$R^*$ component of equilibrium'; ...
    3, '0_2', 3, [0.15 1 -1 10], '$T^*$ component of equilibrium'; ...
    4, '0_2', 4, [0.15 1 -1 25], '$R^*$ component of equilibrium'  ...
};

for p = 1:4
    sp   = panels{p,1};
    wRT  = panels{p,2};
    j    = panels{p,3};
    alim = panels{p,4};
    ylab = panels{p,5};

    subplot(2,2,sp);
    axis(alim);
    ax = gca; ax.FontSize = 12;
    xlabel('Connectivity $w_{TR}$','FontSize',20,'Interpreter','latex');
    ylabel(ylab,'FontSize',20,'Interpreter','latex');
    hold on; box on;

    if strcmp(wRT,'0_4')
        L_list = [62 68 72 73];
    else
        L_list = [57 58 55 52];
    end
    rho_tags  = {'0_5','1','3','5'};
    fillcols  = {[1.00,0.07,0.65], 'c', 'g', [0.93,0.69,0.13]};
    edgecols  = {[0.64,0.08,0.18], 'b', [0.30,0.51,0.02], [1.00,0.41,0.16]};

    for r = 1:4
        rho = rho_tags{r};
        L   = L_list(r);

        % limit-cycle envelope (H1): per-column min/max over the mesh
        load(['H1_wRT_' wRT '_wPP_1_1_wRR_0_2_rho_' rho '.mat'])
        m = zeros(L,1); M = zeros(L,1);
        for k=1:L
            m(k) = min(x(j:8:end-8,k));
            M(k) = max(x(j:8:end-8,k));
        end
        y_original = x(1290,1:L)';
        [y_sorted, sort_indices] = sort(y_original);
        curve1_sorted = M(sort_indices);
        curve2_sorted = m(sort_indices);
        fill_x = [y_sorted; flipud(y_sorted)];
        fill_y = [curve1_sorted; flipud(curve2_sorted)];
        fill(fill_x, fill_y, fillcols{r}, 'EdgeColor', 'none','FaceAlpha',0.3);
        plot(y_sorted, curve1_sorted, 'Color',edgecols{r},'LineWidth', 1);
        plot(y_sorted, curve2_sorted, 'Color',edgecols{r},'LineWidth', 1);
        p1 = x(end,L); q1 = m(L); p2 = x(end,L); q2 = M(L);
        line([p1,p2],[q1,q2],'Linewidth',1,'Color','k');

        % equilibrium branch (EP1)
        load(['EP1_wRT_' wRT '_wPP_1_1_wRR_0_2_rho_' rho '.mat'])
        if strcmp(wRT,'0_2') && strcmp(rho,'0_5')
            % full branch with stability segments (this panel/rho only)
            plot(x(9,1:163),x(j,1:163),'Color','k','Linewidth',1);
            plot(x(9,163:251),x(j,163:251),':','Color','k','Linewidth',1);
            plot(x(9,251:316),x(j,251:316),'--','Color','k','Linewidth',1);
            plot(x(9,316:end),x(j,316:end),':','Color','k','Linewidth',1);
            plot(x(9,163),x(j,163),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',edgecols{r});
            plot(x(9,251),x(j,251),'s','Markersize',8,'Markeredgecolor','k','Markerfacecolor',edgecols{r});
            plot(x(9,316),x(j,316),'s','Markersize',8,'Markeredgecolor','k','Markerfacecolor',edgecols{r});
        elseif strcmp(wRT,'0_4') && strcmp(rho,'0_5')
            plot(x(9,1:197),x(j,1:197),'Color','k','Linewidth',1);
            plot(x(9,197:end),x(j,197:end),':','Color','k','Linewidth',1);
            plot(x(9,197),x(j,197),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',edgecols{r});
        else
            % other rho: single Hopf marker only
            mk = struct('rho',{'1','3','5'},'i04',{193,204,214},'i02',{158,170,180});
            for q = 1:numel(mk)
                if strcmp(rho,mk(q).rho)
                    if strcmp(wRT,'0_4'), idx = mk(q).i04; else, idx = mk(q).i02; end
                    plot(x(9,idx),x(j,idx),'p','Markersize',10,'Markeredgecolor','k','Markerfacecolor',edgecols{r});
                end
            end
        end
    end
end