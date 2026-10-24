% Hopf stability boundary in the (w_TR, w_RT) connectivity plane, traced for
% several values of the mean delay rho. Each curve is a two-parameter Hopf
% continuation computed in MatCont; together they show how increasing delay
% deforms the boundary between stable and oscillatory regimes.
%
% DATA (same folder):
%   hopf_boundary_rho_0_branch1.mat, _branch2.mat
%       rho = 0 reference (no-delay Hopf boundary), two branches.
%       Axes stored in rows 5 (w_TR) and 6 (w_RT).
%   hopf_boundary_rho_0_02 / _0_06 / _0_1 / _0_2 / _0_5 / _5_0.mat
%       rho > 0 boundaries. Axes stored in rows 9 (w_TR) and 10 (w_RT).
%
% OUTPUT:
%   One figure: seven boundary curves over [0,1.4] x [0,1.4], with legend.

figure('Position',[100 50 760 620]);
hold on;
axis([0 1.4 0 1.4]);
axis square;
box on;
ax = gca; ax.FontSize = 12;
xlabel('$w_{TR}$','FontSize',20,'Interpreter','latex');
ylabel('$w_{RT}$','FontSize',20,'Interpreter','latex');

H1 = load('hopf_boundary_rho_0_branch1.mat');
H2 = load('hopf_boundary_rho_0_branch2.mat');
h0 = plot(H1.x(5,:), H1.x(6,:), 'Color','k', 'LineWidth',1.2);
     plot(H2.x(5,:), H2.x(6,:), 'Color','k', 'LineWidth',1.2);

rho_files = { ...
'hopf_boundary_rho_0_02.mat', [0.00 0.20 0.80]; ...
'hopf_boundary_rho_0_06.mat', [0.20 0.40 1.00]; ...
'hopf_boundary_rho_0_1.mat',  [0.45 0.20 0.80]; ...
'hopf_boundary_rho_0_2.mat',  [0.70 0.10 0.45]; ...
'hopf_boundary_rho_0_5.mat',  [0.90 0.30 0.20]; ...
'hopf_boundary_rho_5_0.mat',  [1.00 0.55 0.20]  ...
};
rho_labels = {'$\rho=0.02$','$\rho=0.06$','$\rho=0.1$', ...
'$\rho=0.2$','$\rho=0.5$','$\rho=5$'};

handles = gobjects(1, size(rho_files,1));
for r = 1:size(rho_files,1)
    d = load(rho_files{r,1});
    handles(r) = plot(d.x(9,:), d.x(10,:), 'Color', rho_files{r,2}, 'LineWidth',1.2);
end
legend([h0, handles], [{'$\rho=0$'}, rho_labels], ...
'Interpreter','latex', 'FontSize',12, 'Location','northeast');