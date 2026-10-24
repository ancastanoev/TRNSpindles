% TRN_WC_solver.m — Driver for the no-delay Wilson-Cowan TRN model. Sets the
% parameters (thesis notation), integrates TRN_WC_system with ode45, and plots
% the four population activities (P, I, T, R) over a steady-state window.

global alpha_E theta_E beta_E alpha_I theta_I beta_I alpha_R theta_R beta_R ...
       wPP wPI wPT wIP wII wIT wTP wTR wRP wRT wRR ...
       tau B IC IT IR r_P r_I r_T r_R;

alpha_E = 50;   theta_E = 4;    beta_E = 1.2;
alpha_I = 50;   theta_I = 3.8;  beta_I = 2.5;
alpha_R = 50;   theta_R = 3.8;  beta_R = 2.5;

wPP = 1;    wPI = 0.4;  wPT = 0.4;
wIP = 0.7;  wII = 0.1;  wIT = 0.7;
wTP = 1.25; wTR = 0.8;
wRP = 0.7;  wRT = 0.4;  wRR = 0.2;

tau = 1;
B = 0.6;  IC = 0.1;  IT = 2.8;  IR = 1.5;
r_P = 0.05;  r_I = 0.03;  r_T = 0.03;  r_R = 0.05;

init = [0 0 0 0];

options = odeset('MaxStep', .01);
[T, X] = ode45(@TRN_WC_system, [0 300], init, options);

figure('Color', 'w');
hold on;
plot(T(end-7000:end), X(end-7000:end,1), 'b', 'Linewidth', 1);
plot(T(end-7000:end), X(end-7000:end,2), 'r', 'Linewidth', 1);
plot(T(end-7000:end), X(end-7000:end,3), 'g', 'Linewidth', 1);
plot(T(end-7000:end), X(end-7000:end,4), 'k', 'Linewidth', 1);
axis([285 300 0 15]);
xlabel('Time (ms)');
ylabel('Population activity');
title('Wilson-Cowan CTRC dynamics (no delay)');
legend('P (pyramidal)', 'I (interneurons)', 'T (thalamic relay)', 'R (reticular)');
grid off;