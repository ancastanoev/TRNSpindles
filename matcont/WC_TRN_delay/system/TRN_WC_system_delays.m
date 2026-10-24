% TRN_WC_system_delays.m — Right-hand side of the distributed-delay Wilson-Cowan
% TRN model. The four populations (P, I, T, R) couple through delayed copies
% (P_d, I_d, T_d, R_d) from a first-order linear-chain filter with mean delay rho.
% Variables use thesis notation. Globals set by the solver. ODE function for
% ode45. State x = [P, I, T, R, P_d, I_d, T_d, R_d].

function dx = TRN_WC_system_delays(t, x)

global alpha_E theta_E beta_E alpha_I theta_I beta_I alpha_R theta_R beta_R ...
       wPP wPI wPT wIP wII wIT wTP wTR wRP wRT wRR ...
       tau B IC IT IR r_P r_I r_T r_R rho;

dx = zeros(1,8);
P = x(1); I = x(2); T = x(3); R = x(4);
P_d = x(5); I_d = x(6); T_d = x(7); R_d = x(8);

dP = 1/tau*(-P + (1 - r_P*P)*Sigmo(alpha_E, theta_E, beta_E, wPP*P_d - wPI*I_d + wPT*T_d + B + IC));
dI = 1/tau*(-I + (1 - r_I*I)*Sigmo(alpha_I, theta_I, beta_I, wIP*P_d - wII*I_d + wIT*T_d + B));
dT = 1/tau*(-T + (1 - r_T*T)*Sigmo(alpha_E, theta_E, beta_E, wTP*P_d - wTR*R_d + B + IT));
dR = 1/tau*(-R + (1 - r_R*R)*Sigmo(alpha_R, theta_R, beta_R, wRP*P_d + wRT*T_d - wRR*R_d + B + IR));
dP_d = 1/rho*(P - P_d);
dI_d = 1/rho*(I - I_d);
dT_d = 1/rho*(T - T_d);
dR_d = 1/rho*(R - R_d);

dx(1) = dP;
dx(2) = dI;
dx(3) = dT;
dx(4) = dR;
dx(5) = dP_d;
dx(6) = dI_d;
dx(7) = dT_d;
dx(8) = dR_d;

dx = dx';
end