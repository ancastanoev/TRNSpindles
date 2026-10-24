% TRN_WC_system.m — Right-hand side of the no-delay cortico-thalamo-reticular
% Wilson-Cowan model. Four coupled rate equations for the pyramidal (P),
% interneuron (I), thalamic-relay (T) and reticular (R) populations.
% Variables use thesis notation. Parameters are passed in as globals by the
% solver. ODE function for ode45. State x = [P, I, T, R].
%
% Couplings w_XY (from Y to X), sigmoid params alpha/theta/beta per class,
% drives B/IC/IT/IR, saturation r_P..r_R, time constant tau.

function dx = TRN_WC_system(t, x)

global alpha_E theta_E beta_E alpha_I theta_I beta_I alpha_R theta_R beta_R ...
       wPP wPI wPT wIP wII wIT wTP wTR wRP wRT wRR ...
       tau B IC IT IR r_P r_I r_T r_R;

dx = zeros(1,4);
P = x(1); I = x(2); T = x(3); R = x(4);

dP = 1/tau*(-P + (1 - r_P*P)*Sigmo(alpha_E, theta_E, beta_E, wPP*P - wPI*I + wPT*T + B + IC));
dI = 1/tau*(-I + (1 - r_I*I)*Sigmo(alpha_I, theta_I, beta_I, wIP*P - wII*I + wIT*T + B));
dT = 1/tau*(-T + (1 - r_T*T)*Sigmo(alpha_E, theta_E, beta_E, wTP*P - wTR*R + B + IT));
dR = 1/tau*(-R + (1.5 - r_R*(R-1))*Sigmo(alpha_R, theta_R, beta_R, wRP*P + wRT*T - wRR*R + B + IR));

dx(1) = dP;
dx(2) = dI;
dx(3) = dT;
dx(4) = dR;

dx = dx';
end