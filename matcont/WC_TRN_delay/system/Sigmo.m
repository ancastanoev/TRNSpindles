% Sigmo.m — Baseline-subtracted sigmoid activation for the Wilson-Cowan TRN
% model. Maps net synaptic input x to a firing rate via a logistic curve of
% amplitude A (alpha), threshold theta, and slope b (beta); the constant term
% enforces f = 0 at zero input.

function f = Sigmo(A, theta, b, x)

f = A./(1 + exp(-b*(x - theta))) - A./(1 + exp(b*theta));
end