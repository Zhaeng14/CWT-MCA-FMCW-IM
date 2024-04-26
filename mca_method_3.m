function [result, x] = mca_method_3(x)
% TQWT-MCA Algorithm Implementation
%
% INPUTS:
%   x: Input signal to be processed
%
% OUTPUTS:
%   result: Output of the TQWT-MCA algorithm, representing the suppressed interference
%   x: Updated input signal after interference suppression

%% Parameters for Different Data Applications

% Simulated Data Version

% Q1: High Q-factor (Q1 = 20)
% r1: Redundancy (typically r >= 3) (r1 = 3)
% L1: Number of sub-bands (L1 = 56)

% Q2: Low Q-factor (Q2 = 1)
% r2: Redundancy (typically r >= 3) (r2 = 3)
% L2: Number of sub-bands (L2 = 10)

% Nit: Number of iterations (Nit = 20)
% lam1: Lambda 1 (lam1 = 0.1)
% lam2: Lambda 2 (lam2 = 0.1)
% mu: Mu (mu = 0.04)

% Hella Company's Real 5th Generation Radar Signal Version

% Q1: High Q-factor (Q1 = 20)
% r1: Redundancy (typically r >= 3) (r1 = 3)
% L1: Number of sub-bands (L1 = 56)

% Q2: Low Q-factor (Q2 = 1)
% r2: Redundancy (typically r >= 3) (r2 = 3)
% L2: Number of sub-bands (L2 = 10)

% Nit: Number of iterations (Nit = 10)
% lam1: Lambda 1 (lam1 = 0.1)
% lam2: Lambda 2 (lam2 = 0.1)
% mu: Mu (mu = 0.04)

% Hella 6th Generation Radar with DDMA Modulation Signal Version

% Q1: High Q-factor (Q1 = 20)
% r1: Redundancy (typically r >= 3) (r1 = 3)
% L1: Number of sub-bands (L1 = 56)

% Q2: Low Q-factor (Q2 = 1)
% r2: Redundancy (typically r >= 3) (r2 = 3)
% L2: Number of sub-bands (L2 = 10)

% Nit: Number of iterations (Nit = 10)
% lam1: Lambda 1 (lam1 = 0.1)
% lam2: Lambda 2 (lam2 = 0.1)
% mu: Mu (mu = 0.04)



%% Parameters
Q1 = 20; % Number of subbands for TQWT 1
r1 = 3; % Filter length ratio for TQWT 1
L1 = 56; % Number of filters for TQWT 1
% L1 = 34;

Q2 = 1; % Number of subbands for TQWT 2
r2 = 3; % Filter length ratio for TQWT 2
L2 = 10; % Number of filters for TQWT 2 (Original value)
% L2 = 8;

Nit = 10; % Number of iterations
lam1 = 0.1; % Regularization parameter for TQWT 1
lam2 = 0.1; % Regularization parameter for TQWT 2
mu = 0.04; % Step size parameter
N = length(x); % Signal length
% N = 512;

% Compute wavelets for TQWT 1 and TQWT 2
[wlets1, now1] = ComputeWavelets(N, Q1, r1, L1, 'radix2');
[wlets2, now2] = ComputeWavelets(N, Q2, r2, L2, 'radix2');

% Apply TQWT-MCA algorithm
[y1, y2, w1s, w2s, costfn] = tqwt_mca(x, Q1, r1, L1, Q2, r2, L2, lam1 * now1, lam2 * now2, mu, Nit, 'donotplots');
result = y2; % Final result
x = y2; 
end