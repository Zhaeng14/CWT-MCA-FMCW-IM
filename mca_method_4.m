function [result, y1] = mca_method_4(x)
% MCA method for CWT (Continuous Wavelet Transform)
% Separates target signal from interference based on the sparsity differences
% in different CWT domains.
%
% INPUT
%   x - Input mixed FMCW radar signal
%
% OUTPUT
%   result - Processed signal after interference suppression
%   y1 - First component obtained from the CWT-MCA algorithm

%% Simulated Data
    % Nit = 20;    
    % lam1 = 0.1;         
    % lam2 = 0.1;
    % mu = 0.04;

%% Hella 5th Generation Radar Signal
    % Parameters
    % Nit = 4;     
    % lam1 = 0.1;         
    % lam2 = 0.1;
    % mu = 0.04;
    % N = length(x);
    % dt = 1.2e-5/512;
    % pad = 0;
    % dj1 = 0.2;
    % dj2 = 1.3;

%% Hella 6th Generation Radar Signal with DDMA Modulation
    % Parameters
    % Nit = 4;     
    % lam1 = 0.1;         
    % lam2 = 0.1;
    % mu = 0.04;
    % N = length(x);
    % dt = 1.2e-5/512;
    % pad = 0;
    % dj1 = 0.2;
    % dj2 = 1.3;







%% Parameters for the CWT-MCA algorithm

% Parameters for the first component
Q1 = 20;        % Q-factor
r1 = 3;         % Oversampling rate
L1 = 56;        % Number of levels

% Parameters for the second component
Q2 = 1;         % Q-factor
r2 = 3;         % Oversampling rate
L2 = 10;        % Number of levels

Nit = 4;        % Number of iterations
lam1 = 0.1;     % Regularization parameter for the first component
lam2 = 0.1;     % Regularization parameter for the second component
mu = 0.04;      % SALSA parameter
N = length(x);  % Length of the input signal

% Time step and other CWT parameters
dt = 1.2e-5 / 512;  % Time step
pad = 1;            % Padding
dj1 = 0.2;          % Smaller number gives better resolution for the first component
dj2 = 1.3;          % Dj parameter for the second component

%% Apply CWT-MCA algorithm to real and imaginary parts separately

% Apply CWT-MCA to the real part
[y1_real,y2_real,w1s,w2s,costfn] = cwt_mca(real(x),dt,pad,dj1,dj2,lam1,lam2,mu,Nit,'donotplots');

% Apply CWT-MCA to the imaginary part

[y1_imag,y2_imag,w1s,w2s,costfn] = cwt_mca(imag(x),dt,pad,dj1,dj2,lam1,lam2,mu,Nit,'donotplots');
% [y1_real,y2_real,w1s,w2s,costfn] = fcwt_mca(real(x),dt,pad,dj1,dj2,lam1,lam2,mu,Nit,'donotplots');
% [y1_imag,y2_imag,w1s,w2s,costfn] = fcwt_mca(imag(x),dt,pad,dj1,dj2,lam1,lam2,mu,Nit,'donotplots');

% Combine real and imaginary parts to obtain the final result
y1 = y1_real + 1j * y1_imag;
y2 = y2_real + 1j * y2_imag;
% The result is the first component obtained from the CWT-MCA algorithm
result = y1;

end
