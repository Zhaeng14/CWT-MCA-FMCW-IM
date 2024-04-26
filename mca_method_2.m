function [result, x] = mca_method_2(x_mixed)
% mca_method_2 - MOM-MCA Algorithm Implementation
%
% Syntax:
%   [result, x] = mca_method_2(x_mixed)
%
% Input:
%   x_mixed - Mixed input signal
%
% Output:
%   result - Final result with equal weighting
%   x - Frequency domain components
%
% Description:
%   This function implements the MOM-MCA algorithm for blind source separation.
%   It decomposes the mixed input signal into two frequency domain components,
%   iteratively updating the components using thresholding techniques.
%
% Example:
%   [result, x] = mca_method_2(x_mixed);
%
% See also: transform1, transform2, transform_inv1, transform_inv2, soft, hard

%% Parameters for Different Data Applications

% Simulated Data Version

% N_iter: Number of iterations (N_iter = 30)
% wlen: STFT window length (wlen = floor(Nfft/64))
% overLap: STFT window overlap (overLap = floor(wlen / 2))

% Hella Company's Real 5th Generation Radar Signal Version

% N_iter: Number of iterations (N_iter = 30)
% wlen: STFT window length (wlen = floor(Nfft/64))
% overLap: STFT window overlap (overLap = floor(wlen / 2))

% Hella 6th Generation Radar with DDMA Modulation Signal Version

% N_iter: Number of iterations (N_iter = 15)
% wlen: STFT window length (wlen = floor(Nfft/64))
% overLap: STFT window overlap (overLap = floor(wlen - 3))



%% Initialization
L = length(x_mixed); % Length of the input signal
Nfft = L; % FFT length
wlen = floor(Nfft/64); % Window length for STFT
overLap = floor(wlen/2); % Overlap length for STFT
X_prime = x_mixed; % Store the original mixed signal
X = cell(1, 2); % Cell array to store frequency domain components
alpha1 = 0; % Initialization for frequency domain component 1
alpha2 = 0; % Initialization for frequency domain component 2

X{1} = zeros(L, 1); % Initialization for frequency domain component 1
X{2} = zeros(L, 1); % Initialization for frequency domain component 2
R = x_mixed; % Residual signal
N_iter = 30; % Number of iterations
R_ = zeros(1, N_iter); % Array to store residuals across iterations
i = 0; % Iteration counter

%% Iteration
while((i == 0 || i == 1 || R_(i) < R_(i-1)) && i < N_iter)
    i = i + 1;
    Res{1} = R - X{2}; % Residual for frequency domain component 1
    Res{2} = R - X{1}; % Residual for frequency domain component 2
    
    m1 = norm(transform1(R), inf); % Maximum value of transformed signal 1
    m2 = norm(transform2(R, wlen, overLap, Nfft), inf); % Maximum value of transformed signal 2
    lambda = (m1 + m2) / 2; % Threshold parameter
    
    if(m1 > m2)
        alpha1 = alpha1 + hard(transform1(Res{1}), lambda); % Update frequency domain component 1
        X{1} = transform_inv1(alpha1); % Inverse transform for component 1
    else
        alpha2 = alpha2 + hard(transform2(Res{2}, wlen, overLap, Nfft), lambda); % Update frequency domain component 2
        X{2} = transform_inv2(alpha2, wlen, overLap, Nfft); % Inverse transform for component 2
    end
    
    X_recon = X{1} + X{2}; % Reconstructed signal
    R = X_prime - X_recon; % Update residual
    R_(i) = norm(R); % Store current residual norm
end

result = X{1} + X{2} / 3; % Final result with equal weighting
x = X{1} + X{2} / 3; % Frequency domain components

end

function result = soft(y, T)
    result = wthresh(y, 's', T); % Soft thresholding function
end

function result = hard(y, T)
    result = wthresh(y, 'h', T); % Hard thresholding function
end
