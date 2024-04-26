function [result, x] = mca_method(x_mixed)
% mca_method - MCA Method for Interference Suppression in FMCW Radar Signals
%
% Syntax:
%   [result, x] = mca_method(x_mixed)
%
% Input:
%   x_mixed - Input mixed FMCW radar signal
%
% Output:
%   result - Interference-suppressed signal
%   x - Component of the signal in the DFT domain
%
% Description:
%   This function implements the MCA (Matching Pursuit and Convex Analysis) method
%   for interference suppression in FMCW radar signals. It decomposes the mixed
%   signal into two frequency domain components, applying soft thresholding to
%   iteratively suppress interference.
%
% Example:
%   [result, x] = mca_method(x_mixed);
%
% See also: transform1, transform2, transform_inv1, transform_inv2, soft

%% Parameters for Different Data Applications
% Simulated Data Version

% N_iter: Number of iterations (N_iter = 10)
% wlen: STFT window length (wlen = floor(Nfft/64))
% overLap: STFT window overlap (overLap = floor(wlen / 2))
% lambda_p: Final threshold coefficient (lambda_p = 10000)

% Hella Company's Real 5th Generation Radar Signal Version

% N_iter: Number of iterations (N_iter = 5)
% wlen: STFT window length (wlen = floor(Nfft/64))
% overLap: STFT window overlap (overLap = floor(wlen / 2))
% lambda_p: Final threshold coefficient (lambda_p = 10000)

% Hella 6th Generation Radar with DDMA Modulation Signal Version

% N_iter: Number of iterations (N_iter = 4)
% wlen: STFT window length (wlen = floor(Nfft/32))
% overLap: STFT window overlap (overLap = floor(wlen - 1))
% lambda_p: Final threshold coefficient (lambda_p = 10000)




% Initialization
L = length(x_mixed);
Nfft = L;
wlen = floor(Nfft/64); % Window length for STFT (adjust based on system requirements)
overLap = floor(wlen/2); % Overlap for STFT (adjust based on system requirements)
X_prime = x_mixed;
X = cell(1, 2);
alpha = cell(1, 2);
X{1} = zeros(L, 1);
X{2} = zeros(L, 1);
R = x_mixed;
alpha{1} = transform1(x_mixed);
alpha{2} = transform2(x_mixed, wlen, overLap, Nfft);
a_max1 = max(abs(alpha{1}));
a_max2 = max(abs(alpha{2}));
lambda1 = a_max1;
lambda2 = a_max2;

lambda_init1 = lambda1;
lambda_init2 = lambda2;

lambda_p = 10000;
lambda_min1 = lambda_init1 / lambda_p;
lambda_min2 = lambda_init2 / lambda_p;

N_iter = 5;
R_ = zeros(1, N_iter);

% Iterative process
for i = 1 : N_iter
    Res{1} = R - X{2};
    Res{2} = R - X{1};
    
    % Apply soft thresholding to the transformed components
    alpha{1} = soft(transform1(Res{1}), lambda1);
    alpha{2} = soft(transform2(Res{2}, wlen, overLap, Nfft), lambda2); 
    
    % Inverse transform to obtain updated components
    X{1} = transform_inv1(alpha{1});
    X{2} = transform_inv2(alpha{2}, wlen, overLap, Nfft);
    
    X_recon = X{1} + X{2};
    R = X_prime - X_recon;
    R_(i) = norm(R);
    
    % Update threshold parameters for the next iteration
    lambda1 = lambda1 / 2;
    lambda2 = lambda2 / 2;
end

% Scale the obtained components
X{1} = 2 * X{1};
X{2} = 2 * X{2};


% Output the interference-suppressed signal
x = X{1};
result = X{1};
end

% Soft thresholding function
function result = soft(y, T)
    result = wthresh(y, 's', T);
end

% Hard thresholding function (commented out, as it is not used)
% function result = hard(y, T)
%     result = wthresh(y, 'h', T);
% end
