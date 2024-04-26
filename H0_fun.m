function H0 = H0_fun(w, alpha, beta)
% H0_fun - Compute H0 for given parameters alpha and beta
%
% Syntax:
%   H0 = H0_fun(w, alpha, beta)
%
% Input:
%   w - Frequency values (radians)
%   alpha - Parameter (0.8 as default)
%   beta - Parameter (0.5 as default)
%
% Output:
%   H0 - Computed H0 values
%
% Example:
%   alpha = 0.8;
%   beta = 0.5;
%   w = linspace(0, pi, 100);
%   H0 = H0_fun(w, alpha, beta);
%   H1 = H1_fun(w, alpha, beta);
%   subplot(2, 1, 1), plot(w/pi, H0, w/pi, H1)
%   subplot(2, 1, 2), plot(H0.^2 + H1.^2 - 1)
%
% See also: H1_fun, theta_fun

H0 = zeros(size(w));

w = mod(w + pi, 2*pi) - pi;

H0(abs(w) <= alpha*pi) = 1;

k = (abs(w) >= (1-beta)*pi) & (abs(w) <= alpha*pi);

H0(k) = theta_fun((w(k) + (beta-1)*pi) / (alpha + beta - 1));
