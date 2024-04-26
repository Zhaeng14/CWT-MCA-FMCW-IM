function y = itqwt_radix2(w, Q, r, L)
% itqwt_radix2 - Inverse radix-2 TQWT (Translation-Invariant Wavelet Transform)
%
% Syntax:
%   y = itqwt_radix2(w, Q, r, L)
%
% Input:
%   w - Cell array containing wavelet coefficients
%   Q - Quality factor
%   r - Oversampling factor
%   L - Desired length of the output signal
%
% Output:
%   y - Inverse TQWT transformed signal
%
% Description:
%   This function performs the inverse radix-2 TQWT on the input wavelet
%   coefficients, reconstructing the original signal.
%
% Example:
%   y = itqwt_radix2(w, Q, r, L);
%
% See also: tqwt_radix2

beta = 2 / (Q + 1);
alpha = 1 - beta / r;
J = length(w) - 1;

N = next(L);           % *

Y = uDFT(w{J + 1});

M = 2 * round(alpha^J * N / 2);     % *
Y = lps(Y, M);                     % *

for j = J:-1:1
    W = uDFT(w{j});
    N1 = 2 * round(beta * alpha^(j - 1) * N / 2);     % *
    W = lps(W, N1);                                  % *
    M = 2 * round(alpha^(j - 1) * N / 2);
    Y = sfb(Y, W, M);
end

y = uDFTinv(Y);
y = y(1:L);                     % *

% * : denotes modification for radix-2 case
