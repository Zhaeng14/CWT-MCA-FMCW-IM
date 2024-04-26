function [wlets, now] = ComputeWavelets(N, Q, r, J, radix_flag)
% ComputeWavelets - Compute wavelets for the TQWT (Tunable Q-Factor Wavelet Transform)
%
% Syntax:
%   wlets = ComputeWavelets(N, Q, r, J)
%   Use ComputeWavelets(N, Q, r, J, 'radix2') for radix-2 TQWT
%
% [wlets, now] = ComputeWavelets(...) also returns a vector of
% the norms of the wavelets for each subband.
%
% Input:
%   N - Length of the input signal
%   Q - Quality factor for the TQWT
%   r - Oversampling factor for the TQWT
%   J - Number of decomposition levels
%   radix_flag - Optional flag for using radix-2 TQWT ('radix2' for radix-2)
%
% Output:
%   wlets - Cell array containing wavelets for each subband
%   now - Vector of the norms of the wavelets for each subband
%
% Example:
%   Q = 1; r = 3; J = 10; N = 2^9;   % or
%   Q = 4; r = 3; J = 18; N = 2^9;
%   wlets = ComputeWavelets(N, Q, r, J, 'radix2');
%   figure(1), clf,
%   for j = 1:J+1, line(1:N, wlets{j}/max(abs(wlets{j}))/2 + j); end; xlim([0 N])
%
% See also: tqwt, itqwt, tqwt_radix2, itqwt_radix2

if nargin == 5
    if strcmp(radix_flag, 'radix2')
        xform = @tqwt_radix2;
        inv_xform = @itqwt_radix2;
        C = N / next(N);                      % reduce according to the amount of zero-padding
    else
        disp('Invalid string')
        wlets = [];
        return
    end
else
    xform = @tqwt;
    inv_xform = @itqwt;
    C = 1;
end

z = zeros(1, N);                     % Zero signal

wz = xform(z, Q, r, J);        % All-zero wavelet coefficients

if isempty(wz)
    wlets = [];
    return
end

wlets = cell(1, J+1);

for j = 1:J+1
    w = wz;                         % Set w to all-zero coefficients
    m = round(C * length(w{j}) / 2) + 1;      % m: index of coefficient around the midpoint of the signal
    w{j}(m) = 1;                    % Set single wavelet coeff to 1
    y = inv_xform(w, Q, r, N);     % Inverse rational WT
    wlets{j} = y;
end

now = zeros(1, J+1);
for j = 1:J+1
    now(j) = sqrt(sum(abs(wlets{j}).^2));       % L2 norm of the wavelet
end
