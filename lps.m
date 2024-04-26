function Y = lps(X, N0)
% lps - Low-pass scaling
%
% Syntax:
%   Y = lps(X, N0)
%
% Input:
%   X - Input signal
%   N0 - Desired length of the output signal
%
% Output:
%   Y - Low-pass scaled signal
%
% Notes:
%   - The output Y will be of length N0.
%   - The length of X should be even.
%
% Example:
%   Y = lps(X, N0);


N = length(X);

Y = nan(1, N0);

% Add 1 to indices because Matlab indexing starts at 1 (not 0)

switch 1
    case N0 <= N
        k = 0:N0/2-1;
        Y(k + 1) = X(k + 1);
        Y(N0/2 + 1) = X(N/2 + 1);
        k = 1:N0/2-1;
        Y(N0 - k + 1) = X(N - k + 1);
        
    case N0 >= N
        k = 0:N/2-1;
        Y(k + 1) = X(k + 1);
        k = N/2:N0/2-1;
        Y(k + 1) = 0;
        Y(N0/2 + 1) = X(N/2 + 1);
        Y(N0 - k + 1) = 0;
        k = 1:N/2-1;
        Y(N0 - k + 1) = X(N - k + 1); 
end
