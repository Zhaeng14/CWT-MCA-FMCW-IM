function m = next(k)
% next - Compute the next power of two
%
% Syntax:
%   m = next(k)
%
% Input:
%   k - Input value
%
% Output:
%   m - Next power of two
%
% Description:
%   This function computes the next power of two greater than or equal to
%   the input value 'k'.

m = 2.^ceil(log2(k));
