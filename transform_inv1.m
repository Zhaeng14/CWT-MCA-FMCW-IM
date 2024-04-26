function outputSignal = transform_inv1(inputSignal)
% transform_inv1 - Inverse Transform1 for MCA: Inverse Discrete Fourier Transform (IDFT) with scaling
%
% Syntax:
%   outputSignal = transform_inv1(inputSignal)
%
% Input:
%   inputSignal - Input signal for inverse transform
%
% Output:
%   outputSignal - Inverse transformed signal
%
% Description:
%   This function computes the Inverse Discrete Fourier Transform (IDFT)
%   on the input signal and scales the result by the square root of the input signal length.


    outputSignal = ifft(inputSignal); % Compute the Inverse Discrete Fourier Transform (IDFT)
    outputSignal = outputSignal * sqrt(length(inputSignal)); % Scale the result

end
