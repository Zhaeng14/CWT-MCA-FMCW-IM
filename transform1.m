function outputSignal = transform1(inputSignal)
% transform1 - Transform1 for MCA: Discrete Fourier Transform (DFT) normalization
%
% Syntax:
%   outputSignal = transform1(inputSignal)
%
% Input:
%   inputSignal - Input signal for the transform
%
% Output:
%   outputSignal - Transformed signal
%
% Description:
%   This function computes the Discrete Fourier Transform (DFT) on the input
%   signal and normalizes the result by dividing by the square root of the input signal length.


    outputSignal = fft(inputSignal); % Compute the Discrete Fourier Transform (DFT)
    outputSignal = outputSignal / sqrt(length(inputSignal)); % Normalize the DFT by dividing by the square root of the signal length

end
