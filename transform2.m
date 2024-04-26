function outputSignal = transform2(inputSignal, wlen, L, Nfft)
% transform2 - Transform2 for MCA: Short-Time Fourier Transform (STFT) with normalization
%
% Syntax:
%   outputSignal = transform2(inputSignal, wlen, L, Nfft)
%
% Input:
%   inputSignal - Input signal for the transform
%   wlen - Window length used in the STFT
%   L - Overlap length in the STFT
%   Nfft - FFT length in the STFT
%
% Output:
%   outputSignal - Transformed signal
%
% Description:
%   This function performs the Short-Time Fourier Transform (STFT) on the input
%   signal with the specified window length, overlap length, and FFT length.
%   The result is then flattened to a column vector and normalized by the
%   square root of the window length.


    % Perform STFT on the input signal with specified window, overlap, and FFT length
    outputSignal = stft(inputSignal, 'Window', hamming(wlen), 'OverlapLength', L, 'FFTLength', Nfft);
    
    % Flatten the STFT matrix to a column vector and normalize by the square root of the window length
    outputSignal = outputSignal(:) / sqrt(wlen);

end