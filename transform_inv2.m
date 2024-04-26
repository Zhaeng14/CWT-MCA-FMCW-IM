function outputSignal = transform_inv2(inputSignal, wlen, L, Nfft)
% transform_inv2 - Inverse Transform2 for MCA: Inverse Short-Time Fourier Transform (ISTFT) with scaling
%
% Syntax:
%   outputSignal = transform_inv2(inputSignal, wlen, L, Nfft)
%
% Input:
%   inputSignal - Input signal for inverse transform
%   wlen - Window length used in the STFT
%   L - Overlap length in the STFT
%   Nfft - FFT length in the STFT
%
% Output:
%   outputSignal - Inverse transformed signal
%
% Description:
%   This function performs the Inverse Short-Time Fourier Transform (ISTFT)
%   on the input signal with the specified window length, overlap length,
%   and FFT length. The result is scaled by the square root of the window length.


    k = floor((Nfft - L) / (wlen - L));
    
    % Reshape the input signal
    inputSignal = reshape(inputSignal, Nfft, k);
    
    % Perform ISTFT on the reshaped input signal
    outputSignal = istft(inputSignal, 'Window', hamming(wlen), 'OverlapLength', L, 'FFTLength', Nfft);
    
    % Scale the result by the square root of the window length
    outputSignal = outputSignal * sqrt(wlen);

end
