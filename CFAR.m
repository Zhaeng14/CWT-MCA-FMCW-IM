function result = CFAR(S)
% CFAR (Constant False Alarm Rate) Algorithm for Radar Interference Detection
% Inputs:
%   S: Input signal (magnitude of the FMCW radar signal)
% Output:
%   result: Threshold values calculated using CFAR algorithm


%% Parameters for Different Data Applications

%% Simulated Data Version (Single Chirp Simulation)

% Nc: Window size (Nc = 50)
% Pa: Probability of false alarm (Pfa) (Pa = 0.005)
% alpha: Coefficient (alpha = 1.2)

%% Hella Company's Real 5th Generation Radar Signal Version

% Nc: Window size (Nc = 40)
% alpha: Coefficient (alpha = 1.5)

%% Hella 6th Generation Radar with DDMA Modulation Signal Version

% Nc: Window size (Nc = 40)
% alpha: Coefficient (alpha = 1.4)


Nc = 40; % Number of training cells (adjust based on system requirements)
alpha = 1.4; % Constant false alarm rate parameter

Rsize = length(S);
threshold = zeros(1, Rsize);

% Iterate through each range bin
for i_R = 1:Rsize
    Rleft = i_R - Nc/2;
    Rright = i_R + Nc/2;
    
    % Ensure the range indices are within bounds
    if Rleft < 1 
        Rleft = 1;
    end
    if Rright > Rsize
        Rright = Rsize;
    end
    
    % Calculate the threshold based on the mean of the training cells
    threshold(i_R) = alpha * mean(S(Rleft:Rright));
end

result = threshold;
end