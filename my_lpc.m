function est_x = my_lpc(sb, p)
% LPC Linear Predictive Coding for FMCW Radar Interference Suppression
% Inputs:
%   sb: Input signal containing FMCW radar interference
%   p: Order of the linear predictive coding model
% Output:
%   est_x: Signal with FMCW radar interference suppressed

% Step 1: Apply Constant False Alarm Rate (CFAR) detection
cfar = CFAR(abs(sb)); % Use CFAR to detect interference
outlier = find(abs(sb) > cfar); % Find indices of detected interference

% Step 2: Interference suppression for detected interference
if (~isempty(outlier)) % If interference is detected
    nf = outlier(1);
    if (nf <= 5)
        nf = 5; % Set minimum value for the starting index
    end
    nl = outlier(end);
    if (nl >= length(sb) - 5)
        nl = length(sb) - 5; % Set maximum value for the ending index
    end
    size = nl - nf;
    
    % Step 3: Apply Linear Predictive Coding (LPC) to estimate filter coefficients
    a = lpc([sb(1:nf-1), sb(nl+1:end)], p);
    af = -a(2:end); % Extract coefficients for the forward filter
    ab = conj(flip(af)); % Extract coefficients for the backward filter
    
    % Step 4: Initialize forward and backward signals
    xf = [sb(nf-3:nf-1), zeros(1, size+1)]; % Forward signal initialization
    xb = [zeros(1, size+1), sb(nl+1:nl+3)]; % Backward signal initialization
    
    % Step 5: Apply LPC filtering to forward and backward signals
    for ind = 1:size+1
        xf(3+ind) = af * xf(2+ind:-1:ind).'; % Forward filtering
        xb(end-2-ind) = ab * xb(end+1-ind:-1:end-1-ind).'; % Backward filtering
    end
    
    % Alternatively, you can use the filter function for filtering:
    % xf = filter([0 -af(2:end)], 1, sb(outlier(1)-size:outlier(1)));
    % xb = filter([0 -ab(2:end)], 1, sb(outlier(end):outlier(end)+size));
    
    % Step 6: Combine filtered forward and backward signals to estimate the interference-suppressed signal
    est_x = [sb(1:nf-1), (xf(4:end) + xb(1:end-3))/2, sb(nl+1:end)];
else
    % If no interference is detected, return the original signal
    est_x = sb;
end

end