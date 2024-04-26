function check_params(Q, r, J)
% check_params - Check validity of input parameters for TQWT
%
% Syntax:
%   check_params(Q, r, J)
%
% Input:
%   Q - Quality factor for the TQWT
%   r - Oversampling factor for the TQWT
%   J - Number of decomposition levels for the TQWT
%
% Description:
%   This function checks the validity of the input parameters for the TQWT.
%
% Example:
%   check_params(Q, r, J);
%
% See also: (if applicable)

[st, ~] = dbstack(1);

%--------Check Q-----
if ~isscalar(Q)
    error('Error in %s: Q must be a scalar.\n', st.file)
end

if ~isnumeric(Q)
    error('Error in %s: Q must be numeric.\n', st.file)
end

if Q < 1
    error('Error in %s: Q must be greater than or equal to 1.0.\n', st.file)
end
