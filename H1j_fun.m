function H1 = H1j_fun(w, alpha, beta, j)

% H1 = H1j_fun(w, alpha, beta, j)
%
% j : non-netative integer

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick,  selesi@poly.edu
% Polytechnic Institute of NYU
% November 2010

H1 = H0j_fun(w, alpha, beta, j-1) .* H1_fun(w / alpha^(j-1), alpha, beta);
