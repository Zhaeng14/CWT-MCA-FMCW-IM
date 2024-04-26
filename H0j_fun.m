function H0 = H0j_fun(w, alpha, beta, j)
% H0 = H0j_fun(w, alpha, beta, j)
% 0 < alpha, beta < 1
% j : non-netative integer

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick,  selesi@poly.edu
% Polytechnic Institute of NYU
% November 2010

H0 = ones(size(w));
for m = 0:j-1
    H0 = H0 .* H0_fun(w / alpha^m, alpha, beta);
end


% H0 = ones(size(w));
% for k = 1:j
%     H0 = H0 .* H0_fun(w / alpha^(j-k), alpha, beta);
% end


% H = H_fun(w, alpha, beta);
% 
% for k = 1:j-1
%     H = H .* H_fun(w/alpha^k, alpha, beta);
% end


% H = (abs(w) <= pi * alpha^j);
%
% for k = 1:j
%     H = H .* H_fun(alpha^(k-j)*w, alpha, beta);
% end
