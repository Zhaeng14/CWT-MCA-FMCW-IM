function H1 = H1_fun(w, alpha, beta)

% H1 = H1_fun(w, alpha, beta)

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick,  selesi@poly.edu
% Polytechnic Institute of NYU
% November 2010

H1 = zeros(size(w));

w = mod(w+pi, 2*pi)-pi;

H1(abs(w) >= alpha*pi) = 1;

k = (abs(w) >= (1-beta)*pi) & (abs(w) <= alpha*pi);

H1(k) = theta_fun( (alpha*pi - w(k))/(alpha+beta-1) );

