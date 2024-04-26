function [x1,x2,w1,w2,costfn] = tqwt_mca(x,Q1,r1,J1,Q2,r2,J2,lam1,lam2,mu,Nit,plot_flag)
% [x1,x2,w1,w2] = tqwt_mca(x,Q1,r1,J1,Q2,r2,J2,lam1,lam2,mu,Nit)
% TQWT-MCA decomposition using SALSA
% INPUT
%   x - input signal signal
%   Q1, r1, J1, Q1, r1, J2 - TQWT parameters
%   lam1, lam2 - regularization parameters
%   mu - SALSA parameter
%   Nit - Number of iterations
% OUTPUT
%   x1, x2 - components
%   w1, w2 - transform coefficients of components
%
% Use [x1,x2,w1,w2,costfn] = tqwt_mca(...) to return cost function.
%
% Use [...] = tqwt_mca(...,'plots') to plot progress of algorithm.


RTV = true;       % Return solution after thresholding step
% RTV = false;        % Return solution thresholding after linear projection step
% Note: cost function may be slightly different. Algorithm is simpler
% for RTV = false;

% By default do not compute cost function (to save computation)
if nargout > 4
    COST = true;
    costfn = zeros(1,Nit);     % cost function
else
    COST = false;
    costfn = [];
end

GOPLOTS = false;
if nargin == 12
    if strcmp(plot_flag,'plots')
        GOPLOTS = true;
    end
end

N = length(x);

% Initialize:
w1 = tqwt_radix2(x,Q1,r1,J1);
w2 = tqwt_radix2(x,Q2,r2,J2);
d1 = tqwt_radix2(zeros(size(x)),Q1,r1,J1);
d2 = tqwt_radix2(zeros(size(x)),Q2,r2,J2);


T1 = lam1/(2*mu);
T2 = lam2/(2*mu);

L1 = length(lam1);
L2 = length(lam2);

u1 = cell(1,L1);
u2 = cell(1,L2);

v1 = cell(1,L1);
v2 = cell(1,L2);

N = length(x);
A = 1.1*max(abs(x));

for k = 1:Nit
    
    for j = 1:L1
        u1{j} = soft(w1{j} + d1{j}, T1(j)) - d1{j};
    end
    for j = 1:L2
        u2{j} = soft(w2{j} + d2{j}, T2(j)) - d2{j};
    end
    
    c = x - itqwt_radix2(u1,Q1,r1,N) - itqwt_radix2(u2,Q2,r2,N);
    c = c/(mu+2);
    
    d1 = tqwt_radix2(c,Q1,r1,J1);
    d2 = tqwt_radix2(c,Q2,r2,J2);
    
    for j = 1:L1
        w1{j} = d1{j} + u1{j};
    end
    
    for j = 1:L2
        w2{j} = d2{j} + u2{j};
    end
    
    
    if COST || GOPLOTS
        
        if RTV
            for j = 1:L1
                v1{j} = soft(w1{j} + d1{j}, T1(j));
            end
            for j = 1:L2
                v2{j} = soft(w2{j} + d2{j}, T2(j));
            end
            x1 = itqwt_radix2(v1,Q1,r1,N);
            x2 = itqwt_radix2(v2,Q2,r2,N);
            
            res = x - x1 - x2;
            costfn(k) = sum(abs(res).^2);
            for j = 1:L1
                costfn(k) = costfn(k) + lam1(j)*sum(abs(v1{j}));
            end
            for j = 1:L2
                costfn(k) = costfn(k) + lam2(j)*sum(abs(v2{j}));
            end
            
        else
            x1 = itqwt_radix2(w1,Q1,r1,N);
            x2 = itqwt_radix2(w2,Q2,r2,N);
            
            res = x - x1 - x2;
            costfn(k) = sum(abs(res).^2);
            for j = 1:L1
                costfn(k) = costfn(k) + lam1(j)*sum(abs(w1{j}));
            end
            for j = 1:L2
                costfn(k) = costfn(k) + lam2(j)*sum(abs(w2{j}));
            end
            
        end
        
    end
    
    if GOPLOTS
        figure(gcf)
        clf
        subplot(3,1,1)
        plot(real(x1))
        xlim([0 N])
        ylim([-A A])
        title({sprintf('ITERATION %d',k),'COMPONENT 1'})
        box off
        subplot(3,1,2)
        plot(real(x2))
        xlim([0 N])
        ylim([-A A])
        box off
        title('COMPONENT 2')
        subplot(3,1,3)
        plot(real(res))
        xlim([0 N])
        ylim([-A A])
        title('RESIDUAL')
        box off
        drawnow
    end
    
end

if RTV
    for j = 1:L1
        v1{j} = soft(w1{j} + d1{j}, T1(j));
    end
    for j = 1:L2
        v2{j} = soft(w2{j} + d2{j}, T2(j));
    end
    w1 = v1;
    w2 = v2;
    x1 = itqwt_radix2(w1,Q1,r1,N);
    x2 = itqwt_radix2(w2,Q2,r2,N);
end
end
% --------------- local function: soft ---------------

function y = soft(x,T)
% Soft-threshold function
% y = soft_fun(x,T)
% x : input data
% T : threshold

y = zeros(size(x));
k = (x < -T);
y(k) = x(k) + T;
k = (x > T);
y(k) = x(k) - T;

% following alternative definition works for real and complex data:
g = max(abs(x)-T,0);
y = g./(g+T) .* x;
end
% function result = soft(y, T)
%     result = wthresh(y,'s',T);
% end

