function [x1,x2,w1,w2,costfn] = cwt_mca(x,dt,pad,dj1,dj2,lam1,lam2,mu,Nit,plot_flag)
% CWT-MCA (Continuous Wavelet Transform - Morphological Component Analysis) function
% Separates components from mixed signal based on the sparsity differences
% in different Continuous Wavelet Transform (CWT) domains.
%
% INPUT
%   x - Input mixed FMCW radar signal
%   dt - Time step
%   pad - Padding
%   dj1, dj2 - Dj parameters for CWT
%   lam1, lam2 - Regularization parameters
%   mu - SALSA (Sparse Analysis and Synthesis) parameter
%   Nit - Number of iterations
%   plot_flag - Flag for plotting intermediate results ('donotplots' to disable)
%
% OUTPUT
%   x1, x2 - Processed components obtained from CWT-MCA
%   w1, w2 - CWT coefficients corresponding to x1 and x2
%   costfn - Cost function values at each iteration

RTV = true;       % Return solution after thresholding step

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
% dt = 1.2e-05/256;;
% pad = 0;
% dj1 = 0.1; %smaller number gives better resolution, default = 0.25;
% dj2 = 1.0;
so =  2 * dt ; %default dt
Jfac = 1.0; %Sets the maximum scale to compute at (therefore number of scales). 1 is equiv to default. 
j1 =  round(Jfac*(log2(N*dt/so))/dj1); %default: (log2(N*dt/so))/dj
j2 =  round(Jfac*(log2(N*dt/so))/dj2);

% alpha1 = 1/2;
% alpha2 = 2/3;
% Alpha parameters for thresholding
alpha1 = 0.3;%真实数据0.45   rmse 0.1and 0.4
alpha2 = 0.6;%真实数据0.805
% Morlet wavelet is used as the mother wavelet
mother = 'morlet';
param = 5; %wave number for morlet, see >> help wave_bases for more details   真实数据8
N = length(x);
% Initialize CWT coefficients
[w1, period, scale, coi, dj1, paramout, k1] = contwt(x,dt, pad, dj1, so, j1, mother, param,1/alpha1);
[w2, period2, scale2, coi2, dj2, paramout, k2] = contwt(x,dt, pad, dj2, so, j2, mother, param,1/alpha2);

[L1,~] = size(w1);
[L2,~] = size(w2);
% Initialize auxiliary variables
[d1, period, scale, coi, dj1, paramout, k1] = contwt(zeros(1,N),dt, pad, dj1, so, j1, mother, param,1/alpha1);
[d2, period2, scale2, coi2, dj2, paramout, k2] = contwt(zeros(1,N),dt, pad, dj2, so, j2, mother, param,1/alpha2);

% Thresholding parameters
T1 = lam1/(2*mu);
T2 = lam2/(2*mu);

u1 = zeros(L1,N);
u2 = zeros(L2,N);

v1 = zeros(L1,N);
v2 = zeros(L2,N);

% N = length(x);
A = 1.1*max(abs(x));

for k = 1:Nit
    
    for j = 1:L1
            % Soft-thresholding for the first component
        u1(j,:) = soft(w1(j,:) + d1(j,:), T1) - d1(j,:);
    end
    for j = 1:L2
           % Soft-thresholding for the second component
        u2(j,:) = soft(w2(j,:) + d2(j,:), T2) - d2(j,:);
    end
        % Update residuals
    c = x - invcwt(u1, mother, scale, param, k1) - invcwt(u2, mother, scale2, param, k2);
    c = c/(mu+2);
      % Compute CWT of the residuals
    [d1, period, scale, coi, dj1, paramout, k1] = contwt(c,dt, pad, dj1, so, j1, mother, param,1/alpha1);
    [d2, period2, scale2, coi2, dj2, paramout, k2] = contwt(c,dt, pad, dj2, so, j2, mother, param,1/alpha2);
      % Update CWT coefficients
    for j = 1:L1
        w1(j,:) = d1(j,:) + u1(j,:);
    end
    
    for j = 1:L2
        w2(j,:) = d2(j,:) + u2(j,:);
    end
    
       % Compute cost function if needed
    if COST || GOPLOTS
        
        if RTV
            for j = 1:L1
                v1(j,:) = soft(w1(j,:) + d1(j,:), T1);
            end
            for j = 1:L2
                v2(j,:) = soft(w2(j,:) + d2(j,:), T2);
            end
            x1 = invcwt(w1, mother, scale, param, k1);
            x2 = invcwt(w2, mother, scale2, param, k2);
            
            res = x - x1 - x2;
            % costfn(k) = sum(abs(res).^2);
            % for j = 1:L1
            %     costfn(k) = costfn(k) + lam1*sum(abs(v1(j,:)));
            % end
            % for j = 1:L2
            %     costfn(k) = costfn(k) + lam2*sum(abs(v2(j,:)));
            % end
            
        else
            x1 = invcwt(w1, mother, scale, param, k1);
            x2 = invcwt(w2, mother, scale2, param, k2);
            
            res = x - x1 - x2;
            % costfn(k) = sum(abs(res).^2);
            % for j = 1:L1
            %     costfn(k) = costfn(k) + lam1*sum(abs(w1(j,:)));
            % end
            % for j = 1:L2
            %     costfn(k) = costfn(k) + lam2*sum(abs(w2(j,:)));
            % end
            
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
        v1(j,:) = soft(w1(j,:) + d1(j,:), T1);
    end
    for j = 1:L2
        v2(j,:) = soft(w2(j,:) + d2(j,:), T2);
    end
    w1 = v1;
    w2 = v2;
            x1 = invcwt(w1, mother, scale, param, k1);
            x2 = invcwt(w2, mother, scale2, param, k2);
end
end
% --------------- local function: soft ---------------

% function y = soft(x,T)
% % Soft-threshold function
% % y = soft_fun(x,T)
% % x : input data
% % T : threshold
% 
% y = zeros(size(x));
% k = (x < -T);
% y(k) = x(k) + T;
% k = (x > T);
% y(k) = x(k) - T;
% 
% % following alternative definition works for real and complex data:
% % g = max(abs(x)-T,0);
% % y = g./(g+T) .* x;
% end
% function result = soft(y, T)
%     result = wthresh(y,'s',T);
% end

function y = soft(x, T)
% 硬阈值函数
% y = hard(x, T)
% x : 输入数据
% T : 阈值

y = x;                  % 创建输出变量y并初始化为输入数据x
y(abs(x) <= T) = 0;     % 将绝对值小于等于阈值T的元素归零
end

