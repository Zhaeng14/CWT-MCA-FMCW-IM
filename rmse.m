clc;
clear;
close all;

%load("data/rmse_10_new.mat")   %Load the previous statistical dataset for plotting

load('data/500_10.mat')%Load simulation dataset
%FMCW radar parameters
Tr=1.2e-05; % Period
bw=360e6; % Bandwidth
slope=bw/Tr; % slope(Hz/s)
c0=3e+8; 

N = 512; % number of samples
Fs = N/Tr;% sampling frequency

t=0:1/Fs:Tr-1/Fs; %intervalul de timp(s) al semnalului transmis

Nfft = N;% number of FFT points

F = (0:1/Nfft:1-1/Nfft)*Fs; % frequency vector
r_axis = c0/2/slope*F;
%%
SNR = -20:2:20;
rmse_sb = zeros(1,length(SNR));
rmse_sb0 = zeros(1,length(SNR));
rmse_lpc = zeros(1,length(SNR));
rmse_mca = zeros(1,length(SNR));
rmse_mca_2 = zeros(1,length(SNR));
rmse_mca_3 = zeros(1,length(SNR));
rmse_mca_4 = zeros(1,length(SNR));


sinr_sb = zeros(1,length(SNR));
sinr_sb0 = zeros(1,length(SNR));
sinr_lpc = zeros(1,length(SNR));
sinr_mca = zeros(1,length(SNR));
sinr_mca2 = zeros(1,length(SNR));
sinr_mca3 = zeros(1,length(SNR));
sinr_mca4 = zeros(1,length(SNR));


cdf_sb = zeros(1,length(SNR));
cdf_sb0 = zeros(1,length(SNR));
cdf_lpc = zeros(1,length(SNR));
cdf_mca = zeros(1,length(SNR));
cdf_mca_2 = zeros(1,length(SNR));
cdf_mca_3 = zeros(1,length(SNR));
cdf_mca_4 = zeros(1,length(SNR));


sir = 5;
interfer_coef = 0.5;
mtkl = 500;
[~,len] = size(sb_mat);


for i_SNR = 1:length(SNR)
    for i_mtkl = 1:mtkl
        i_ind= i_mtkl + (i_SNR-1)*mtkl;

        mca_prediction_4 = mca_method_4((sb_mat(i_ind,:))) /len;
        [~, mca_index_4] = max(mca_prediction_4);


        target_index = find(abs(amplitude_mat(i_ind,:)) > 0);
        target = r_axis(target_index);
        target = r_true(i_ind);

        mca_prediction_3 = mca_method_3((sb_mat(i_ind,:))) /len;
        [~, mca_index_3] = max(mca_prediction_3);


        mca_prediction_2 = mca_method_2((sb_mat(i_ind,:)).') /len;
        [~, mca_index_2] = max(mca_prediction_2);


        mca_prediction = mca_method((sb_mat(i_ind,:)).')/len;
        mca_abs = abs(mca_prediction);
        [~, mca_index] = max(mca_prediction);


        sb0_fft = abs(fft(sb0_mat(i_ind, :)))/len;
        [~, sb0_index] = max(sb0_fft);

        sb_fft = abs(fft(sb_mat(i_ind, :)))/len;
        [~, sb_index] = max(sb_fft);

        linear_prediction = abs(fft(my_lpc(sb_mat(i_ind, :), 3)))/len; 
        [~, lpc_index] = max(linear_prediction);

        rmse_sb(i_SNR) = rmse_sb(i_SNR) + (abs(r_axis(sb_index) - target))^2;
        rmse_sb0(i_SNR) = rmse_sb0(i_SNR) + (abs(r_axis(sb0_index) - target))^2;
        rmse_lpc(i_SNR) = rmse_lpc(i_SNR) + (abs(r_axis(lpc_index) - target))^2;
        rmse_mca(i_SNR) = rmse_mca(i_SNR) + (abs(r_axis(mca_index) - target))^2;
        rmse_mca_2(i_SNR) = rmse_mca_2(i_SNR) + (abs(r_axis(mca_index_2) - target))^2;
        rmse_mca_3(i_SNR) = rmse_mca_3(i_SNR) + (abs(r_axis(mca_index_3) - target))^2;
        rmse_mca_4(i_SNR) = rmse_mca_4(i_SNR) + (abs(r_axis(mca_index_4) - target))^2;



        sinr_sb(i_SNR) = sinr_sb(i_SNR) + mean(sb_fft(target_index).^2) / mean(sb_fft.^2);
        sinr_sb0(i_SNR) = sinr_sb0(i_SNR) + mean(sb0_fft(target_index).^2) / mean(sb0_fft.^2);
        sinr_lpc(i_SNR) = sinr_lpc(i_SNR) + mean(linear_prediction(target_index).^2) / mean(linear_prediction.^2);
        sinr_mca(i_SNR) = sinr_mca(i_SNR) + mean(mca_prediction(target_index).^2) / mean(mca_prediction.^2);
        sinr_mca2(i_SNR) = sinr_mca2(i_SNR) + mean(mca_prediction_2(target_index).^2) / (mean(mca_prediction_2.^2)+1e-8);
        sinr_mca3(i_SNR) = sinr_mca3(i_SNR) + mean(mca_prediction_3(target_index).^2) / mean(mca_prediction_3.^2);
        sinr_mca4(i_SNR) = sinr_mca4(i_SNR) + mean(mca_prediction_4(target_index).^2) / mean(mca_prediction_4.^2);



        if(target_index == sb0_index)
            cdf_sb0(i_SNR) = cdf_sb0(i_SNR) + 1;
        end
        if(target_index == sb_index)
            cdf_sb(i_SNR) = cdf_sb(i_SNR) + 1;
        end 
        if(abs(target_index - lpc_index) < 3)
            cdf_lpc(i_SNR) = cdf_lpc(i_SNR) + 1;
        end
        if(abs(target_index - mca_index) < 3)
            cdf_mca(i_SNR) = cdf_mca(i_SNR) + 1;
        end
        if(abs(target_index - mca_index_2) < 3)
            cdf_mca_2(i_SNR) = cdf_mca_2(i_SNR) + 1;
        end
        if(abs(target_index - mca_index_3) < 3)
            cdf_mca_3(i_SNR) = cdf_mca_3(i_SNR) + 1;
        end
        if(abs(target_index - mca_index_4) < 3)
            cdf_mca_4(i_SNR) = cdf_mca_4(i_SNR) + 1;
        end






    end
end
rmse_sb0 = sqrt(rmse_sb0/mtkl);
rmse_sb = sqrt(rmse_sb/mtkl);
rmse_lpc = sqrt(rmse_lpc/mtkl);
rmse_mca = sqrt(rmse_mca/mtkl);
rmse_mca_2 = sqrt(rmse_mca_2/mtkl);
rmse_mca_3 = sqrt(rmse_mca_3/mtkl);
rmse_mca_4 = sqrt(rmse_mca_4/mtkl);



sinr_sb0 = 10*log10(sinr_sb0/mtkl);
sinr_sb = 10*log10(sinr_sb/mtkl);
sinr_lpc = 10*log10(sinr_lpc/mtkl);
sinr_mca = 10*log10(sinr_mca/mtkl);
sinr_mca2 = 10*log10(sinr_mca2/mtkl);
sinr_mca3 = 10*log10(sinr_mca3/mtkl);
sinr_mca4 = 10*log10(sinr_mca4/mtkl);




cdf_sb0 = cdf_sb0/mtkl;
cdf_sb = cdf_sb/mtkl;
cdf_lpc = cdf_lpc/mtkl;
cdf_mca = cdf_mca/mtkl;
cdf_mca_2 = cdf_mca_2/mtkl;
cdf_mca_3 = cdf_mca_3/mtkl;
cdf_mca_4 = cdf_mca_4/mtkl;



    figure(1)
    semilogy(SNR, rmse_sb0, 'linewidth', 1.5);hold on;
    % semilogy(SNR, rmse_sb, 'linewidth', 1.5);hold on;
    semilogy(SNR, rmse_lpc, 'linewidth', 1.5);hold on;
    semilogy(SNR, rmse_mca, 'linewidth', 1.5);grid on;
    semilogy(SNR, rmse_mca_2, 'linewidth', 1.5);grid on;
    semilogy(SNR, rmse_mca_3, 'linewidth', 1.5);grid on;
    semilogy(SNR, rmse_mca_4, 'linewidth', 1.5);grid on;



    legend('clear','LPC', 'DFT-STFT-MCA','MOM-MCA','TQWT-MCA','CWT-MCA');
    % axis([-20 20 1e-5 1e2]);
    xlabel('SNR(dB)');ylabel('RMSE(m)');
    %     title('sir=5 dB');


    figure(2)
          plot(SNR, sinr_sb, 'linewidth', 1.5);hold on;
          plot(SNR, sinr_mca, 'linewidth', 1.5);grid on;
          plot(SNR, sinr_lpc, 'linewidth', 1.5);grid on;
          plot(SNR, sinr_mca2, 'linewidth', 1.5);grid on;
          plot(SNR, sinr_mca3, 'linewidth', 1.5);grid on;
          plot(SNR, sinr_mca4, 'linewidth', 1.5);grid on;


          legend('BEFORE IM', 'LPC','DFT-STFT-MCA','MOM-MCA','TQWT-MCA','CWT-MCA','location','southeast');
          xlabel('SNR(dB)');ylabel('SINR(dB)');




    figure(3)
    plot(SNR, cdf_sb0, 'linewidth', 2);hold on;
    % plot(SNR, cdf_sb, 'linewidth', 2);hold on;
    plot(SNR, cdf_lpc, 'linewidth', 2);hold on;
    plot(SNR, cdf_mca, 'linewidth', 2);grid on;
    plot(SNR, cdf_mca_2, 'linewidth', 2);grid on;
    plot(SNR, cdf_mca_3, 'linewidth', 2);grid on;
    plot(SNR, cdf_mca_4, 'linewidth', 2);grid on;
    legend('CLEAR', 'LPC', 'DFT-STFT-MCA', 'MOM-MCA','TQWT-MCA','CWT-MCA');
    xlabel('SNR(dB)');ylabel('Detection Probability');

  % save('data/rmse_10_new.mat', 'rmse_sb0', 'rmse_lpc', 'rmse_mca', 'rmse_mca_2', "rmse_mca_3", "rmse_mca_4");
  % save('data/100_10_z.mat', 'rmse_sb0','rmse_sb','rmse_lpc', 'rmse_zero', 'rmse_IMAT','rmse_mca','rmse_mca_2','rmse_mca_3',"rmse_mca_4");
