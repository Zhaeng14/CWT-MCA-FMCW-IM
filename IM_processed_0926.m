clc
clear all
close all

p0 = [];
p123_1 = [];
p456_1 = [];
sb0 = [];
sb = [];

FFT_LEN = 512;
CHIRP_NUM = 512;
fidimag = fopen('Rx0_IFFT_imag0824_RX1_145244_float.bin','rb');
[data_imag, count_imag] = fread(fidimag,'float');
fidreal = fopen('Rx0_IFFT_real0824_RX1_145244_float.bin','rb');;
[data_real, count_real] = fread(fidreal,'float');
fclose(fidimag);
fclose(fidreal);

frame_cnt = count_real/FFT_LEN/CHIRP_NUM;

 % for framenum = 1:frame_cnt
 for framenum = 450

    datablock_startIdx = CHIRP_NUM*FFT_LEN*(framenum-1);

    for chirpnum=1:512
        for rangebin = 1:512           
            data(chirpnum,rangebin) = data_real(datablock_startIdx + (chirpnum-1)*FFT_LEN + rangebin)+sqrt(-1)*data_imag(datablock_startIdx + (chirpnum-1)*FFT_LEN + rangebin);    
        end

             data2(chirpnum,:) = data(chirpnum,:);  
             data_abs(chirpnum,:) = abs(data2(chirpnum,:));

%          data2(chirpnum,:) = add_interference(data(chirpnum,:), -15);

           % data_processed_by_lpc(chirpnum,:) = my_lpc(data2(chirpnum,:),3); 
           data_processed_by_dft(chirpnum,:) = mca_method(data2(chirpnum,:).');
           % data_processed_by_mom(chirpnum,:) = mca_method_2(data2(chirpnum,:).');
           % [data_processed_by_tqwt(chirpnum,:),data_processed_by_tqwt2(chirpnum,:)] = mca_method_3(data2(chirpnum,:)); 
           % data_processed_by_mca_cwt(chirpnum,:) = mca_method_4(data2(chirpnum,:));


           % dataRx0Tx1_fft(chirpnum,:) = fft(data2(chirpnum,:)); 
           % dataRx0Tx1_fft(chirpnum,:) = fft(data_processed_by_lpc(chirpnum,:));           
           dataRx0Tx1_fft(chirpnum,:) = fft(data_processed_by_dft(chirpnum,:)); 
           % dataRx0Tx1_fft(chirpnum,:) = fft(data_processed_by_mom(chirpnum,:)); 
           % dataRx0Tx1_fft(chirpnum,:) = fft(data_processed_by_tqwt(chirpnum,:));  dataRx0Tx1_fft2(chirpnum,:) = fft(data_processed_by_tqwt2(chirpnum,:));
           % dataRx0Tx1_fft(chirpnum,:) = fft(data_processed_by_mca_cwt(chirpnum,:)); 
    end
      
    for rangebin = 1:128
        test_data = dataRx0Tx1_fft(:,rangebin);
        dataRx0Tx1_2dfft(:,rangebin) = fft(test_data);

        % test_data2 = dataRx0Tx1_fft2(:,rangebin);
        % dataRx0Tx1_2dfft2(:,rangebin) = fft(test_data2);
    end
    
    dataRx0Tx1_2dfft_abs = abs(dataRx0Tx1_2dfft).*abs(dataRx0Tx1_2dfft);
    dataRx0Tx1_2dfft_abs = (mag2db(dataRx0Tx1_2dfft_abs));
    dataRx0Tx1_2dfft_abs_y = flipud(dataRx0Tx1_2dfft_abs);
    dataRx0Tx1_2dfft_abs_x = fliplr(dataRx0Tx1_2dfft_abs);
    dataRx0Tx1_2dfft_abs_mirror = rot90(dataRx0Tx1_2dfft_abs, 2);
    dataRx0Tx1_2dfft_abs_reshape = permute(dataRx0Tx1_2dfft_abs, [3, 1, 2]);

    cycleCount = framenum;
    fig = figure(123);
    if isempty(p123_1)  

            p123_1 = imagesc(dataRx0Tx1_2dfft_abs,[0,180]); 
        set(gca,'ydir','normal')
            grid on
            xlabel('range bin []')
            ylabel('kappa bin []')
            title(['Kappa-Image - cycleCount: ' num2str(cycleCount)]) 
            colorbar   
    else
        set(p123_1, 'CData', dataRx0Tx1_2dfft_abs);
        title(['RV-Image - cycleCount: ' num2str(cycleCount)]) 
    end


%%无干扰数据集
    sb0(framenum,:,:) = dataRx0Tx1_2dfft_abs;
    sb0(framenum + 300,:,:) = dataRx0Tx1_2dfft_abs_y;
    sb0(framenum + 600,:,:) = dataRx0Tx1_2dfft_abs_x;    
    sb0(framenum + 900,:,:) = dataRx0Tx1_2dfft_abs_mirror; 
 end




