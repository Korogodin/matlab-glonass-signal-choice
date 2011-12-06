clear 
close all
clc

path_to_results_back_inter = [pwd '/results/back_intersystem_L2'];
path_to_results_inter = [pwd '/results/intersystem_L2'];

load([path_to_results_back_inter '/BackInterSysJam_BoCsin_L2_BoC_10_5.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BoCcos_L2_BoC_10_5.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BPSK_L2_BoC_10_5.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BoCsin_L2_BoC_0_1.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BoCcos_L2_BoC_0_1.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BPSK_L2_BoC_0_1.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BoCsin_L2_BoC_0_10.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BoCcos_L2_BoC_0_10.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BPSK_L2_BoC_0_10.mat']);
load([path_to_results_inter '/InterSysJam_BoCsin_L2_GloST.mat'], 'InterSysJam_BoCsin_L2_GloST');
load([path_to_results_inter '/InterSysJam_BoCcos_L2_GloST.mat'], 'InterSysJam_BoCcos_L2_GloST');
load([path_to_results_inter '/InterSysJam_BoCsin_L2_GloVT.mat'], 'InterSysJam_BoCsin_L2_GloVT');
load([path_to_results_inter '/InterSysJam_BoCcos_L2_GloVT.mat'], 'InterSysJam_BoCcos_L2_GloVT');
load([path_to_results_inter '/InterSysJam_BPSK_L2_GloST.mat'], 'InterSysJam_BPSK_L2_GloST');
load([path_to_results_inter '/InterSysJam_BPSK_L2_GloVT.mat'], 'InterSysJam_BPSK_L2_GloVT');
load([path_to_results_inter '/InterSysJam_BoCsin_L2_GloST_mean.mat'], 'InterSysJam_BoCsin_L2_GloST_mean');
load([path_to_results_inter '/InterSysJam_BoCcos_L2_GloST_mean.mat'], 'InterSysJam_BoCcos_L2_GloST_mean');
load([path_to_results_inter '/InterSysJam_BoCsin_L2_GloVT_mean.mat'], 'InterSysJam_BoCsin_L2_GloVT_mean');
load([path_to_results_inter '/InterSysJam_BoCcos_L2_GloVT_mean.mat'], 'InterSysJam_BoCcos_L2_GloVT_mean');
load([path_to_results_inter '/InterSysJam_BPSK_L2_GloST_mean.mat'], 'InterSysJam_BPSK_L2_GloST_mean');
load([path_to_results_inter '/InterSysJam_BPSK_L2_GloVT_mean.mat'], 'InterSysJam_BPSK_L2_GloVT_mean');

k_JN0_GPS_L2_P = 0;
k_JN0_GPS_L2_M = 0;
k_JN0_GPS_L2_L2C = 0;
k_JN0_GLO_L2_L2OF = 0;
k_JN0_GLO_L2_L2SF = 0;

Ml_GPS = 1/4;
Ml_GLO = 1/4;

n8max = 80;
m8max = 80;
farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты

BackInterSysJam_BoCsin_L2 = nan(m8max, n8max, fmax);
BackInterSysJam_BoCcos_L2 = nan(m8max, n8max, fmax);
BackInterSysJam_BPSK_L2 = nan(n8max, fmax);
%     fprintf('|+ Intersystem jammer for BoCsin L2: GPS and GLONASS FDMA to GLONASS CDMA, dB\n');

for freq_index = 1:fmax
    fprintf('sin %f\n', freq_index);
    for n8 = 1:n8max
        for m8 = 1:m8max

            if m8 < n8
                    continue;
            end            
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (18-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(BackInterSysJam_BoCsin_L2_BoC_0_1(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + (10^((BackInterSysJam_BoCsin_L2_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(BackInterSysJam_BoCsin_L2_BoC_0_10(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L2_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L2_P/10);
                
            end
            if isnan(BackInterSysJam_BoCsin_L2_BoC_10_5(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L2_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L2_M/10);
            end
            if isnan(InterSysJam_BoCsin_L2_GloST_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L2_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L2_L2OF/10);
            end            
            if isnan(InterSysJam_BoCsin_L2_GloVT_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L2_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L2_L2SF/10);
            end                  
            if sum_dB == 0
                fprintf('NaN\n');
            else
                BackInterSysJam_BoCsin_L2(m8, n8, freq_index) = 10*log10(sum_dB);                
            end
        end
    end
end
save([path_to_results_back_inter '/Common/BackInterSysJam_BoCsin_L2.mat'], 'BackInterSysJam_BoCsin_L2');        


for freq_index = 1:fmax
    for n8 = 1:n8max
        for m8 = 1:m8max

            if m8 < n8
                    continue;
            end            
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (18-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(BackInterSysJam_BoCcos_L2_BoC_0_1(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + (10^((BackInterSysJam_BoCcos_L2_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(BackInterSysJam_BoCcos_L2_BoC_0_10(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L2_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L2_P/10);
                
            end
            if isnan(BackInterSysJam_BoCcos_L2_BoC_10_5(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L2_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L2_M/10);
            end
            if isnan(InterSysJam_BoCcos_L2_GloST_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L2_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L2_L2OF/10);
            end            
            if isnan(InterSysJam_BoCcos_L2_GloVT_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L2_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L2_L2SF/10);
            end                  
            if sum_dB == 0
                fprintf('NaN\n');
            else
                BackInterSysJam_BoCcos_L2(m8, n8, freq_index) = 10*log10(sum_dB);                
            end
        end
    end
end
save([path_to_results_back_inter '/Common/BackInterSysJam_BoCcos_L2.mat'], 'BackInterSysJam_BoCcos_L2');        


for freq_index = 1:fmax
    for n8 = 1:n8max
            m8 = 0;
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (18-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(BackInterSysJam_BPSK_L2_BoC_0_1(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + (10^((BackInterSysJam_BPSK_L2_BoC_0_1(n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(BackInterSysJam_BPSK_L2_BoC_0_10(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L2_BoC_0_10(n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L2_P/10);
                
            end
            if isnan(BackInterSysJam_BPSK_L2_BoC_10_5(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L2_BoC_10_5(n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L2_M/10);
            end
            if isnan(InterSysJam_BPSK_L2_GloST_mean(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L2_GloST_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L2_L2OF/10);
            end            
            if isnan(InterSysJam_BPSK_L2_GloVT_mean(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L2_GloVT_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L2_L2SF/10);
            end                  
            if sum_dB == 0
                fprintf('NaN\n');
            else
                BackInterSysJam_BPSK_L2(n8, freq_index) = 10*log10(sum_dB);                
            end
    end
end
save([path_to_results_back_inter '/Common/BackInterSysJam_BPSK_L2.mat'], 'BackInterSysJam_BPSK_L2');   