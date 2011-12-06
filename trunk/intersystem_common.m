clear 
close all
clc

path_to_results_inter = [pwd '/results/intersystem_L1'];

load([path_to_results_inter '/InterSysJam_BoCsin_L1_BoC_1_1.mat']); 
load([path_to_results_inter '/InterSysJam_BoCcos_L1_BoC_1_1.mat']);
load([path_to_results_inter '/InterSysJam_BPSK_L1_BoC_1_1.mat']);
load([path_to_results_inter '/InterSysJam_BoCsin_L1_BoC_6_1.mat']); 
load([path_to_results_inter '/InterSysJam_BoCcos_L1_BoC_6_1.mat']);
load([path_to_results_inter '/InterSysJam_BPSK_L1_BoC_6_1.mat']);
load([path_to_results_inter '/InterSysJam_BoCsin_L1_BoC_10_5.mat']);
load([path_to_results_inter '/InterSysJam_BoCcos_L1_BoC_10_5.mat']);
load([path_to_results_inter '/InterSysJam_BPSK_L1_BoC_10_5.mat']);
load([path_to_results_inter '/InterSysJam_BoCsin_L1_BoC_0_1.mat']);
load([path_to_results_inter '/InterSysJam_BoCcos_L1_BoC_0_1.mat']);
load([path_to_results_inter '/InterSysJam_BPSK_L1_BoC_0_1.mat']);
load([path_to_results_inter '/InterSysJam_BoCsin_L1_BoC_0_10.mat']);
load([path_to_results_inter '/InterSysJam_BoCcos_L1_BoC_0_10.mat']);
load([path_to_results_inter '/InterSysJam_BPSK_L1_BoC_0_10.mat']);
load([path_to_results_inter '/InterSysJam_BoCsin_L1_GloST.mat'], 'InterSysJam_BoCsin_L1_GloST');
load([path_to_results_inter '/InterSysJam_BoCcos_L1_GloST.mat'], 'InterSysJam_BoCcos_L1_GloST');
load([path_to_results_inter '/InterSysJam_BoCsin_L1_GloVT.mat'], 'InterSysJam_BoCsin_L1_GloVT');
load([path_to_results_inter '/InterSysJam_BoCcos_L1_GloVT.mat'], 'InterSysJam_BoCcos_L1_GloVT');
load([path_to_results_inter '/InterSysJam_BPSK_L1_GloST.mat'], 'InterSysJam_BPSK_L1_GloST');
load([path_to_results_inter '/InterSysJam_BPSK_L1_GloVT.mat'], 'InterSysJam_BPSK_L1_GloVT');
load([path_to_results_inter '/InterSysJam_BoCsin_L1_GloST_mean.mat'], 'InterSysJam_BoCsin_L1_GloST_mean');
load([path_to_results_inter '/InterSysJam_BoCcos_L1_GloST_mean.mat'], 'InterSysJam_BoCcos_L1_GloST_mean');
load([path_to_results_inter '/InterSysJam_BoCsin_L1_GloVT_mean.mat'], 'InterSysJam_BoCsin_L1_GloVT_mean');
load([path_to_results_inter '/InterSysJam_BoCcos_L1_GloVT_mean.mat'], 'InterSysJam_BoCcos_L1_GloVT_mean');
load([path_to_results_inter '/InterSysJam_BPSK_L1_GloST_mean.mat'], 'InterSysJam_BPSK_L1_GloST_mean');
load([path_to_results_inter '/InterSysJam_BPSK_L1_GloVT_mean.mat'], 'InterSysJam_BPSK_L1_GloVT_mean');

k_JN0_GPS_L1_P = -3;
k_JN0_GPS_L1_M = +0.5;
k_JN0_GPS_L1_L1C = +1.5;
k_JN0_GLO_L1_L1OF = -2.5;
k_JN0_GLO_L1_L1SF = -2.5;

Ml_GPS = 10;
Ml_GLO = 8;

n8max = 80;
m8max = 80;
farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты

InterSysJam_BoCsin_L1 = nan(m8max, n8max, fmax);
InterSysJam_BoCcos_L1 = nan(m8max, n8max, fmax);
InterSysJam_BPSK_L1 = nan(n8max, fmax);
%     fprintf('|+ Intersystem jammer for BoCsin L1: GPS and GLONASS FDMA to GLONASS CDMA, dB\n');

for freq_index = 1:fmax
    fprintf('sin %f\n', freq_index);
    for n8 = 1:n8max
        for m8 = 1:m8max

            if m8 < n8
                    continue;
            end            
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (16-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + (10^((InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('NaN\n');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BoCsin_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if isnan(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end                  
            if sum_dB == 0
                fprintf('NaN\n');
            else
                InterSysJam_BoCsin_L1(m8, n8, freq_index) = 10*log10(sum_dB);                
            end
        end
    end
end
save([path_to_results_inter '/Common/InterSysJam_BoCsin_L1.mat'], 'InterSysJam_BoCsin_L1');        


for freq_index = 1:fmax
    for n8 = 1:n8max
        for m8 = 1:m8max

            if m8 < n8
                    continue;
            end            
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (16-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + (10^((InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('NaN\n');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BoCcos_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if isnan(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end                  
            if sum_dB == 0
                fprintf('NaN\n');
            else
                InterSysJam_BoCcos_L1(m8, n8, freq_index) = 10*log10(sum_dB);                
            end
        end
    end
end
save([path_to_results_inter '/Common/InterSysJam_BoCcos_L1.mat'], 'InterSysJam_BoCcos_L1');        


for freq_index = 1:fmax
    for n8 = 1:n8max
            m8 = 0;
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (16-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + (10^((InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))
                fprintf('NaN\n');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BPSK_L1_BoC_6_1(n8, freq_index))/10);
                
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if isnan(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end                  
            if sum_dB == 0
                fprintf('NaN\n');
            else
                InterSysJam_BPSK_L1(n8, freq_index) = 10*log10(sum_dB);                
            end
    end
end
save([path_to_results_inter '/Common/InterSysJam_BPSK_L1.mat'], 'InterSysJam_BPSK_L1');   