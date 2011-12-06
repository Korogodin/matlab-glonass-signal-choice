clear 
close all
clc

path_to_results_back_inter = [pwd '/results/back_intersystem_L3'];
path_to_results_inter = [pwd '/results/intersystem_L3'];

load([path_to_results_back_inter '/BackInterSysJam_BoCsin_L3_BoC_0_10.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BoCcos_L3_BoC_0_10.mat']);
load([path_to_results_back_inter '/BackInterSysJam_BPSK_L3_BoC_0_10.mat']);

k_JN0_GPS_L3_P = 0;
k_JN0_GPS_L3_M = 0;
k_JN0_GPS_L3_L3C = 0;
k_JN0_GLO_L3_L3OF = 0;
k_JN0_GLO_L3_L3SF = 0;

Ml_GPS = 1/1;
Ml_GLO = 1/1;

n8max = 80;
m8max = 80;
farr = 1164:1184; fmax = length(farr); % Нормированный центральные частоты
Signal_Band = 20;

BackInterSysJam_BoCsin_L3 = nan(m8max, n8max, fmax);
BackInterSysJam_BoCcos_L3 = nan(m8max, n8max, fmax);
BackInterSysJam_BPSK_L3 = nan(n8max, fmax);
%     fprintf('|+ Intersystem jammer for BoCsin L3: GPS and GLONASS FDMA to GLONASS CDMA, dB\n');

% BOC sin
for freq_index = 1:fmax
    fprintf('sin %f\n', freq_index);
    for n8 = 1:n8max
        for m8 = 1:m8max

            if m8 < n8
                continue;
            end            
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (Signal_Band-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(BackInterSysJam_BoCsin_L3_BoC_0_10(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L3_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L3_P/10);
            end
            
            if sum_dB == 0
                fprintf('NaN\n');
            else
                BackInterSysJam_BoCsin_L3(m8, n8, freq_index) = 10*log10(sum_dB);                
            end
        end
    end
end
save([path_to_results_back_inter '/Common/BackInterSysJam_BoCsin_L3.mat'], 'BackInterSysJam_BoCsin_L3');        


for freq_index = 1:fmax
    for n8 = 1:n8max
        for m8 = 1:m8max

            if m8 < n8
                    continue;
            end            
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (Signal_Band-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(BackInterSysJam_BoCcos_L3_BoC_0_10(m8, n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L3_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L3_P/10);
            end
            
            if sum_dB == 0
                fprintf('NaN\n');
            else
                BackInterSysJam_BoCcos_L3(m8, n8, freq_index) = 10*log10(sum_dB);                
            end
        end
    end
end
save([path_to_results_back_inter '/Common/BackInterSysJam_BoCcos_L3.mat'], 'BackInterSysJam_BoCcos_L3');        


for freq_index = 1:fmax
    for n8 = 1:n8max
            m8 = 0;
            
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (Signal_Band-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            
            sum_dB = 0;

            if isnan(BackInterSysJam_BPSK_L3_BoC_0_10(n8, freq_index))
                fprintf('NaN\n');
            else
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L3_BoC_0_10(n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L3_P/10);
            end
            
            if sum_dB == 0
                fprintf('NaN\n');
            else
                BackInterSysJam_BPSK_L3(n8, freq_index) = 10*log10(sum_dB);                
            end
    end
end
save([path_to_results_back_inter '/Common/BackInterSysJam_BPSK_L3.mat'], 'BackInterSysJam_BPSK_L3');   