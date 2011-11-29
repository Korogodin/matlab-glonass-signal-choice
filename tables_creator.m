clear 
close all
clc

Table_num = 3.5;
% 1 - BoCsin intersys GPS
% 1.1 - BoCsin intersys GLONASS FDMA
% 1.2 - BoCsin intersys GPS+GLONASS FDMA
% 1.3 - BoCsin Back intersys GPS 
% 1.4 - BoCsin Back intersys GLO FDMA 
% 1.5 - BoCsin Back intersys GPS+GLO FDMA
% 2 - BoCcos intersys GPS
% 2.1 - BoCcos intersys GLONASS FDMA
% 2.2 - BoCcos intersys GPS+GLONASS FDMA
% 2.3 - BoCcos Back intersys GPS 
% 2.4 - BoCcos Back intersys GLO FDMA 
% 2.5 - BoCcos Back intersys GPS+GLO FDMA
% 3 - BPSK intersys
% 3.1 - BPSK intersys GLONASS FDMA
% 3.2 - BPSK intersys GPS+GLONASS FDMA
% 3.3 - BPSK Back intersys GPS 
% 3.4 - BPSK Back intersys GLO FDMA 
% 3.5 - BPSK Back intersys GPS+GLO FDMA
% 4 - BoC sin intrasys
% 5 - BoC cos intrasys
% 6 - BPSK cos intrasys
% 7 - BoC sin by FDMA
% 8 - BoC cos by FDMA
% 9 - BPSK cos by FDMA


BoCsincos_L1; % Список сигналов, удовлетворяющих РА

load('results/InterSysJam_BoCsin_L1_BoC_1_1.mat'); 
load('results/InterSysJam_BoCcos_L1_BoC_1_1.mat');
load('results/InterSysJam_BPSK_L1_BoC_1_1.mat');
load('results/InterSysJam_BoCsin_L1_BoC_6_1.mat'); 
load('results/InterSysJam_BoCcos_L1_BoC_6_1.mat');
load('results/InterSysJam_BPSK_L1_BoC_6_1.mat');
load('results/InterSysJam_BoCsin_L1_BoC_10_5.mat');
load('results/InterSysJam_BoCcos_L1_BoC_10_5.mat');
load('results/InterSysJam_BPSK_L1_BoC_10_5.mat');
load('results/InterSysJam_BoCsin_L1_BoC_0_1.mat');
load('results/InterSysJam_BoCcos_L1_BoC_0_1.mat');
load('results/InterSysJam_BPSK_L1_BoC_0_1.mat');
load('results/InterSysJam_BoCsin_L1_BoC_0_10.mat');
load('results/InterSysJam_BoCcos_L1_BoC_0_10.mat');
load('results/InterSysJam_BPSK_L1_BoC_0_10.mat');
load('InSysJam_BoCsin.mat');
load('InSysJam_BoCcos.mat');
load('InSysJam_BPSK.mat');
load('results/InterSysJam_BoCsin_L1_GloST.mat', 'InterSysJam_BoCsin_L1_GloST');
load('results/InterSysJam_BoCcos_L1_GloST.mat', 'InterSysJam_BoCcos_L1_GloST');
load('results/InterSysJam_BoCsin_L1_GloVT.mat', 'InterSysJam_BoCsin_L1_GloVT');
load('results/InterSysJam_BoCcos_L1_GloVT.mat', 'InterSysJam_BoCcos_L1_GloVT');
load('results/InterSysJam_BPSK_L1_GloST.mat', 'InterSysJam_BPSK_L1_GloST');
load('results/InterSysJam_BPSK_L1_GloVT.mat', 'InterSysJam_BPSK_L1_GloVT');
load('results/InterSysJam_BoCsin_L1_GloST_mean.mat', 'InterSysJam_BoCsin_L1_GloST_mean');
load('results/InterSysJam_BoCcos_L1_GloST_mean.mat', 'InterSysJam_BoCcos_L1_GloST_mean');
load('results/InterSysJam_BoCsin_L1_GloVT_mean.mat', 'InterSysJam_BoCsin_L1_GloVT_mean');
load('results/InterSysJam_BoCcos_L1_GloVT_mean.mat', 'InterSysJam_BoCcos_L1_GloVT_mean');
load('results/InterSysJam_BPSK_L1_GloST_mean.mat', 'InterSysJam_BPSK_L1_GloST_mean');
load('results/InterSysJam_BPSK_L1_GloVT_mean.mat', 'InterSysJam_BPSK_L1_GloVT_mean');

load('results/back/BackInterSysJam_BoCsin_L1_BoC_1_1.mat'); 
load('results/back/BackInterSysJam_BoCcos_L1_BoC_1_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_1_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_6_1.mat'); 
load('results/back/BackInterSysJam_BoCcos_L1_BoC_6_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_6_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_10_5.mat');
load('results/back/BackInterSysJam_BoCcos_L1_BoC_10_5.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_10_5.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_0_1.mat');
load('results/back/BackInterSysJam_BoCcos_L1_BoC_0_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_0_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_0_10.mat');
load('results/back/BackInterSysJam_BoCcos_L1_BoC_0_10.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_0_10.mat');
% load('results/back/BackInterSysJam_BoCsin_L1_GloST.mat', 'InterSysJam_BoCsin_L1_GloST');
% load('results/back/BackInterSysJam_BoCcos_L1_GloST.mat', 'InterSysJam_BoCcos_L1_GloST');
% load('results/back/BackInterSysJam_BoCsin_L1_GloVT.mat', 'InterSysJam_BoCsin_L1_GloVT');
% load('results/back/BackInterSysJam_BoCcos_L1_GloVT.mat', 'InterSysJam_BoCcos_L1_GloVT');
% load('results/back/BackInterSysJam_BPSK_L1_GloST.mat', 'InterSysJam_BPSK_L1_GloST');
% load('results/back/BackInterSysJam_BPSK_L1_GloVT.mat', 'InterSysJam_BPSK_L1_GloVT');
% load('results/back/BackInterSysJam_BoCsin_L1_GloST_mean.mat', 'InterSysJam_BoCsin_L1_GloST_mean');
% load('results/back/BackInterSysJam_BoCcos_L1_GloST_mean.mat', 'InterSysJam_BoCcos_L1_GloST_mean');
% load('results/back/BackInterSysJam_BoCsin_L1_GloVT_mean.mat', 'InterSysJam_BoCsin_L1_GloVT_mean');
% load('results/back/BackInterSysJam_BoCcos_L1_GloVT_mean.mat', 'InterSysJam_BoCcos_L1_GloVT_mean');
% load('results/back/BackInterSysJam_BPSK_L1_GloST_mean.mat', 'InterSysJam_BPSK_L1_GloST_mean');
% load('results/back/BackInterSysJam_BPSK_L1_GloVT_mean.mat', 'InterSysJam_BPSK_L1_GloVT_mean');

k_JN0_GPS_L1_P = -3;
k_JN0_GPS_L1_M = +0.5;
k_JN0_GPS_L1_L1C = +1.5;
k_JN0_GLO_L1_L1OF = -2.5;
k_JN0_GLO_L1_L1SF = -2.5;

Ml_GPS = 10;
Ml_GLO = 8;

if Table_num == 1
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCsin L1: GPS to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB\n', Ml_GPS, 0);
    fprintf('!GPS L1 P(Y)<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_P);
    fprintf('!GPS L1C MBOC<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_L1C);
    fprintf('!GPS L1 M<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_M);
    fprintf('!k<sub>intersys, GPS</sub> ');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)

            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BoCsin_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 1.1
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCsin L1: GLONASS FDMA to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GLONASS L1OF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1OF);    
    fprintf('!GLONASS L1SF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1SF);        
    fprintf('!k<sub>intersys, GLONASS</sub> ');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)

            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end            
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 1.2
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCsin L1: GPS and GLONASS FDMA to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB\n', Ml_GPS, 0);
    fprintf('!GPS L1 P(Y)<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_P);
    fprintf('!GPS L1C MBOC<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_L1C);
    fprintf('!GPS L1 M<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_M);
    fprintf('!GLONASS L1OF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1OF);    
    fprintf('!GLONASS L1SF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1SF);     
    fprintf('!k<sub>intersys</sub> ');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)

            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((InterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BoCsin_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if isnan(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end                  
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 1.3
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCsin L1: GLONASS CDMA to GPS, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A \n');
    fprintf('!GPS L1 P(Y) \n');
    fprintf('!GPS L1C MBOC \n');
    fprintf('!GPS L1 M \n');
    fprintf('!Mean');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)

            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))/10));
            end
            if isnan(BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))/10);
                
            end
            if isnan(BackInterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((BackInterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((BackInterSysJam_BoCsin_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC;
            end
            if isnan(BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))/10);
            end
            sum_dB = sum_dB / 4;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 1.4
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCsin L1: GLONASS CDMA to GLONASS FDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GLONASS L1OF\n');    
    fprintf('!GLONASS L1SF\n');   
    fprintf('!Mean');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)

            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))/10);
            end            
            if isnan(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))/10);
            end                 
            sum_dB = sum_dB / 2;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 1.5
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCsin L1: GLONASS CDMA to GPS and GLONASS FDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A \n');
    fprintf('!GPS L1 P(Y) \n');
    fprintf('!GPS L1C MBOC \n');
    fprintf('!GPS L1 M \n');
    fprintf('!GLONASS L1OF\n');    
    fprintf('!GLONASS L1SF\n');       
    fprintf('!Mean');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)

            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, freq_index))/10));
            end
            if isnan(BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, freq_index))/10);
                
            end
            if isnan(BackInterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((BackInterSysJam_BoCsin_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((BackInterSysJam_BoCsin_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC;
            end
            if isnan(BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, freq_index))/10);
            end
            if isnan(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))/10);
            end            
            if isnan(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))/10);
            end              
            sum_dB = sum_dB / 6;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 2
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCcos L1: GPS to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB\n', Ml_GPS, 0);
    fprintf('!GPS L1 P(Y)<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_P);
    fprintf('!GPS L1C MBOC<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_L1C);
    fprintf('!GPS L1 M<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_M);
    fprintf('!k<sub>intersys, GPS</sub> ');    
        for i = 1:size(BoCcos_Freq_L1_num, 1)

            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BoCcos_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 2.1
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCcos L1: GLONASS FDMA to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GLONASS L1OF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1OF);    
    fprintf('!GLONASS L1SF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1SF);  
    fprintf('!k<sub>intersys, GLONASS</sub> ');    
        for i = 1:size(BoCcos_Freq_L1_num, 1)

            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end               
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 2.2
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCcos L1: GPS and GLONASS FDMA to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB\n', Ml_GPS, 0);
    fprintf('!GPS L1 P(Y)<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_P);
    fprintf('!GPS L1C MBOC<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_L1C);
    fprintf('!GPS L1 M<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_M);
    fprintf('!GLONASS L1OF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1OF);    
    fprintf('!GLONASS L1SF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1SF);      
    fprintf('!k<sub>intersys</sub> ');    
        for i = 1:size(BoCcos_Freq_L1_num, 1)

            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((InterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((InterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BoCcos_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);
            end
            if isnan(InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if isnan(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end               
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 2.3
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCcos L1: GLONASS CDMA to GPS, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A \n');
    fprintf('!GPS L1 P(Y) \n');
    fprintf('!GPS L1C MBOC \n');
    fprintf('!GPS L1 M \n');
    fprintf('!Mean');    
        for i = 1:size(BoCcos_Freq_L1_num, 1)

            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))/10));
            end
            if isnan(BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))/10);
                
            end
            if isnan(BackInterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((BackInterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((BackInterSysJam_BoCcos_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC;
            end
            if isnan(BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))/10);
            end
            sum_dB = sum_dB / 4;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 2.4
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCсos L1: GLONASS CDMA to GLONASS FDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GLONASS L1OF\n');    
    fprintf('!GLONASS L1SF\n');   
    fprintf('!Mean');    
        for i = 1:size(BoCcos_Freq_L1_num, 1)

            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))/10);
            end            
            if isnan(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))/10);
            end                 
            sum_dB = sum_dB / 2;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 2.5
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BoCcos L1: GLONASS CDMA to GPS and GLONASS FDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A \n');
    fprintf('!GPS L1 P(Y) \n');
    fprintf('!GPS L1C MBOC \n');
    fprintf('!GPS L1 M \n');
    fprintf('!GLONASS L1OF\n');    
    fprintf('!GLONASS L1SF\n');       
    fprintf('!Mean');    
        for i = 1:size(BoCcos_Freq_L1_num, 1)

            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            if isnan(BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, freq_index))/10));
            end
            if isnan(BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, freq_index))/10);
                
            end
            if isnan(BackInterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((BackInterSysJam_BoCcos_L1_BoC_1_1(m8, n8, freq_index))/10) + ...
                          1/11  * 10^((BackInterSysJam_BoCcos_L1_BoC_6_1(m8, n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC;
            end
            if isnan(BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, freq_index))/10);
            end
            if isnan(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))/10);
            end            
            if isnan(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))/10);
            end              
            sum_dB = sum_dB / 6;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 3
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BPSK L1: GPS to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB\n', Ml_GPS, 0);
    fprintf('!GPS L1 P(Y)<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_P);
    fprintf('!GPS L1C MBOC<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_L1C);
    fprintf('!GPS L1 M<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_M);
    fprintf('!k<sub>intersys,GPS</sub> ');    
        for i = 1:size(BPSK_Freq_L1_num, 1)

            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            if isnan(InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))
                fprintf('|| n/d ');
            else
                IS_TMBOC = 10/11 * 10^((InterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BPSK_L1_BoC_6_1(n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);                
            end
            if isnan(InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 3.1
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BPSK L1: GLONASS FDMA to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GLONASS L1OF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1OF);    
    fprintf('!GLONASS L1SF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1SF);  
    fprintf('!k<sub>intersys, GLONASS</sub> ');    
        for i = 1:size(BPSK_Freq_L1_num, 1)

            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            if isnan(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end               
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 3.2
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BPSK L1: GPS and GLONASS FDMA to GLONASS CDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB\n', Ml_GPS, 0);
    fprintf('!GPS L1 P(Y)<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_P);
    fprintf('!GPS L1C MBOC<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_L1C);
    fprintf('!GPS L1 M<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GPS, k_JN0_GPS_L1_M);
    fprintf('!GLONASS L1OF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1OF);    
    fprintf('!GLONASS L1SF<br>M<sub>l</sub> = %.0f<br>k<sub>J/N0,l</sub>=%.1f dB \n', Ml_GLO, k_JN0_GLO_L1_L1SF);  
    fprintf('!k<sub>intersys</sub> ');    
        for i = 1:size(BPSK_Freq_L1_num, 1)

            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            if isnan(InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((InterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))/10))...
                    *Ml_GPS;
            end
            if isnan(InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))/10)...
                    *Ml_GPS*10^(k_JN0_GPS_L1_P/10);
                
            end
            if isnan(InterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))
                fprintf('|| n/d ');
            else
                IS_TMBOC = 10/11 * 10^((InterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))/10) + ...
                          1/11  * 10^((InterSysJam_BPSK_L1_BoC_6_1(n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC...
                    *Ml_GPS*10^(k_JN0_GPS_L1_L1C/10);                
            end
            if isnan(InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))/10)...
                     *Ml_GPS*10^(k_JN0_GPS_L1_M/10);
            end
            if isnan(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1OF/10);
            end            
            if isnan(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))/10)...
                     *Ml_GLO*10^(k_JN0_GLO_L1_L1SF/10);
            end              
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 3.3
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BPSK L1: GLONASS CDMA to GPS, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A \n');
    fprintf('!GPS L1 P(Y) \n');
    fprintf('!GPS L1C MBOC \n');
    fprintf('!GPS L1 M \n');
    fprintf('!Mean');    
        for i = 1:size(BPSK_Freq_L1_num, 1)

            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = 0;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            if isnan(BackInterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BPSK_L1_BoC_0_1(n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((BackInterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))/10));
            end
            if isnan(BackInterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BPSK_L1_BoC_0_10(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))/10);
                
            end
            if isnan(BackInterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((BackInterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))/10) + ...
                          1/11  * 10^((BackInterSysJam_BPSK_L1_BoC_6_1(n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC;
            end
            if isnan(BackInterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BPSK_L1_BoC_10_5(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))/10);
            end
            sum_dB = sum_dB / 4;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 3.4
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BPSK L1: GLONASS CDMA to GLONASS FDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GLONASS L1OF\n');    
    fprintf('!GLONASS L1SF\n');  
    fprintf('!Mean');    
        for i = 1:size(BPSK_Freq_L1_num, 1)

            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = 0;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            if isnan(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))/10);
            end            
            if isnan(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))/10);
            end                 
            sum_dB = sum_dB / 2;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 3.5
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Intersystem jammer for BPSK L1: GLONASS CDMA to GPS and GLONASS FDMA, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!GPS L1 C/A \n');
    fprintf('!GPS L1 P(Y) \n');
    fprintf('!GPS L1C MBOC \n');
    fprintf('!GPS L1 M \n');
    fprintf('!GLONASS L1OF\n');    
    fprintf('!GLONASS L1SF\n');     
    fprintf('!Mean');    
        for i = 1:size(BPSK_Freq_L1_num, 1)

            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = 0;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            sum_dB = 0;


            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            if isnan(BackInterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BPSK_L1_BoC_0_1(n8, freq_index)*10)/10);
                sum_dB = sum_dB + (10^((BackInterSysJam_BPSK_L1_BoC_0_1(n8, freq_index))/10));
            end
            if isnan(BackInterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BPSK_L1_BoC_0_10(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L1_BoC_0_10(n8, freq_index))/10);
                
            end
            if isnan(BackInterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))
                fprintf('|| n/d ');
            else
               IS_TMBOC = 10/11 * 10^((BackInterSysJam_BPSK_L1_BoC_1_1(n8, freq_index))/10) + ...
                          1/11  * 10^((BackInterSysJam_BPSK_L1_BoC_6_1(n8, freq_index))/10);
                
                fprintf('|| %.1f ', round( 10*log10(IS_TMBOC) *10)/10);
                sum_dB = sum_dB + IS_TMBOC;
            end
            if isnan(BackInterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(BackInterSysJam_BPSK_L1_BoC_10_5(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((BackInterSysJam_BPSK_L1_BoC_10_5(n8, freq_index))/10);
            end
            if isnan(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))/10);
            end            
            if isnan(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index)*10)/10);
                sum_dB = sum_dB + 10^((InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))/10);
            end                      
            sum_dB = sum_dB / 6;
            if sum_dB == 0
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(10*log10(sum_dB)*10)/10);
            end
        end
    fprintf('\n|} \n');
end


if Table_num == 4
%     fprintf('|+ Intrasystem jammer for BoCsin L1, dB\n');

        min_v = 0; min_n = 0; min_m = 0;
        first_flag = 1;
        n_in_row = 3;
        arr = BoCsin_Freq_L1_num;
        arr2 = InSysJam_BoCsin;
        for i = 1:size(arr, 1)

            if ~mod(i-1, n_in_row)
                if first_flag == 1
                    first_flag = 0;
                else
                    fprintf('%s\n|- align="center"\n%s\n', s1, s2);
                    fprintf('|} \n');
                end
                s1 = ['{| class="wikitable" border="1"'  sprintf('\n|-\n') '| ''''Signal''''  '];
                s2 = '|<math>k_{cd}, dB</math> ';
            end

            n = arr(i, 2); n8 = n*8;
            m = arr(i, 1); m8 = m*8;
           
            s1 = [s1 sprintf('|| BoCsin(%.3f, %.3f) ', m, n)];                
            s2 = [s2 '|| '];
            
            if isnan(arr2(m8, n8))
                s2 = [s2 'n/d '];
            else
                s2 = [s2 sprintf('%.1f ', round(arr2(m8, n8)*10)/10)];
                if arr2(m8, n8) < min_v
                    min_v = arr2(m8, n8);
                    min_n = n;
                    min_m = m;
                end
            end
            
            
            if i == size(arr, 1)
                fprintf('%s\n|- align="center"\n%s\n', s1, s2);
                fprintf('|} \n'); 
            end
        end
end


if Table_num == 5
%     fprintf('|+ Intrasystem jammer for BoCcos L1, dB\n');

        min_v = 0; min_n = 0; min_m = 0;
        first_flag = 1;
        n_in_row = 3;
        arr = BoCcos_Freq_L1_num;
        arr2 = InSysJam_BoCcos;
        for i = 1:size(arr, 1)

            if ~mod(i-1, n_in_row)
                if first_flag == 1
                    first_flag = 0;
                else
                    fprintf('%s\n|- align="center"\n%s\n', s1, s2);
                    fprintf('|} \n');
                end
                s1 = ['{| class="wikitable" border="1"'  sprintf('\n|-\n') '| ''''Signal''''  '];
                s2 = '|<math>k_{cd}, dB</math> ';
            end

            n = arr(i, 2); n8 = n*8;
            m = arr(i, 1); m8 = m*8;

            s1 = [s1 sprintf('|| BoCcos(%.3f, %.3f) ', m, n)];                
            s2 = [s2 '|| '];
            
            if isnan(arr2(m8, n8))
                s2 = [s2 'n/d '];
            else
                s2 = [s2 sprintf('%.1f ', round(arr2(m8, n8)*10)/10)];
                if arr2(m8, n8) < min_v
                    min_v = arr2(m8, n8);
                    min_n = n;
                    min_m = m;
                end
            end
            
            
            if i == size(arr, 1)
                fprintf('%s\n|- align="center"\n%s\n', s1, s2);
                fprintf('|} \n'); 
            end
        end
end


if Table_num == 6
%     fprintf('|+ Intrasystem jammer for BoCcos L1, dB\n');

        min_v = 0; min_n = 0; min_m = 0;
        first_flag = 1;
        n_in_row = 5;
        arr = BPSK_Freq_L1_num;
        arr2 = InSysJam_BPSK;
        for i = 1:size(arr, 1)

            if ~mod(i-1, n_in_row)
                if first_flag == 1
                    first_flag = 0;
                else
                    fprintf('%s\n|- align="center"\n%s\n', s1, s2);
                    fprintf('|} \n');
                end
                s1 = ['{| class="wikitable" border="1"'  sprintf('\n|-\n') '| ''''Signal''''  '];
                s2 = '|<math>k_{cd}, dB</math> ';
            end

            n = arr(i, 1); n8 = n*8;

            s1 = [s1 sprintf('|| BPSK(%.3f) ', n)];                
            s2 = [s2 '|| '];
            
            if isnan(arr2(n8))
                s2 = [s2 'n/d '];
            else
                s2 = [s2 sprintf('%.1f ', round(arr2(n8)*10)/10)];
                if arr2(n8) < min_v
                    min_v = arr2(n8);
                    min_n = n;
                end
            end
            
            
            if i == size(arr, 1)
                fprintf('%s\n|- align="center"\n%s\n', s1, s2);
                fprintf('|} \n'); 
            end
        end
end


if Table_num == 7
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Interference of BoC<sub>sin</sub>(m,n) and existing civilian GLONASS L1OF signals, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    
    lit_arr = -7:6;
    lit_size = length(lit_arr);
    
    for i = 1:lit_size
        fprintf('!l = %.0f\n', lit_arr(i));
    end
    fprintf('!Mean\n');   
    fprintf('!Max');    
        for i = 1:size(BoCsin_Freq_L1_num, 1)
            
            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            
            max_val = -1000;
            fprintf('\n|- align="center"\n');
            fprintf('|BoC<sub>sin</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            
            for lit = 1:lit_size    
                if isnan(InterSysJam_BoCsin_L1_GloST(m8, n8, freq_index, lit))
                    fprintf('|| n/d ');
                else
                    fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloST(m8, n8, freq_index, lit)*10)/10);
                    if max_val < round(InterSysJam_BoCsin_L1_GloST(m8, n8, freq_index, lit)*10)/10
                        max_val = round(InterSysJam_BoCsin_L1_GloST(m8, n8, freq_index, lit)*10)/10;
                    end
                end
            end
            
            if isnan(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloST_mean(m8, n8, freq_index)*10)/10);
            end
            fprintf('|| %.1f ', max_val);            
           
        end
    fprintf('\n|} \n \n \n');
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Interference of BoC<sub>sin</sub>(m,n) and existing military GLONASS L1SF signals, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    
    lit_arr = -7:6;
    lit_size = length(lit_arr);
    
    for i = 1:lit_size
        fprintf('!l = %.0f\n', lit_arr(i));
    end
    fprintf('!Mean\n');    
    fprintf('!Max');        
        for i = 1:size(BoCsin_Freq_L1_num, 1)
            
            n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCsin_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            
            max_val = -1000;
            fprintf('\n|- align="center"\n');
            fprintf('|BoC<sub>sin</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            
            for lit = 1:lit_size    
                if isnan(InterSysJam_BoCsin_L1_GloVT(m8, n8, freq_index, lit))
                    fprintf('|| n/d ');
                else
                    fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloVT(m8, n8, freq_index, lit)*10)/10);
                    if max_val < round(InterSysJam_BoCsin_L1_GloVT(m8, n8, freq_index, lit)*10)/10
                        max_val = round(InterSysJam_BoCsin_L1_GloVT(m8, n8, freq_index, lit)*10)/10;
                    end
                end
            end
            
            if isnan(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
            end
            fprintf('|| %.1f ', max_val);            
           
        end
    fprintf('\n|} \n');    
end


if Table_num == 8
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Interference of BoC<sub>cos</sub>(m,n) and existing civilian GLONASS L1OF signals, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    
    lit_arr = -7:6;
    lit_size = length(lit_arr);
    
    for i = 1:lit_size
        fprintf('!l = %.0f\n', lit_arr(i));
    end
    fprintf('!Mean\n');    
    fprintf('!Max');        
        for i = 1:size(BoCcos_Freq_L1_num, 1)
            
            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            
            max_val = -1000;
            fprintf('\n|- align="center"\n');
            fprintf('|BoC<sub>cos</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            
            for lit = 1:lit_size    
                if isnan(InterSysJam_BoCcos_L1_GloST(m8, n8, freq_index, lit))
                    fprintf('|| n/d ');
                else
                    fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloST(m8, n8, freq_index, lit)*10)/10);
                    if max_val < round(InterSysJam_BoCcos_L1_GloST(m8, n8, freq_index, lit)*10)/10
                        max_val = round(InterSysJam_BoCcos_L1_GloST(m8, n8, freq_index, lit)*10)/10;
                    end
                end
            end
            
            if isnan(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloST_mean(m8, n8, freq_index)*10)/10);
            end
            fprintf('|| %.1f ', max_val);            
           
        end
    fprintf('\n|} \n \n \n');
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Interference of BoC<sub>cos</sub>(m,n) and existing military GLONASS L1SF signals, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    
    lit_arr = -7:6;
    lit_size = length(lit_arr);
    
    for i = 1:lit_size
        fprintf('!l = %.0f\n', lit_arr(i));
    end
    fprintf('!Mean\n');    
    fprintf('!Max');        
        for i = 1:size(BoCcos_Freq_L1_num, 1)
            
            n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
            m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
            freq = BoCcos_Freq_L1_num(i, 3);
            freq_index = freq - 1558 + 1;
            
            max_val = -1000;
            fprintf('\n|- align="center"\n');
            fprintf('|BoC<sub>cos</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            
            for lit = 1:lit_size    
                if isnan(InterSysJam_BoCcos_L1_GloVT(m8, n8, freq_index, lit))
                    fprintf('|| n/d ');
                else
                    fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloVT(m8, n8, freq_index, lit)*10)/10);
                    if max_val < round(InterSysJam_BoCcos_L1_GloVT(m8, n8, freq_index, lit)*10)/10
                        max_val = round(InterSysJam_BoCcos_L1_GloVT(m8, n8, freq_index, lit)*10)/10;
                    end
                end
            end
            
            if isnan(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, freq_index)*10)/10);
            end
            fprintf('|| %.1f ', max_val);            
           
        end
    fprintf('\n|} \n');    
end


if Table_num == 9
    
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Interference of BPSK(n) and existing civilian GLONASS L1OF signals, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    
    lit_arr = -7:6;
    lit_size = length(lit_arr);
    
    for i = 1:lit_size
        fprintf('!l = %.0f\n', lit_arr(i));
    end
    fprintf('!Mean\n');   
    fprintf('!Max');        
        for i = 1:size(BPSK_Freq_L1_num, 1)
            
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            
            max_val = -1000;
            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            
            for lit = 1:lit_size    
                if isnan(InterSysJam_BPSK_L1_GloST(n8, freq_index, lit))
                    fprintf('|| n/d ');
                else
                    fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloST(n8, freq_index, lit)*10)/10);
                    if max_val < round(InterSysJam_BPSK_L1_GloST(n8, freq_index, lit)*10)/10
                        max_val = round(InterSysJam_BPSK_L1_GloST(n8, freq_index, lit)*10)/10;
                    end
                end
            end
            
            if isnan(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloST_mean(n8, freq_index)*10)/10);
            end
            
            fprintf('|| %.1f ', max_val);
           
        end
    fprintf('\n|} \n \n \n');
    
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Interference of BPSK(n) and existing military GLONASS L1SF signals, dB\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    
    lit_arr = -7:6;
    lit_size = length(lit_arr);
    
    for i = 1:lit_size
        fprintf('!l = %.0f\n', lit_arr(i));
    end
    fprintf('!Mean\n');    
    fprintf('!Max');    
    
        for i = 1:size(BPSK_Freq_L1_num, 1)
            
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;
            
            max_val = -1000;
            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            
            for lit = 1:lit_size    
                if isnan(InterSysJam_BPSK_L1_GloVT(n8, freq_index, lit))
                    fprintf('|| n/d ');
                else
                    fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloVT(n8, freq_index, lit)*10)/10);
                    if max_val < round(InterSysJam_BPSK_L1_GloVT(n8, freq_index, lit)*10)/10
                        max_val = round(InterSysJam_BPSK_L1_GloVT(n8, freq_index, lit)*10)/10;
                    end                                        
                end
            end
                        
            if isnan(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index))
                fprintf('|| n/d ');
            else
                fprintf('|| %.1f ', round(InterSysJam_BPSK_L1_GloVT_mean(n8, freq_index)*10)/10);
            end
            
            fprintf('|| %.1f ', max_val);
           
        end
    fprintf('\n|} \n');    
end