%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results_back = [pwd '/results/back_intersystem_L1'];
path_to_results = [pwd '/results/intersystem_L1'];
path_to_pics = [pwd '/k_intersys/L1/BackGPS'];

load([path_to_results_back '/BackInterSysJam_BoCsin_L1_BoC_1_1.mat']); 
load([path_to_results_back '/BackInterSysJam_BoCcos_L1_BoC_1_1.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L1_BoC_1_1.mat']);
load([path_to_results_back '/BackInterSysJam_BoCsin_L1_BoC_6_1.mat']); 
load([path_to_results_back '/BackInterSysJam_BoCcos_L1_BoC_6_1.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L1_BoC_6_1.mat']);
load([path_to_results_back '/BackInterSysJam_BoCsin_L1_BoC_10_5.mat']);
load([path_to_results_back '/BackInterSysJam_BoCcos_L1_BoC_10_5.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L1_BoC_10_5.mat']);
load([path_to_results_back '/BackInterSysJam_BoCsin_L1_BoC_0_1.mat']);
load([path_to_results_back '/BackInterSysJam_BoCcos_L1_BoC_0_1.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L1_BoC_0_1.mat']);
load([path_to_results_back '/BackInterSysJam_BoCsin_L1_BoC_0_10.mat']);
load([path_to_results_back '/BackInterSysJam_BoCcos_L1_BoC_0_10.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L1_BoC_0_10.mat']);

load([path_to_results '/InterSysJam_BoCsin_L1_GloST_mean.mat'], 'InterSysJam_BoCsin_L1_GloST_mean');
load([path_to_results '/InterSysJam_BoCcos_L1_GloST_mean.mat'], 'InterSysJam_BoCcos_L1_GloST_mean');
load([path_to_results '/InterSysJam_BPSK_L1_GloST_mean.mat'], 'InterSysJam_BPSK_L1_GloST_mean');

load([path_to_results '/InterSysJam_BoCsin_L1_GloVT_mean.mat'], 'InterSysJam_BoCsin_L1_GloVT_mean');
load([path_to_results '/InterSysJam_BoCcos_L1_GloVT_mean.mat'], 'InterSysJam_BoCcos_L1_GloVT_mean');
load([path_to_results '/InterSysJam_BPSK_L1_GloVT_mean.mat'], 'InterSysJam_BPSK_L1_GloVT_mean');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
Signal_Band = 16;
m8max = 80; n8max = 80;

DoPlot = 1;

L1_mn = [1 1 1 1 1 1];
filename_pref = 'All_';
CLimm = [-91 -73];
    
Ml_GPS = 1/6;
Ml_GLO = 1/6;

k_JN0_GPS_L1_P = 0;
k_JN0_GPS_L1_M = 0;
k_JN0_GPS_L1_L1C = 0;
k_JN0_GLO_L1_L1OF = 0;
k_JN0_GLO_L1_L1SF = 0;

NaNder = nan(80, 80, fmax);
for freq_index = 1:fmax
    for n8 = 1:80
        for m8 = 1:80
            if m8 < n8
                continue;
            end
            if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (Signal_Band-freq_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end        
            NaNder(m8, n8, freq_index) = 0;
        end
    end
end

BPSK_NaN = nan(80, fmax);
for freq_index = 1:fmax
    for n8 = 1:80
        m8 = 0;
        if ((m8+n8)/8 > freq_index + 2) || ((m8+n8)/8 > (Signal_Band-freq_index) +2) % Если этот сигнал не влазиет в полосу
            continue;
        end        
        BPSK_NaN(n8, freq_index) = 0;
    end
end

for Signal_Type = 1:3
sum_dB_min = 999; n8_min = -1; sum_dB_max = -999; n8_max = -1; 
if Signal_Type == BOCsin
    
    L1C = 10*log10(10/11 * 10.^(BackInterSysJam_BoCsin_L1_BoC_1_1/10) + 1/11 * 10.^(BackInterSysJam_BoCsin_L1_BoC_6_1/10) ) ;
    L1CA = BackInterSysJam_BoCsin_L1_BoC_0_1;
    L1P = BackInterSysJam_BoCsin_L1_BoC_0_10;
    L1M = BackInterSysJam_BoCsin_L1_BoC_10_5;
    L1OF = InterSysJam_BoCsin_L1_GloST_mean;
    L1SF = InterSysJam_BoCsin_L1_GloVT_mean;    
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L1_mn(3)*10.^((L1C(:,:,f_in) + k_JN0_GPS_L1_L1C)/10) + ...
            L1_mn(1)*10.^((L1CA(:,:,f_in) + 0)/10) + ...
            L1_mn(4)*10.^((L1M(:,:,f_in) + k_JN0_GPS_L1_M)/10) + ...
            L1_mn(2)*10.^((L1P(:,:,f_in) + k_JN0_GPS_L1_P)/10) ...
            ) * Ml_GPS + ...
            ( ...
            L1_mn(5)*10.^((L1OF(:,:,f_in) + k_JN0_GLO_L1_L1OF)/10) + ...
            L1_mn(6)*10.^((L1SF(:,:,f_in) + k_JN0_GLO_L1_L1SF)/10) ...
            ) * Ml_GLO ...
            );        
        sum_dB = sum_dB + NaNder(:, :, f_in);
            if DoPlot
                hF = figure(hF+1);
                pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
                xlabel('n', 'FontSize', 10)
                ylabel('m', 'FontSize', 10)      
                colorbar
                set(gca, 'CLim', CLimm)
                set(gca, 'FontSize', 10)
                set(hF, 'PaperPosition', [0 0 14 11]);    
                title(sprintf('Mean k_{sd}: BOC_{sin} -> GPS L1 and GLONASS FDMA L1, f_n = %.0f, dB', farr(f_in)), 'FontSize', 10);
                saveas(hF, [path_to_pics '/png/' filename_pref 'k_sd_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.png']);
                saveas(hF, [path_to_pics '/fig/' filename_pref 'k_sd_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.fig']);
            end
            [a b] = min(sum_dB);
            [c d] = min(min(sum_dB));
            if sum_dB_min > c
                sum_dB_min = c;
                m8_min = b(d);
                n8_min = d;
                freq_min = f_in;
            end
            [a b] = max(sum_dB);
            [c d] = max(max(sum_dB));
            if sum_dB_max < c
                sum_dB_max = c;
                m8_max = b(d);
                n8_max = d;
                freq_max = f_in;
            end            
    end
    fprintf('Maximum: %.1f for %s\n', round(10*sum_dB_max)/10, ['BoCsin(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.1f for %s\n', round(10*sum_dB_min)/10, ['BoCsin(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    

elseif Signal_Type == BOCcos

    L1C = 10*log10(10/11 * 10.^(BackInterSysJam_BoCcos_L1_BoC_1_1/10) + 1/11 * 10.^(BackInterSysJam_BoCcos_L1_BoC_6_1/10) ) ;
    L1CA = BackInterSysJam_BoCcos_L1_BoC_0_1;
    L1P = BackInterSysJam_BoCcos_L1_BoC_0_10;
    L1M = BackInterSysJam_BoCcos_L1_BoC_10_5;
    L1OF = InterSysJam_BoCcos_L1_GloST_mean;
    L1SF = InterSysJam_BoCcos_L1_GloVT_mean;    
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L1_mn(3)*10.^((L1C(:,:,f_in) + k_JN0_GPS_L1_L1C)/10) + ...
            L1_mn(1)*10.^((L1CA(:,:,f_in) + 0)/10) + ...
            L1_mn(4)*10.^((L1M(:,:,f_in) + k_JN0_GPS_L1_M)/10) + ...
            L1_mn(2)*10.^((L1P(:,:,f_in) + k_JN0_GPS_L1_P)/10) ...
            ) * Ml_GPS + ...
            ( ...
            L1_mn(5)*10.^((L1OF(:,:,f_in) + k_JN0_GLO_L1_L1OF)/10) + ...
            L1_mn(6)*10.^((L1SF(:,:,f_in) + k_JN0_GLO_L1_L1SF)/10) ...
            ) * Ml_GLO ...
            );        
        sum_dB = sum_dB + NaNder(:, :, f_in);
            if DoPlot
                hF = figure(hF+1);
                pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
                xlabel('n', 'FontSize', 10)
                ylabel('m', 'FontSize', 10)      
                colorbar
                set(gca, 'CLim', CLimm)
                set(gca, 'FontSize', 10)
                set(hF, 'PaperPosition', [0 0 14 11]);    
                title(sprintf('Mean k_{sd}: BOC_{cos} -> GPS L1 and GLONASS FDMA L1, f_n = %.0f, dB', farr(f_in)), 'FontSize', 10);
                saveas(hF, [path_to_pics '/png/' filename_pref 'k_sd_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.png']);
                saveas(hF, [path_to_pics '/fig/' filename_pref 'k_sd_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.fig']);
            end
            [a b] = min(sum_dB);
            [c d] = min(min(sum_dB));
            if sum_dB_min > c
                sum_dB_min = c;
                m8_min = b(d);
                n8_min = d;
                freq_min = f_in;
            end
            [a b] = max(sum_dB);
            [c d] = max(max(sum_dB));
            if sum_dB_max < c
                sum_dB_max = c;
                m8_max = b(d);
                n8_max = d;
                freq_max = f_in;
            end            
    end
    fprintf('Maximum: %.1f for %s\n', round(10*sum_dB_max)/10, ['BoCcos(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.1f for %s\n', round(10*sum_dB_min)/10, ['BoCcos(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    

elseif Signal_Type == BPSK
    
    L1C = 10*log10(10/11 * 10.^(BackInterSysJam_BPSK_L1_BoC_1_1/10) + 1/11 * 10.^(BackInterSysJam_BPSK_L1_BoC_6_1/10) ) ;
    L1CA = BackInterSysJam_BPSK_L1_BoC_0_1;
    L1P = BackInterSysJam_BPSK_L1_BoC_0_10;
    L1M = BackInterSysJam_BPSK_L1_BoC_10_5;
    L1OF = InterSysJam_BPSK_L1_GloST_mean;
    L1SF = InterSysJam_BPSK_L1_GloVT_mean;    
    
    hF = 0;
    sum_dB = 10*log10(...
        ( ...
        L1_mn(3)*10.^((L1C(:,:) + k_JN0_GPS_L1_L1C)/10) + ...
        L1_mn(1)*10.^((L1CA(:,:) + 0)/10) + ...
        L1_mn(4)*10.^((L1M(:,:) + k_JN0_GPS_L1_M)/10) + ...
        L1_mn(2)*10.^((L1P(:,:) + k_JN0_GPS_L1_P)/10) ...
        ) * Ml_GPS + ...
        ( ...
        L1_mn(5)*10.^((L1OF(:,:) + k_JN0_GLO_L1_L1OF)/10) + ...
        L1_mn(6)*10.^((L1SF(:,:) + k_JN0_GLO_L1_L1SF)/10) ...
        ) * Ml_GLO ...
        );            
    sum_dB = sum_dB + BPSK_NaN;
    if DoPlot
        hF = figure(hF+1);
        pcolor((1:80)/8, [farr farr(fmax)+1], [sum_dB(1:80,1:fmax), nan(80,1)]');
        xlabel('n', 'FontSize', 10)
        ylabel('f_{n}', 'FontSize', 10)      
        colorbar
        set(gca, 'CLim', CLimm)
        set(gca, 'FontSize', 10)
        set(hF, 'PaperPosition', [0 0 14 11]);             
        title(sprintf('Mean k_{sd}: BPSK -> GPS L1 and GLONASS FDMA L1, dB'), 'FontSize', 10);
        saveas(hF, [path_to_pics '/png/' filename_pref 'k_sd_BPSK.png']);
        saveas(hF, [path_to_pics '/fig/' filename_pref 'k_sd_BPSK.fig']);
    end
    [a b] = min(sum_dB);
    [c d] = min(min(sum_dB));
    if sum_dB_min > c
        sum_dB_min = c;
        n8_min = b(d);
        freq_min = d;
    end
    [a b] = max(sum_dB);
    [c d] = max(max(sum_dB));
    if sum_dB_max < c
        sum_dB_max = c;
        n8_max = b(d);
        freq_max = d;
    end    
    fprintf('Maximum: %.1f for %s\n', round(10*sum_dB_max)/10, ['BPSK(' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.1f for %s\n', round(10*sum_dB_min)/10, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
end
end