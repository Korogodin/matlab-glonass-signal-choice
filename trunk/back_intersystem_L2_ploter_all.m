%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results_back = [pwd '/results/back_intersystem_L2'];
path_to_results = [pwd '/results/intersystem_L2'];
path_to_pics = [pwd '/k_intersys/L2/BackGPS'];

load([path_to_results_back '/BackInterSysJam_BoCsin_L2_BoC_10_5.mat']);
load([path_to_results_back '/BackInterSysJam_BoCcos_L2_BoC_10_5.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L2_BoC_10_5.mat']);
load([path_to_results_back '/BackInterSysJam_BoCsin_L2_BoC_0_1.mat']);
load([path_to_results_back '/BackInterSysJam_BoCcos_L2_BoC_0_1.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L2_BoC_0_1.mat']);
load([path_to_results_back '/BackInterSysJam_BoCsin_L2_BoC_0_10.mat']);
load([path_to_results_back '/BackInterSysJam_BoCcos_L2_BoC_0_10.mat']);
load([path_to_results_back '/BackInterSysJam_BPSK_L2_BoC_0_10.mat']);

load([path_to_results '/InterSysJam_BoCsin_L2_GloST_mean.mat'], 'InterSysJam_BoCsin_L2_GloST_mean');
load([path_to_results '/InterSysJam_BoCcos_L2_GloST_mean.mat'], 'InterSysJam_BoCcos_L2_GloST_mean');
load([path_to_results '/InterSysJam_BPSK_L2_GloST_mean.mat'], 'InterSysJam_BPSK_L2_GloST_mean');

load([path_to_results '/InterSysJam_BoCsin_L2_GloVT_mean.mat'], 'InterSysJam_BoCsin_L2_GloVT_mean');
load([path_to_results '/InterSysJam_BoCcos_L2_GloVT_mean.mat'], 'InterSysJam_BoCcos_L2_GloVT_mean');
load([path_to_results '/InterSysJam_BPSK_L2_GloVT_mean.mat'], 'InterSysJam_BPSK_L2_GloVT_mean');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
Signal_Band = 18;
m8max = 80; n8max = 80;

DoPlot = 1;

L2_mn = [1 1 0 1 1 1];
filename_pref = 'All_';
CLimm = [-90 -73];
    
Ml_GPS = 1/5;
Ml_GLO = 1/5;

k_JN0_GPS_L2_P = 0;
k_JN0_GPS_L2_M = 0;
k_JN0_GPS_L2_L2C = 0;
k_JN0_GLO_L2_L2OF = 0;
k_JN0_GLO_L2_L2SF = 0;

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
    
    L2C = 0*BackInterSysJam_BoCsin_L2_BoC_0_1;
    L2CA = BackInterSysJam_BoCsin_L2_BoC_0_1;
    L2P = BackInterSysJam_BoCsin_L2_BoC_0_10;
    L2M = BackInterSysJam_BoCsin_L2_BoC_10_5;
    L2OF = InterSysJam_BoCsin_L2_GloST_mean;
    L2SF = InterSysJam_BoCsin_L2_GloVT_mean;    
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L2_mn(3)*10.^((L2C(:,:,f_in) + k_JN0_GPS_L2_L2C)/10) + ...
            L2_mn(1)*10.^((L2CA(:,:,f_in) + 0)/10) + ...
            L2_mn(4)*10.^((L2M(:,:,f_in) + k_JN0_GPS_L2_M)/10) + ...
            L2_mn(2)*10.^((L2P(:,:,f_in) + k_JN0_GPS_L2_P)/10) ...
            ) * Ml_GPS + ...
            ( ...
            L2_mn(5)*10.^((L2OF(:,:,f_in) + k_JN0_GLO_L2_L2OF)/10) + ...
            L2_mn(6)*10.^((L2SF(:,:,f_in) + k_JN0_GLO_L2_L2SF)/10) ...
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
                title(sprintf('Mean k_{sd}: BOC_{sin} -> GPS L2 and GLONASS FDMA L2, f_n = %.0f, dB', farr(f_in)), 'FontSize', 10);
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

    L2C = 0*BackInterSysJam_BoCcos_L2_BoC_0_1;
    L2CA = BackInterSysJam_BoCcos_L2_BoC_0_1;
    L2P = BackInterSysJam_BoCcos_L2_BoC_0_10;
    L2M = BackInterSysJam_BoCcos_L2_BoC_10_5;
    L2OF = InterSysJam_BoCcos_L2_GloST_mean;
    L2SF = InterSysJam_BoCcos_L2_GloVT_mean;    
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L2_mn(3)*10.^((L2C(:,:,f_in) + k_JN0_GPS_L2_L2C)/10) + ...
            L2_mn(1)*10.^((L2CA(:,:,f_in) + 0)/10) + ...
            L2_mn(4)*10.^((L2M(:,:,f_in) + k_JN0_GPS_L2_M)/10) + ...
            L2_mn(2)*10.^((L2P(:,:,f_in) + k_JN0_GPS_L2_P)/10) ...
            ) * Ml_GPS + ...
            ( ...
            L2_mn(5)*10.^((L2OF(:,:,f_in) + k_JN0_GLO_L2_L2OF)/10) + ...
            L2_mn(6)*10.^((L2SF(:,:,f_in) + k_JN0_GLO_L2_L2SF)/10) ...
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
                title(sprintf('Mean k_{sd}: BOC_{cos} -> GPS L2 and GLONASS FDMA L2, f_n = %.0f, dB', farr(f_in)), 'FontSize', 10);
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
    
    L2C = 0*BackInterSysJam_BPSK_L2_BoC_0_1;
    L2CA = BackInterSysJam_BPSK_L2_BoC_0_1;
    L2P = BackInterSysJam_BPSK_L2_BoC_0_10;
    L2M = BackInterSysJam_BPSK_L2_BoC_10_5;
    L2OF = InterSysJam_BPSK_L2_GloST_mean;
    L2SF = InterSysJam_BPSK_L2_GloVT_mean;    
    
    hF = 0;
    sum_dB = 10*log10(...
        ( ...
        L2_mn(3)*10.^((L2C(:,:) + k_JN0_GPS_L2_L2C)/10) + ...
        L2_mn(1)*10.^((L2CA(:,:) + 0)/10) + ...
        L2_mn(4)*10.^((L2M(:,:) + k_JN0_GPS_L2_M)/10) + ...
        L2_mn(2)*10.^((L2P(:,:) + k_JN0_GPS_L2_P)/10) ...
        ) * Ml_GPS + ...
        ( ...
        L2_mn(5)*10.^((L2OF(:,:) + k_JN0_GLO_L2_L2OF)/10) + ...
        L2_mn(6)*10.^((L2SF(:,:) + k_JN0_GLO_L2_L2SF)/10) ...
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
        title(sprintf('Mean k_{sd}: BPSK -> GPS L2 and GLONASS FDMA L2, dB'), 'FontSize', 10);
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