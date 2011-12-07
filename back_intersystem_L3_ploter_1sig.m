%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/back_intersystem_L3'];
path_to_pics = [pwd '/k_intersys/L3/BackGPS'];

load([path_to_results '/BackInterSysJam_BoCsin_L3_BoC_0_10.mat']);
load([path_to_results '/BackInterSysJam_BoCcos_L3_BoC_0_10.mat']);
load([path_to_results '/BackInterSysJam_BPSK_L3_BoC_0_10.mat']);

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1164:1184; fmax = length(farr); % Нормированный центральные частоты
Signal_Band = 20;
m8max = 80; n8max = 80;

k_JN0_GPS_L3_P = 0;
k_JN0_GPS_L3_L3C = 0;
k_JN0_GPS_L3_M = 0;

Out_Type = 2;
% 1 - С/A
% 2 - P(Y)
% 3 - MBOC
% 4 - M-code

DoPlot = 1;

if Out_Type == 1
    L3_mn = [1 0 0 0];
    filename_pref = 'CA_';
    signal_str = 'C';
    CLimm = [-135 -73];
elseif Out_Type == 2
    L3_mn = [0 1 0 0];
    filename_pref = 'PY_';
    signal_str = '';
    CLimm = [-113 -83];
elseif Out_Type == 3
%     L3_mn = [0 0 1 0]; 
%     filename_pref = 'MBOC_';
%     signal_str = 'MBOC';
%     CLimm = [-113 -77];
elseif Out_Type == 4
    L3_mn = [0 0 0 1];
    filename_pref = 'MCode_';
    signal_str = 'M-Code';
    CLimm = [-111 -71];
end
    
Ml_GPS = 1;
Ml_GLO = 0;

% CLimm = [-100 -65];

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
    
    L3C = 0*BackInterSysJam_BoCsin_L3_BoC_0_10;
    L3CA = 0*BackInterSysJam_BoCsin_L3_BoC_0_10;
    L3P = BackInterSysJam_BoCsin_L3_BoC_0_10;
    L3M = 0*BackInterSysJam_BoCsin_L3_BoC_0_10;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L3_mn(3)*10.^((L3C(:,:,f_in) + k_JN0_GPS_L3_L3C)/10) + ...
            L3_mn(1)*10.^((L3CA(:,:,f_in) + 0)/10) + ...
            L3_mn(4)*10.^((L3M(:,:,f_in) + k_JN0_GPS_L3_M)/10) + ...
            L3_mn(2)*10.^((L3P(:,:,f_in) + k_JN0_GPS_L3_P)/10) ...
            ) * Ml_GPS);
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
                title(sprintf('k_{sd} between BOC_{sin} and received GPS L5%s, f_n = %.0f, dB', signal_str, farr(f_in)), 'FontSize', 10);
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

    L3C = 0*BackInterSysJam_BoCcos_L3_BoC_0_10;
    L3CA = 0*BackInterSysJam_BoCcos_L3_BoC_0_10;
    L3P = BackInterSysJam_BoCcos_L3_BoC_0_10;
    L3M = 0*BackInterSysJam_BoCcos_L3_BoC_0_10;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L3_mn(3)*10.^((L3C(:,:,f_in) + k_JN0_GPS_L3_L3C)/10) + ...
            L3_mn(1)*10.^((L3CA(:,:,f_in) + 0)/10) + ...
            L3_mn(4)*10.^((L3M(:,:,f_in) + k_JN0_GPS_L3_M)/10) + ...
            L3_mn(2)*10.^((L3P(:,:,f_in) + k_JN0_GPS_L3_P)/10) ...
            ) * Ml_GPS);
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
                title(sprintf('k_{sd} between BOC_{cos} and received GPS L5%s, f_n = %.0f, dB', signal_str, farr(f_in)), 'FontSize', 10);
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
    
    L3C = 0*BackInterSysJam_BPSK_L3_BoC_0_10;
    L3CA = 0*BackInterSysJam_BPSK_L3_BoC_0_10;
    L3P = BackInterSysJam_BPSK_L3_BoC_0_10;
    L3M = 0*BackInterSysJam_BPSK_L3_BoC_0_10;
    
    hF = 0;
    sum_dB = 10*log10(...
        ( ...
        L3_mn(3)*10.^((L3C(:,:) + k_JN0_GPS_L3_L3C)/10) + ...
        L3_mn(1)*10.^((L3CA(:,:) + 0)/10) + ...
        L3_mn(4)*10.^((L3M(:,:) + k_JN0_GPS_L3_M)/10) + ...
        L3_mn(2)*10.^((L3P(:,:) + k_JN0_GPS_L3_P)/10) ...
        ) * Ml_GPS);    
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
        title(sprintf('k_{sd} between BPSK and received GPS L5%s, dB', signal_str), 'FontSize', 10);
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