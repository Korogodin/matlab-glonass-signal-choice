%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/intersystem_L3'];
path_to_pics = [pwd '/k_intersys/L3'];

load([path_to_results '/InterSysJam_BoCsin_L3_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BoCcos_L3_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BPSK_L3_BoC_0_10.mat']);

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

DoPlot = 1;

farr = 1164:1184; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

k_JN0_GPS_L3_P = 0;

Ml_GPS = 10;

for Signal_Type = 1:3
sum_dB_min = 999; sum_dB_max = -999; n8_min = -1; n8_max = -1;
if Signal_Type == BOCsin
    
    L3P = InterSysJam_BoCsin_L3_BoC_0_10;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            10.^((L3P(:,:,f_in) + k_JN0_GPS_L3_P)/10) ...
            ) * Ml_GPS ...
            );
        if DoPlot
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 10)
            ylabel('m', 'FontSize', 10)      
            colorbar
            set(gca, 'CLim', [-103 -73])
            set(gca, 'FontSize', 10)
            set(hF, 'PaperPosition', [0 0 14 11]);          
            title(sprintf('k_{intersys} for BOC_{sin} L3, f_n = %.0f', farr(f_in)), 'FontSize', 10);
            saveas(hF, [path_to_pics '/png/k_intersys_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/k_intersys_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    fprintf('Maximum: %.3f for %s\n', sum_dB_max, ['BoCsin(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BoCsin(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
elseif Signal_Type == BOCcos
    L3P = InterSysJam_BoCcos_L3_BoC_0_10;
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            10.^((L3P(:,:,f_in) + k_JN0_GPS_L3_P)/10) ...
            ) * Ml_GPS ...
            );

        if DoPlot
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 10)
            ylabel('m', 'FontSize', 10)      
            colorbar
            set(gca, 'CLim', [-103 -73])
            set(gca, 'FontSize', 10)
            set(hF, 'PaperPosition', [0 0 14 11]);         
            title(sprintf('k_{intersys} for BOC_{cos} L3, f_n = %.0f', farr(f_in)), 'FontSize', 10);
            saveas(hF, [path_to_pics '/png/k_intersys_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/k_intersys_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    fprintf('Maximum: %.3f for %s\n', sum_dB_max, ['BoCcos(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );        
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BoCcos(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
elseif Signal_Type == BPSK
    
    L3P = InterSysJam_BPSK_L3_BoC_0_10;
    
    hF = 0;
        sum_dB = 10*log10(...
            ( ...
            10.^((L3P(:,:) + k_JN0_GPS_L3_P)/10) ...
            ) * Ml_GPS ...
            );

    if DoPlot
        hF = figure(hF+1);
        pcolor((1:80)/8, [farr farr(fmax)+1], [sum_dB(1:80,1:fmax), nan(80,1)]');
        xlabel('n', 'FontSize', 10)
        ylabel('f_{n}', 'FontSize', 10)      
        colorbar
        set(gca, 'CLim', [-103 -73])
        set(gca, 'FontSize', 10)
        set(hF, 'PaperPosition', [0 0 14 11]);              
        title('k_{intersys} for BPSK(n)', 'FontSize', 10);
        saveas(hF, [path_to_pics '/png/k_intersys_BPSK.png']);
        saveas(hF, [path_to_pics '/fig/k_intersys_BPSK.fig']);
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
    fprintf('Maximum: %.3f for %s\n', sum_dB_max, ['BPSK(' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
end
end