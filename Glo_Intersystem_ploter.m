%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/intersystem_L1'];
path_to_pics = [pwd '/k_intersys/L1/Glo'];


load([path_to_results '/InterSysJam_BoCsin_L1_GloST_mean.mat'], 'InterSysJam_BoCsin_L1_GloST_mean');
load([path_to_results '/InterSysJam_BoCcos_L1_GloST_mean.mat'], 'InterSysJam_BoCcos_L1_GloST_mean');
load([path_to_results '/InterSysJam_BPSK_L1_GloST_mean.mat'], 'InterSysJam_BPSK_L1_GloST_mean');

load([path_to_results '/InterSysJam_BoCsin_L1_GloVT_mean.mat'], 'InterSysJam_BoCsin_L1_GloVT_mean');
load([path_to_results '/InterSysJam_BoCcos_L1_GloVT_mean.mat'], 'InterSysJam_BoCcos_L1_GloVT_mean');
load([path_to_results '/InterSysJam_BPSK_L1_GloVT_mean.mat'], 'InterSysJam_BPSK_L1_GloVT_mean');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

k_JN0_GLO_L1_L1OF = 0;
k_JN0_GLO_L1_L1SF = 0;

Out_Type = 2;
% 1 - OF
% 2 - SF
% 3 - OF+SF

DoPlot = 1;

if Out_Type == 1
    L1OF_mn = 1;
    L1SF_mn = 0;
    filename_pref = 'OF_';
    signal_str = ['L1OF'];
elseif Out_Type == 2
    L1OF_mn = 0;
    L1SF_mn = 1;
    filename_pref = 'SF_';
    signal_str = 'L1SF';
elseif Out_Type == 3
    L1OF_mn = 1;
    L1SF_mn = 1;
    filename_pref = 'Both_';
    signal_str = 'L1OF+L1SF';
end
    
Ml_GPS = 0;
Ml_GLO = 1;

% CLimm = [-100 -65];
CLimm = [-87 -69];
sum_dB_min = 999; n8_min = -1; sum_dB_max = -999; n8_max = -1; 
if Signal_Type == BOCsin
    
    L1OF = InterSysJam_BoCsin_L1_GloST_mean;
    L1SF = InterSysJam_BoCsin_L1_GloVT_mean;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L1OF_mn*10.^((L1OF(:,:,f_in) + k_JN0_GLO_L1_L1OF)/10) + ...
            L1SF_mn*10.^((L1SF(:,:,f_in) + k_JN0_GLO_L1_L1SF)/10) ...
            ) * Ml_GLO ...
            );
            if DoPlot
                hF = figure(hF+1);
                pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
                xlabel('n', 'FontSize', 10)
                ylabel('m', 'FontSize', 10)      
                colorbar
                set(gca, 'CLim', CLimm)
                set(gca, 'FontSize', 10)
                set(hF, 'PaperPosition', [0 0 14 11]);    
                title(sprintf('Mean k_{sd} between BOC_{sin} and GLONASS %s, f_n = %.0f, dB', signal_str, farr(f_in)), 'FontSize', 10);
                saveas(hF, [path_to_pics '/png/' filename_pref 'k_intersys_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.png']);
                saveas(hF, [path_to_pics '/fig/' filename_pref 'k_intersys_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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

    L1OF = InterSysJam_BoCcos_L1_GloST_mean;
    L1SF = InterSysJam_BoCcos_L1_GloVT_mean;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            L1OF_mn*10.^((L1OF(:,:,f_in) + k_JN0_GLO_L1_L1OF)/10) + ...
            L1SF_mn*10.^((L1SF(:,:,f_in) + k_JN0_GLO_L1_L1SF)/10) ...
            ) * Ml_GLO ...
            );
            if DoPlot
                hF = figure(hF+1);
                pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
                xlabel('n', 'FontSize', 10)
                ylabel('m', 'FontSize', 10)      
                colorbar
                set(gca, 'CLim', CLimm)
                set(gca, 'FontSize', 10)
                set(hF, 'PaperPosition', [0 0 14 11]);    
                title(sprintf('Mean k_{sd} between BOC_{cos} and GLONASS %s, f_n = %.0f, dB', signal_str, farr(f_in)), 'FontSize', 10);
                saveas(hF, [path_to_pics '/png/' filename_pref 'k_intersys_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.png']);
                saveas(hF, [path_to_pics '/fig/' filename_pref 'k_intersys_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    
    L1OF = InterSysJam_BPSK_L1_GloST_mean;
    L1SF = InterSysJam_BPSK_L1_GloVT_mean;
    
    hF = 0;
    sum_dB = 10*log10(...
        ( ...
        L1OF_mn*10.^((L1OF(:,:) + k_JN0_GLO_L1_L1OF)/10) + ...
        L1SF_mn*10.^((L1SF(:,:) + k_JN0_GLO_L1_L1SF)/10) ...
        ) * Ml_GLO ...
        );
    if DoPlot
        hF = figure(hF+1);
        pcolor((1:80)/8, [farr farr(fmax)+1], [sum_dB(1:80,1:fmax), nan(80,1)]');
        xlabel('n', 'FontSize', 10)
        ylabel('f_{n}', 'FontSize', 10)      
        colorbar
        set(gca, 'CLim', CLimm)
        set(gca, 'FontSize', 10)
        set(hF, 'PaperPosition', [0 0 14 11]);             
        title(sprintf('Mean k_{sd} between BPSK and GLONASS %s, dB', signal_str), 'FontSize', 10);
        saveas(hF, [path_to_pics '/png/' filename_pref 'k_intersys_BPSK.png']);
        saveas(hF, [path_to_pics '/fig/' filename_pref 'k_intersys_BPSK.fig']);
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
    fprintf('Maximun: %.3f for %s\n', sum_dB_max, ['BPSK(' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
end