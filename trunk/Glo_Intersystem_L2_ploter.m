%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/intersystem_L2'];
path_to_pics = [pwd '/k_intersys/L2/Glo'];


load([path_to_results '/InterSysJam_BoCsin_L2_GloST_mean.mat'], 'InterSysJam_BoCsin_L2_GloST_mean');
load([path_to_results '/InterSysJam_BoCcos_L2_GloST_mean.mat'], 'InterSysJam_BoCcos_L2_GloST_mean');
load([path_to_results '/InterSysJam_BPSK_L2_GloST_mean.mat'], 'InterSysJam_BPSK_L2_GloST_mean');

load([path_to_results '/InterSysJam_BoCsin_L2_GloVT_mean.mat'], 'InterSysJam_BoCsin_L2_GloVT_mean');
load([path_to_results '/InterSysJam_BoCcos_L2_GloVT_mean.mat'], 'InterSysJam_BoCcos_L2_GloVT_mean');
load([path_to_results '/InterSysJam_BPSK_L2_GloVT_mean.mat'], 'InterSysJam_BPSK_L2_GloVT_mean');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

k_JN0_GLO_L2_L2OF = 0;
k_JN0_GLO_L2_L2SF = 0;

Out_Type = 2;
% 1 - OF
% 2 - SF
% 3 - OF+SF

DoPlot = 1;

if Out_Type == 1
    L2OF_mn = 1;
    L2SF_mn = 0;
    filename_pref = 'OF_';
    signal_str = 'L2OF';
    CLimm = [-100 -65];
elseif Out_Type == 2
    L2OF_mn = 0;
    L2SF_mn = 1;
    filename_pref = 'SF_';
    signal_str = 'L2SF';
    CLimm = [-87 -69];
elseif Out_Type == 3
    L2OF_mn = 1;
    L2SF_mn = 1;
    filename_pref = 'Both_';
    signal_str = 'L2OF+L2SF';
end
    
Ml_GPS = 0;
Ml_GLO = 1;

for Signal_Type = 1:3
    sum_dB_min = 999; n8_min = -1; sum_dB_max = -999; n8_max = -1; 
    if Signal_Type == BOCsin

        L2OF = InterSysJam_BoCsin_L2_GloST_mean;
        L2SF = InterSysJam_BoCsin_L2_GloVT_mean;

        hF = 0;
        for f_in = 1:fmax
            sum_dB = 10*log10(...
                ( ...
                L2OF_mn*10.^((L2OF(:,:,f_in) + k_JN0_GLO_L2_L2OF)/10) + ...
                L2SF_mn*10.^((L2SF(:,:,f_in) + k_JN0_GLO_L2_L2SF)/10) ...
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
        fprintf('Maximum: %.1f for %s\n', round(10*sum_dB_max)/10, ['BoCsin(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
        fprintf('Minimum: %.1f for %s\n', round(10*sum_dB_min)/10, ['BoCsin(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    

    elseif Signal_Type == BOCcos

        L2OF = InterSysJam_BoCcos_L2_GloST_mean;
        L2SF = InterSysJam_BoCcos_L2_GloVT_mean;

        hF = 0;
        for f_in = 1:fmax
            sum_dB = 10*log10(...
                ( ...
                L2OF_mn*10.^((L2OF(:,:,f_in) + k_JN0_GLO_L2_L2OF)/10) + ...
                L2SF_mn*10.^((L2SF(:,:,f_in) + k_JN0_GLO_L2_L2SF)/10) ...
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
        fprintf('Maximum: %.1f for %s\n', round(10*sum_dB_max)/10, ['BoCcos(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
        fprintf('Minimum: %.1f for %s\n', round(10*sum_dB_min)/10, ['BoCcos(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    

    elseif Signal_Type == BPSK

        L2OF = InterSysJam_BPSK_L2_GloST_mean;
        L2SF = InterSysJam_BPSK_L2_GloVT_mean;

        hF = 0;
        sum_dB = 10*log10(...
            ( ...
            L2OF_mn*10.^((L2OF(:,:) + k_JN0_GLO_L2_L2OF)/10) + ...
            L2SF_mn*10.^((L2SF(:,:) + k_JN0_GLO_L2_L2SF)/10) ...
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
        fprintf('Maximun: %.1f for %s\n', round(10*sum_dB_max)/10, ['BPSK(' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
        fprintf('Minimum: %.1f for %s\n', round(10*sum_dB_min)/10, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
    end
end