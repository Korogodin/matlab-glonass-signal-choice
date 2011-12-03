%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/intersystem_L1'];
path_to_pics = [pwd '/k_intersys/L1'];

load([path_to_results '/InterSysJam_BoCsin_L1_BoC_1_1.mat']); 
load([path_to_results '/InterSysJam_BoCcos_L1_BoC_1_1.mat']);
load([path_to_results '/InterSysJam_BPSK_L1_BoC_1_1.mat']);
load([path_to_results '/InterSysJam_BoCsin_L1_BoC_6_1.mat']); 
load([path_to_results '/InterSysJam_BoCcos_L1_BoC_6_1.mat']);
load([path_to_results '/InterSysJam_BPSK_L1_BoC_6_1.mat']);
load([path_to_results '/InterSysJam_BoCsin_L1_BoC_10_5.mat']);
load([path_to_results '/InterSysJam_BoCcos_L1_BoC_10_5.mat']);
load([path_to_results '/InterSysJam_BPSK_L1_BoC_10_5.mat']);
load([path_to_results '/InterSysJam_BoCsin_L1_BoC_0_1.mat']);
load([path_to_results '/InterSysJam_BoCcos_L1_BoC_0_1.mat']);
load([path_to_results '/InterSysJam_BPSK_L1_BoC_0_1.mat']);
load([path_to_results '/InterSysJam_BoCsin_L1_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BoCcos_L1_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BPSK_L1_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BoCsin_L1_GloST_mean.mat'], 'InterSysJam_BoCsin_L1_GloST_mean');
load([path_to_results '/InterSysJam_BoCcos_L1_GloST_mean.mat'], 'InterSysJam_BoCcos_L1_GloST_mean');
load([path_to_results '/InterSysJam_BoCsin_L1_GloVT_mean.mat'], 'InterSysJam_BoCsin_L1_GloVT_mean');
load([path_to_results '/InterSysJam_BoCcos_L1_GloVT_mean.mat'], 'InterSysJam_BoCcos_L1_GloVT_mean');
load([path_to_results '/InterSysJam_BPSK_L1_GloST_mean.mat'], 'InterSysJam_BPSK_L1_GloST_mean');
load([path_to_results '/InterSysJam_BPSK_L1_GloVT_mean.mat'], 'InterSysJam_BPSK_L1_GloVT_mean');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 2; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

k_JN0_GPS_L1_P = -3;
k_JN0_GPS_L1_M = +0.5;
k_JN0_GPS_L1_L1C = +1.5;
k_JN0_GLO_L1_L1OF = -2.5;
k_JN0_GLO_L1_L1SF = -2.5;

Ml_GPS = 10;
Ml_GLO = 8;

sum_dB_min = 999; n8_min = -1; m8_min = -1; freq_min = -10;

if Signal_Type == BOCsin
    L1C = 10*log10(10/11 * 10.^(InterSysJam_BoCsin_L1_BoC_1_1/10) + 1/11 * 10.^(InterSysJam_BoCsin_L1_BoC_6_1/10) ) ;
    L1CA = InterSysJam_BoCsin_L1_BoC_0_1;
    L1P = InterSysJam_BoCsin_L1_BoC_0_10;
    L1M = InterSysJam_BoCsin_L1_BoC_10_5;
    L1OF = InterSysJam_BoCsin_L1_GloST_mean;
    L1SF = InterSysJam_BoCsin_L1_GloVT_mean;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            10.^((L1C(:,:,f_in) + k_JN0_GPS_L1_L1C)/10) + ...
            10.^((L1CA(:,:,f_in) + 0)/10) + ...
            10.^((L1M(:,:,f_in) + k_JN0_GPS_L1_M)/10) + ...
            10.^((L1P(:,:,f_in) + k_JN0_GPS_L1_P)/10) ...
            ) * Ml_GPS + ...
            ( ...
            10.^((L1OF(:,:,f_in) + k_JN0_GLO_L1_L1OF)/10) + ...
            10.^((L1SF(:,:,f_in) + k_JN0_GLO_L1_L1SF)/10) ...
            ) * Ml_GLO ...
            );
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 14)
            ylabel('m', 'FontSize', 14)      
            colorbar
            set(gca, 'CLim', [-80 -60])
            set(gca, 'FontSize', 14)
            set(hF, 'Position', [300 300 558 436])%496 388]);             
            title(sprintf('k_{intersys} for BOC_{sin}, f_n = %.0f', farr(f_in)), 'FontSize', 14);
            saveas(hF, [path_to_pics '/png/k_intersys_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/k_intersys_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.fig']);
            [a b] = min(sum_dB);
            [c d] = min(min(sum_dB));
            if sum_dB_min > c
                sum_dB_min = c;
                m8_min = b(d);
                n8_min = d;
                freq_min = f_in;
            end
    end
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BoCsin(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
elseif Signal_Type == BOCcos
    L1C = 10*log10(10/11 * 10.^(InterSysJam_BoCcos_L1_BoC_1_1/10) + 1/11 * 10.^(InterSysJam_BoCcos_L1_BoC_6_1/10) ) ;
    L1CA = InterSysJam_BoCcos_L1_BoC_0_1;
    L1P = InterSysJam_BoCcos_L1_BoC_0_10;
    L1M = InterSysJam_BoCcos_L1_BoC_10_5;
    L1OF = InterSysJam_BoCcos_L1_GloST_mean;
    L1SF = InterSysJam_BoCcos_L1_GloVT_mean;
    
    hF = 0;
    for f_in = 1:fmax
        sum_dB = 10*log10(...
            ( ...
            10.^((L1C(:,:,f_in) + k_JN0_GPS_L1_L1C)/10) + ...
            10.^((L1CA(:,:,f_in) + 0)/10) + ...
            10.^((L1M(:,:,f_in) + k_JN0_GPS_L1_M)/10) + ...
            10.^((L1P(:,:,f_in) + k_JN0_GPS_L1_P)/10) ...
            ) * Ml_GPS + ...
            ( ...
            10.^((L1OF(:,:,f_in) + k_JN0_GLO_L1_L1OF)/10) + ...
            10.^((L1SF(:,:,f_in) + k_JN0_GLO_L1_L1SF)/10) ...
            ) * Ml_GLO ...
            );
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 14)
            ylabel('m', 'FontSize', 14)      
            colorbar
            set(gca, 'CLim', [-80 -60])
            set(gca, 'FontSize', 14)
            set(hF, 'Position', [300 300 558 436])%496 388]);             
            title(sprintf('k_{intersys} for BOC_{cos}, f_n = %.0f', farr(f_in)), 'FontSize', 14);
            saveas(hF, [path_to_pics '/png/k_intersys_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/k_intersys_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.fig']);
            [a b] = min(sum_dB);
            [c d] = min(min(sum_dB));
            if sum_dB_min > c
                sum_dB_min = c;
                m8_min = b(d);
                n8_min = d;
                freq_min = f_in;
            end            
    end
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BoCcos(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
elseif Signal_Type == BPSK
    
    L1C = 10*log10(10/11 * 10.^(InterSysJam_BPSK_L1_BoC_1_1/10) + 1/11 * 10.^(InterSysJam_BPSK_L1_BoC_6_1/10) ) ;
    L1CA = InterSysJam_BPSK_L1_BoC_0_1;
    L1P = InterSysJam_BPSK_L1_BoC_0_10;
    L1M = InterSysJam_BPSK_L1_BoC_10_5;
    L1OF = InterSysJam_BPSK_L1_GloST_mean;
    L1SF = InterSysJam_BPSK_L1_GloVT_mean;
    
    hF = 0;
    sum_dB = 10*log10(...
        ( ...
        10.^((L1C(:,:) + k_JN0_GPS_L1_L1C)/10) + ...
        10.^((L1CA(:,:) + 0)/10) + ...
        10.^((L1M(:,:) + k_JN0_GPS_L1_M)/10) + ...
        10.^((L1P(:,:) + k_JN0_GPS_L1_P)/10) ...
        ) * Ml_GPS + ...
        ( ...
        10.^((L1OF(:,:) + k_JN0_GLO_L1_L1OF)/10) + ...
        10.^((L1SF(:,:) + k_JN0_GLO_L1_L1SF)/10) ...
        ) * Ml_GLO ...
        );
    hF = figure(hF+1);
    pcolor((1:80)/8, [farr farr(fmax)+1], [sum_dB(1:80,1:fmax), nan(80,1)]');
    xlabel('n', 'FontSize', 14)
    ylabel('f_{n}', 'FontSize', 14)      
    colorbar
    set(gca, 'CLim', [-80 -60])
    set(gca, 'FontSize', 14)
    set(hF, 'Position', [300 300 833 554])%496 388]);             
    title('k_{intersys} for BPSK(n)', 'FontSize', 14);
    saveas(hF, [path_to_pics '/png/k_intersys_BPSK.png']);
    saveas(hF, [path_to_pics '/fig/k_intersys_BPSK.fig']);
    [a b] = min(sum_dB);
    [c d] = min(min(sum_dB));
    if sum_dB_min > c
        sum_dB_min = c;
        n8_min = b(d);
        freq_min = d;
    end
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
end