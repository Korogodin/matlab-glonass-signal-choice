%/**
% Скрипт отрисовки графики для излучения в РА
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/radioastronomy'];
path_to_pics = [pwd '/stuff/radioastronomy_stuff'];

load([path_to_results '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
load([path_to_results '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
load([path_to_results '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

n_min = -1; m_min = -1; freq_min = -1; sum_dB_min = 999; 
if Signal_Type == BOCsin
    
    hF = 0;
    for f_in = 1:fmax
            hF = figure(hF+1);
            sum_dB = Radioastronomy_BoCsin(:, :, f_in);
            pcolor((1:80)/8, (1:80)/8, sum_dB);
            xlabel('n', 'FontSize', 10)
            ylabel('m', 'FontSize', 10)      
            colorbar
            set(gca, 'CLim', [-115 -55])
            set(gca, 'FontSize', 10)
            set(hF, 'PaperPosition', [0 0 13 11]);
            title(sprintf('PSDmax in RA for BOC_{sin}, f_n = %.0f', farr(f_in)), 'FontSize', 10);
            drawnow
            saveas(hF, [path_to_pics '/png/PSD_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/PSD_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    hF = 0;
    for f_in = 1:fmax
            sum_dB = Radioastronomy_BoCcos(:, :, f_in);
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 10)
            ylabel('m', 'FontSize', 10)      
            colorbar
            set(gca, 'CLim', [-115 -55])
            set(gca, 'FontSize', 10)
            set(hF, 'PaperPosition', [0 0 13 11]);           
            title(sprintf('PSDmax in RA for BOC_{cos}, f_n = %.0f', farr(f_in)), 'FontSize', 10);
            drawnow
            saveas(hF, [path_to_pics '/png/PSD_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/PSD_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    
    hF = 0;
    sum_dB = Radioastronomy_BPSK(:, :);
    hF = figure(hF+1);
    pcolor((1:80)/8, [farr farr(fmax)+1], [sum_dB(1:80,1:fmax), nan(80,1)]');
    xlabel('n', 'FontSize', 10)
    ylabel('f_{n}', 'FontSize', 10)      
    colorbar
    set(gca, 'CLim', [-115 -55])
    set(gca, 'FontSize', 10)
    set(hF, 'PaperPosition', [0 0 13 11]);          
    title('PSDmax in RA for BPSK(n)', 'FontSize', 10);
    drawnow
    saveas(hF, [path_to_pics '/png/PSD_BPSK.png']);
    saveas(hF, [path_to_pics '/fig/PSD_BPSK.fig']);
    [a b] = min(sum_dB);
    [c d] = min(min(sum_dB));
    if sum_dB_min > c
        sum_dB_min = c;
        n8_min = b(d);
        freq_min = d;
    end
    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
end