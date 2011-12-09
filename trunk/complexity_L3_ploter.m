%/**
% Скрипт отрисовки графики для межсистемных помех
%*/

clear 
close all
clc

path_to_pics = [pwd '/stuff/complexity_stuff/L3']; 
% Путь должен существовать, а в нем каталоги png и fig
% Туда будем класть картинки, туда же сдедует положить html
path_to_results = [pwd '/results/complexity_L3'];
path_to_ro = [pwd '/ro'];

% С помощью данной секции можно создать пустые массивы
load([path_to_results '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
load([path_to_results '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
load([path_to_results '/Complexity_BPSK.mat'], 'Complexity_BPSK');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1164:1184; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

DoPlot = 0;

CLimm = [0 200];
sum_dB_min = 999; sum_dB_max = -999; n8_min = -1; n8_max = -1; 
if Signal_Type == BOCsin
    hF = 0;
    for f_in = 1:fmax
        sum_dB = Complexity_BoCsin(:, :, f_in);
        if DoPlot
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 10)
            ylabel('m', 'FontSize', 10)      
            colorbar
            set(gca, 'CLim', CLimm)
            set(gca, 'FontSize', 10)
            set(hF, 'PaperPosition', [0 0 14 11]);   
            title(sprintf('Complexity for BOC_{sin} L3 signals at f_n = %.0f, scores', farr(f_in)), 'FontSize', 10);
            saveas(hF, [path_to_pics '/png/Compl_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/Compl_BoCsin_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    fprintf('Maximum: %.1f for %s\n', sum_dB_max, ['BoCsin(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.1f for %s\n', sum_dB_min, ['BoCsin(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
elseif Signal_Type == BOCcos
    hF = 0;
    for f_in = 1:fmax
        sum_dB = Complexity_BoCcos(:, :, f_in);
        if DoPlot
            hF = figure(hF+1);
            pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
            xlabel('n', 'FontSize', 10)
            ylabel('m', 'FontSize', 10)      
            colorbar
            set(gca, 'CLim', CLimm)
            set(gca, 'FontSize', 10)
            set(hF, 'PaperPosition', [0 0 14 11]);   
            title(sprintf('Complexity for BOC_{cos} L3 signals at f_n = %.0f, scores', farr(f_in)), 'FontSize', 10);
            saveas(hF, [path_to_pics '/png/Compl_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.png']);
            saveas(hF, [path_to_pics '/fig/Compl_BoCcos_at_' sprintf('%.0f', farr(f_in)) '.fig']);
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
    fprintf('Maximum: %.1f for %s\n', sum_dB_max, ['BoCcos(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.1f for %s\n', sum_dB_min, ['BoCcos(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    

elseif Signal_Type == BPSK
    
    sum_dB = Complexity_BPSK;
    
    if DoPlot
        hF = 0; hF = figure(hF+1);
        pcolor((1:80)/8, [farr farr(fmax)+1], [sum_dB(1:80,1:fmax), nan(80,1)]');
        xlabel('n', 'FontSize', 10)
        ylabel('f_{n}', 'FontSize', 10)      
        colorbar
        set(gca, 'CLim', CLimm)
        set(gca, 'FontSize', 10)
        set(hF, 'PaperPosition', [0 0 14 11]);   
        title(sprintf('Complexity for BPSK L3 signals, scores'), 'FontSize', 10);
        saveas(hF, [path_to_pics '/png/Compl_BPSK.png']);
        saveas(hF, [path_to_pics '/fig/Compl_BPSK.fig']);
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
    fprintf('Maximum: %.1f for %s\n', sum_dB_max, ['BPSK(' sprintf('%.3f', n8_max/8) ') at ' sprintf('%.0f', farr(freq_max))] );    
    fprintf('Minimum: %.1f for %s\n', sum_dB_min, ['BPSK(' sprintf('%.3f', n8_min/8) ') at ' sprintf('%.0f', farr(freq_min))] );    
end