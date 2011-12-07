%/**
% Скрипт отрисовки графики для параметра 1/beta
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/beta'];
path_to_pics = [pwd '/stuff/beta'];

load([path_to_results '/Beta_BoCsin.mat'], 'Beta_BoCsin');
load([path_to_results '/Beta_BoCcos.mat'], 'Beta_BoCcos');
load([path_to_results '/Beta_BPSK.mat'], 'Beta_BPSK');

BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 2; %1 - BOCsin, 2 - BOCcos, 3 - BPSK

% farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

Color_Lim = [-78 -55];

n8_min = -1; m8_min = -1; sum_dB_min = 1e9; 
n8_max = -1; m8_max = -1; sum_dB_max = -1e9; 
if Signal_Type == BOCsin
    hF = 0;
        hF = figure(hF+1);
        sum_dB = 1./Beta_BoCsin(:, :)*1e9;
%         sum_dB = 10*log10(1./Beta_BoCsin(:, :));
        pcolor((1:n8max)/8, (1:m8max)/8, sum_dB(:,:));
        xlabel('n', 'FontSize', 10)
        ylabel('m', 'FontSize', 10)      
        colorbar
%         set(gca, 'CLim', Color_Lim)
        set(gca, 'FontSize', 10)
        set(hF, 'PaperPosition', [0 0 13 11]);
        title('1 / \beta for BOC_{sin}, ns', 'FontSize', 10);
        drawnow
%         saveas(hF, [path_to_pics '/png/OneDevBeta_BoCsin.png']);
%         saveas(hF, [path_to_pics '/fig/OneDevBeta_BoCsin.fig']);
        hF = figure(hF+1);
        mesh((1:n8max)/8, (1:m8max)/8, sum_dB(:,:));
        xlabel('n', 'FontSize', 10)
        ylabel('m', 'FontSize', 10)      
        title('1 / \beta for BOC_{sin}, ns', 'FontSize', 10);
        zlim([min(min(sum_dB)) - 10 max(max(sum_dB)) - 10])
        saveas(hF, [path_to_pics '/png/OneDevBeta_Mesh_BoCsin.png']);
        saveas(hF, [path_to_pics '/fig/OneDevBeta_Mesh_BoCsin.fig']);
        
        [a b] = min(sum_dB);
        [c d] = min(min(sum_dB));
            sum_dB_min = c;
            m8_min = b(d);
            n8_min = d;

        [a b] = max(sum_dB);
        [c d] = max(max(sum_dB));
            sum_dB_max = c;
            m8_max = b(d);
            n8_max = d;

    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BoCsin(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ')'] );    
    fprintf('Maximum: %.3f for %s\n', sum_dB_max, ['BoCsin(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ')'] );    

elseif Signal_Type == BOCcos
    hF = 0;
        sum_dB = 1./Beta_BoCcos(:, :)*1e9;
%         sum_dB = 10*log10(1./Beta_BoCcos(:, :));
        hF = figure(hF+1);
        pcolor((1:80)/8, (1:80)/8, sum_dB(1:80,1:80));
        xlabel('n', 'FontSize', 10)
        ylabel('m', 'FontSize', 10)      
        colorbar
%         set(gca, 'CLim', Color_Lim)
        set(gca, 'FontSize', 10)
        set(hF, 'PaperPosition', [0 0 14 11]);           
        title('1 / {\beta} for BOC_{cos}, ns', 'FontSize', 10);
        drawnow
%         saveas(hF, [path_to_pics '/png/OneDevBeta_BoCcos.png']);
%         saveas(hF, [path_to_pics '/fig/OneDevBeta_BoCcos.fig']);

        hF = figure(hF+1);
        mesh((1:n8max)/8, (1:m8max)/8, sum_dB(:,:));
        xlabel('n', 'FontSize', 10)
        ylabel('m', 'FontSize', 10)      
        title('1 / \beta for BOC_{cos}, ns', 'FontSize', 10);
        zlim([min(min(sum_dB)) - 10 max(max(sum_dB)) - 10])
        saveas(hF, [path_to_pics '/png/OneDevBeta_Mesh_BoCcos.png']);
        saveas(hF, [path_to_pics '/fig/OneDevBeta_Mesh_BoCcos.fig']);
        [a b] = min(sum_dB);
        [c d] = min(min(sum_dB));
            sum_dB_min = c;
            m8_min = b(d);
            n8_min = d;
    
        [a b] = max(sum_dB);
        [c d] = max(max(sum_dB));
            sum_dB_max = c;
            m8_max = b(d);
            n8_max = d;

    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BoCcos(' sprintf('%.3f', m8_min/8) ', ' sprintf('%.3f', n8_min/8) ')'] );    
    fprintf('Maximum: %.3f for %s\n', sum_dB_max, ['BoCcos(' sprintf('%.3f', m8_max/8) ', ' sprintf('%.3f', n8_max/8) ')'] );    

elseif Signal_Type == BPSK

    hF = 0;
    sum_dB = 1./Beta_BPSK(:)*1e9;
    hF = figure(hF+1);
    plot( (1:n8max)/8, sum_dB(1:n8max) );
    xlabel('n', 'FontSize', 10)
    ylabel('1/\beta, ns', 'FontSize', 10)      
    grid on
    set(gca, 'FontSize', 10)
    set(hF, 'PaperPosition', [0 0 14 11]);          
    title('1 / {\beta} for BPSK(n)', 'FontSize', 10);
    drawnow
    saveas(hF, [path_to_pics '/png/OneDevBeta_BPSK.png']);
    saveas(hF, [path_to_pics '/fig/OneDevBeta_BPSK.fig']);
    [a b] = min(sum_dB);
        sum_dB_min = a;
        n8_min = b;

    [a b] = max(sum_dB);
        sum_dB_max = a;
        n8_max = b;

    fprintf('Minimum: %.3f for %s\n', sum_dB_min, ['BPSK(' sprintf('%.3f', n8_min/8) ')'] );    
    fprintf('Maximum: %.3f for %s\n', sum_dB_max, ['BPSK(' sprintf('%.3f', n8_max/8) ')'] );    
end