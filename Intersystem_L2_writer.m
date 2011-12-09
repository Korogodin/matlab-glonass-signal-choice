%/**
% Скрипт создает html-документ с полученными картинками
%*/

clear 
close all
clc

BOCsin = 1; BOCcos = 2; BPSK = 3;

farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

    fprintf('<!DOCTYPE html> \n');
    fprintf('<html><head><title>!DOCTYPE</title> \n');
    fprintf('<meta charset="utf-8"> </head> <body> \n');   

pic_num = 633;

for Signal_Type = 1:3
    for f_in = 1:fmax
        if (Signal_Type == BOCsin) 
            fprintf('<img src="%s" /><br> \n', ['png/k_intersys_BoCsin_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Figure 1.7.%.0f</b> -  <i>Coefficient k<sub>intersys</sub> for BoC<sub>sin</sub>(m, n) L2 signals<br>at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, farr(f_in));
        elseif (Signal_Type == BOCcos)
            fprintf('<img src="%s" /><br> \n', ['png/k_intersys_BoCcos_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Figure 1.7.%.0f</b> -  <i>Coefficient k<sub>intersys</sub> for BoC<sub>cos</sub>(m, n) L2 signals<br>at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, farr(f_in));
        elseif (Signal_Type == BPSK)
            fprintf('<img src="%s" /><br> \n', 'png/k_intersys_BPSK.png');
            fprintf('<b>Figure 1.7.%.0f</b> -  <i>Coefficient k<sub>intersys</sub> for BPSK(n) signals<br>and GPS L2 signals</i><br><br> \n', pic_num);
        end
        pic_num = pic_num + 1;
        if Signal_Type == BPSK
            break;
        end
    end
end
fprintf('</body> </html> \n');