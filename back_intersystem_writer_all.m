%/**
% Скрипт создает html-документ с полученными картинками
%*/

clear 
close all
clc

BOCsin = 1; BOCcos = 2; BPSK = 3;

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

    fprintf('<!DOCTYPE html> \n');
    fprintf('<html><head><title>!DOCTYPE</title> \n');
    fprintf('<meta charset="utf-8"> </head> <body> \n');   

Out_Type = 1;    
pic_num = 416;
pic_num = pic_num + 33*(Out_Type-1);
% fprintf('%.0f - %.0f : %.0f - %.0f : %.0f \n', pic_num, pic_num+15, pic_num+16, pic_num+31, pic_num + 32);
% return

filename_pref = 'All_';


for Signal_Type = 1:3
    for f_in = 1:fmax
        if (Signal_Type == BOCsin) 
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_sd_BoCsin_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Mean separation factor between BoC<sub>sin</sub>(m, n) signals<br>and received GPS L1 and GLONASS FDMA L1 signals at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, farr(f_in));
        elseif (Signal_Type == BOCcos)
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_sd_BoCcos_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Mean separation factor between BoC<sub>cos</sub>(m, n) signals<br>and received GPS L1 and GLONASS FDMA L1 signals at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, farr(f_in));
        elseif (Signal_Type == BPSK)
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_sd_BPSK.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Mean separation factor between BPSK(n) signals<br>and received GPS L1 and GLONASS FDMA L1 signals</i><br><br> \n', pic_num);
        end
        pic_num = pic_num + 1;
        if Signal_Type == BPSK
            break;
        end
    end
end
fprintf('</body> </html> \n');