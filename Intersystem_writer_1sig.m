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

Out_Type = 4;
% 1 - С/A
% 2 - P(Y)
% 3 - MBOC
% 4 - M-code
pic_num = 119;
pic_num = pic_num + 33*(Out_Type-1);
if Out_Type == 1
    filename_pref = 'CA_';
    signal_str = 'C/A';
elseif Out_Type == 2
    filename_pref = 'PY_';
    signal_str = 'P(Y)';
elseif Out_Type == 3
    filename_pref = 'MBOC_';
    signal_str = 'MBOC';
elseif Out_Type == 4
    filename_pref = 'MCode_';
    signal_str = 'M-Code';
end

for Signal_Type = 1:3
    for f_in = 1:fmax
        if (Signal_Type == BOCsin) 
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_sd_BoCsin_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Separation factor between received BoC<sub>sin</sub>(m, n) signals<br>and GPS L1 %s signals at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, signal_str, farr(f_in));
        elseif (Signal_Type == BOCcos)
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_sd_BoCcos_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Separation factor between received BoC<sub>cos</sub>(m, n) signals<br>and GPS L1 %s signals at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, signal_str, farr(f_in));
        elseif (Signal_Type == BPSK)
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_sd_BPSK.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Separation factor between received BPSK(n) signals<br>and GPS L1 %s signals</i><br><br> \n', pic_num, signal_str);
        end
        pic_num = pic_num + 1;
        if Signal_Type == BPSK
            break;
        end
    end
end
fprintf('</body> </html> \n');