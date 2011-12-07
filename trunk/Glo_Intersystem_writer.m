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

pic_num = 85;
    
Out_Type = 2;
% 1 - OF
% 2 - SF
% 3 - OF+SF
if Out_Type == 1
    filename_pref = 'OF_';
    signal_str = ['L1OF'];
elseif Out_Type == 2
    filename_pref = 'SF_';
    signal_str = 'L1SF';
elseif Out_Type == 3
    filename_pref = 'Both_';
    signal_str = 'L1OF+L1SF';
end

for Signal_Type = 1:3
    for f_in = 1:fmax
        if (Signal_Type == BOCsin) 
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_intersys_BoCsin_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Mean separation factor between BoC<sub>sin</sub>(m, n) signals<br>and GLONASS %s signals at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, signal_str, farr(f_in));
        elseif (Signal_Type == BOCcos)
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_intersys_BoCcos_at_' sprintf('%3.0f', farr(f_in)) '.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Mean separation factor between BoC<sub>cos</sub>(m, n) signals<br>and GLONASS %s signals at f<sub>n</sub>=%.0f</i><br><br> \n', pic_num, signal_str, farr(f_in));
        elseif (Signal_Type == BPSK)
            fprintf('<img src="%s" /><br> \n', ['png/' filename_pref 'k_intersys_BPSK.png']);
            fprintf('<b>Рисунок 1.7.%.0f</b> -  <i>Mean separation factor between BPSK(n) signals<br>and GLONASS %s signals</i><br><br> \n', pic_num, signal_str);
        end
        pic_num = pic_num + 1;
        if Signal_Type == BPSK
            break;
        end
    end
end
fprintf('</body> </html> \n');