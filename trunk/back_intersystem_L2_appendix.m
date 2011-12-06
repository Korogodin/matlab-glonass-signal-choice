%/**
% Скрипт составления приложения для межсистемных помех L2 - они нам
%*/

clear 
close all
clc

path_to_results = [pwd '/results/back_intersystem_L2/Common'];

load([path_to_results '/BackInterSysJam_BoCsin_L2.mat'], 'BackInterSysJam_BoCsin_L2');
load([path_to_results '/BackInterSysJam_BoCcos_L2.mat'], 'BackInterSysJam_BoCcos_L2');
load([path_to_results '/BackInterSysJam_BPSK_L2.mat'], 'BackInterSysJam_BPSK_L2');


fid_out = 1;
fid_out = fopen('temp.html', 'w');

BOCsin = 1; BOCcos = 2; BPSK = 3;
% Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

max_n_on_page = 27;
max_m_on_page = 40;
max_xn = 3;
max_xm = 2;
pril_num = 1;
fprintf(fid_out, '<!DOCTYPE html> \n');
fprintf(fid_out, '<html><head><title>InterIm L2 appendix</title> \n');
fprintf(fid_out, '<meta charset="utf-8"> </head> <body> \n');
for f_in = 1:fmax
    % 81 хорошо делится на 3, 80 на 2
%    1 2 3
%    4 5 6
%    При этом не стоит отображать пустые таблицы, если m,n выходит за
%    пределы полосы
    sum_dB_sin = BackInterSysJam_BoCsin_L2(:, :, f_in);
    sum_dB_cos = BackInterSysJam_BoCcos_L2(:, :, f_in);
    freq_index = f_in;
    for xm = 1:max_xm
        for xn = 1:max_xn
            
            m8_min_sin = 1 + max_m_on_page*max_xm - (max_m_on_page*(xm-1) + max_m_on_page);
            n8_min_cos = m8_min_sin;
            n8_min_sin = max_n_on_page*(xn-1) + 1;
            m8_min_cos = max_n_on_page*(xn-1) + 1 - 1;
            min_sum_m8_and_n8 =  min([n8_min_sin + m8_min_sin n8_min_cos + m8_min_cos]);
            
            if ( min_sum_m8_and_n8/8 > (freq_index + 2) ) || ( min_sum_m8_and_n8/8 > ((18-freq_index) +2) ) % Полоса
                continue;
            end
                fprintf(fid_out, '<b>Table I.%.0f</b>. <i>Mean separation factor for BOC<sub>sin</sub>(m, n) (blue) and BOC<sub>cos</sub>(m, n) (pink) L2 signals with f<sub>n</sub>=%.0f when receiving other systems signals, dB</i><br> \n', pril_num, farr(f_in));
                pril_num = pril_num + 1;
                fprintf(fid_out, '<table width="1000" cellpadding="0" cellspacing="0" border="0" bgcolor="#000000"><tr><td>\n'); 
                fprintf(fid_out, '<TABLE width="100%%" cellpadding="0" cellspacing="1" border="0" bgcolor="#000000">\n');
                fprintf(fid_out, '<FONT SIZE="1">');
                fprintf(fid_out, '<tr>');
                fprintf(fid_out, '<td bgcolor="#ffffff"><b>m \\ n</b></td>');
                for jj = 1:(max_n_on_page)
                    if (jj == (max_n_on_page))&&(xn == max_xn)
                        fprintf(fid_out, '<td bgcolor="#C0C0C0">');
                        fprintf(fid_out, '&nbsp;');
                    else
                        fprintf(fid_out, '<td bgcolor="#CCECFF">');
                        fprintf(fid_out, '<b>%.3f</b>', (max_n_on_page*(xn-1) + jj)/8);
                    end
                    fprintf(fid_out, '</td>');
                end
                fprintf(fid_out, '<td bgcolor="#FFFFFF">');
                fprintf(fid_out, '<b>n / n</b></td>');
                fprintf(fid_out, '</tr>\n');
                for j_m = 1:max_m_on_page    
                    fprintf(fid_out, '<tr>');
                    m8_sin = (1 + max_m_on_page*max_xm - (max_m_on_page*(xm-1) + j_m));
                    n8_cos = m8_sin;
                    fprintf(fid_out, '<td bgcolor="#CCECFF"><b>%.3f</b></td>', m8_sin/8); % m убывает от 10 до 0.125
                    for j_n = 1:max_n_on_page
                        n8_sin = max_n_on_page*(xn-1) + j_n;
                        m8_cos = max_n_on_page*(xn-1) + j_n - 1;
                        if n8_sin <= m8_sin
                            if isnan(sum_dB_sin(m8_sin, n8_sin))
                                fprintf(fid_out, '<td bgcolor="#CCECFF">');
                                fprintf(fid_out, '<center>-</center>');                            
                            else
                                fprintf(fid_out, '<td bgcolor="#CCECFF">');
                                fprintf(fid_out, '%.1f', sum_dB_sin(m8_sin, n8_sin));                        
                            end
                        else
                            if isnan(sum_dB_cos(m8_cos, n8_cos))
                                fprintf(fid_out, '<td bgcolor="#FFCCCC">');
                                fprintf(fid_out, '<center>-</center>');                            
                            else
                                fprintf(fid_out, '<td bgcolor="#FFCCCC">');
                                fprintf(fid_out, '%.1f', sum_dB_cos(m8_cos, n8_cos));                        
                            end
                        end
                        fprintf(fid_out, '</td>');
                    end
                    fprintf(fid_out, '<td bgcolor="FFCCCC"><b>%.3f</b></td>', n8_cos/8); 
                    fprintf(fid_out, '</tr>\n');                
                end

                fprintf(fid_out, '<tr>');
                fprintf(fid_out, '<td bgcolor="#ffffff">m / m</td>');
                for jj = 1:(max_n_on_page)
                    if (jj == (1))&&(xn == 1)
                        fprintf(fid_out, '<td bgcolor="#C0C0C0">');
                        fprintf(fid_out, '&nbsp;');
                    else
                        fprintf(fid_out, '<td bgcolor="#FFCCCC">');
                        fprintf(fid_out, '<b>%.3f</b>', ((max_n_on_page*(xn-1) + jj) - 1)/8);
                    end
                    fprintf(fid_out, '</td>');
                end
                fprintf(fid_out, '<td bgcolor="#FFFFFF">');
                fprintf(fid_out, '<b>m \\ n</b> </td>');
                fprintf(fid_out, '</tr>\n');            
                fprintf(fid_out, '</FONT>');
                fprintf(fid_out, '</TABLE>\n');
                fprintf(fid_out, '</td></tr></table><br>\n<br>\n');
    %             fprintf(fid_out, '&nbsp;<br>&nbsp;<br>&nbsp;<br>');
        end
    end
   f_in 
end

% Теперь для BPSK
sum_dB_BPSK = BackInterSysJam_BPSK_L2(:, :);
max_n_on_page = 40;
max_xn = 2;

for xn = 1:max_xn
% по x - частоты, по y - n
%    1 
%    2 
    fprintf(fid_out, '<b>Table I.%.0f</b>. <i>Mean separation factor for BPSK(n) L2 signals  when receiving other systems signals, dB</i><br> \n', pril_num);
    pril_num = pril_num + 1;
    fprintf(fid_out, '<table width="1000" cellpadding="0" cellspacing="0" border="0" bgcolor="#000000"><tr><td>\n'); 
    fprintf(fid_out, '<TABLE width="100%%" cellpadding="0" cellspacing="1" border="0" bgcolor="#000000">\n');
    fprintf(fid_out, '<FONT SIZE="1">\n');

    fprintf(fid_out, '<tr><td bgcolor="#ffffff"><b>n \\ f<sub>n</sub></b></td>');
    for f_in = 1:fmax    
        fprintf(fid_out, '<td bgcolor="#FFFFFF"><b>%.0f</b></td>', farr(f_in));
    end
    fprintf(fid_out, '<td bgcolor="#ffffff"><b>f<sub>n</sub> / n</b></td></tr>\n');

    for j_n = 1:max_n_on_page
        n8 = (xn-1)*max_n_on_page + j_n;
        fprintf(fid_out, '<tr><td bgcolor="#ffffff"><b>%.3f</b></td>', n8/8);
        for f_in = 1:fmax    
            if isnan(sum_dB_BPSK(n8, f_in))
                fprintf(fid_out, '<td bgcolor="#FFFFFF">');
                fprintf(fid_out, '<center>-</center>');                            
            else
                fprintf(fid_out, '<td bgcolor="#FFFFFF">%.1f</td>', sum_dB_BPSK(n8, f_in));
            end
        end
        fprintf(fid_out, '<td bgcolor="#ffffff"><b>%.3f</b></td></tr>\n', n8/8);
    end
    
    fprintf(fid_out, '<tr><td bgcolor="#ffffff"><b>n \\ f<sub>n</sub></b></td>');
    for f_in = 1:fmax    
        fprintf(fid_out, '<td bgcolor="#FFFFFF"><b>%.0f</b></td>', farr(f_in));
    end
    fprintf(fid_out, '<td bgcolor="#ffffff"><b>f<sub>n</sub> / n</b></td></tr>\n');
                    
            fprintf(fid_out, '</FONT>');
            fprintf(fid_out, '</TABLE>\n');
            fprintf(fid_out, '</td></tr></table><br>\n<br>\n');
%             fprintf(fid_out, '&nbsp;<br>&nbsp;<br>&nbsp;<br>');
end


fprintf(fid_out, '</body> </html> \n');

if fid_out ~= 1
    fclose(fid_out);
end
    