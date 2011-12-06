%/**
% Скрипт составления приложения для Точности
%*/

clear 
close all
clc

path_to_results = [pwd '/results/beta'];

load([path_to_results '/Beta_BoCsin.mat'], 'Beta_BoCsin');
load([path_to_results '/Beta_BoCcos.mat'], 'Beta_BoCcos');
load([path_to_results '/Beta_BPSK.mat'], 'Beta_BPSK');


fid_out = 1;
fid_out = fopen('temp.html', 'w');

BOCsin = 1; BOCcos = 2; BPSK = 3;
% Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

% farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

max_n_on_page = 27;
max_m_on_page = 40;
max_xn = 3;
max_xm = 2;
pril_num = 1;
fprintf(fid_out, '<!DOCTYPE html> \n');
fprintf(fid_out, '<html><head><title>Radioastronomy appendix</title> \n');
fprintf(fid_out, '<meta charset="utf-8"> </head> <body> \n');
for f_in = 1:1
    % 81 хорошо делится на 3, 80 на 2
%    1 2 3
%    4 5 6
    sum_dB_sin = 1./Beta_BoCsin(:, :)*1e9;
    sum_dB_cos = 1./Beta_BoCcos(:, :)*1e9;
    for xm = 1:max_xm
        for xn = 1:max_xn
            fprintf(fid_out, '<b>Table B.%.0f</b>. <i>Inverse effective width of signal spectrum for BOC<sub>sin</sub>(m, n) (blue) and BOC<sub>cos</sub>(m, n) (pink) signals, ns</i><br> \n', pril_num);
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
                        fprintf(fid_out, '<td bgcolor="#CCECFF">');
                        fprintf(fid_out, '%.1f', sum_dB_sin(m8_sin, n8_sin));                        
                    else
                        fprintf(fid_out, '<td bgcolor="#FFCCCC">');
                        fprintf(fid_out, '%.1f', sum_dB_cos(m8_cos, n8_cos));                        
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
end

% Теперь для BPSK
sum_dB_BPSK = 1./Beta_BPSK(:)*1e9;
max_n_on_row = 20;
xn_max = 4;

fprintf(fid_out, '<b>Table B.%.0f</b>. <i>Inverse effective width of signal spectrum for BPSK(n) signals, ns</i><br> \n', pril_num);
pril_num = pril_num + 1;
fprintf(fid_out, '<table width="1000" cellpadding="0" cellspacing="0" border="0" bgcolor="#000000"><tr><td>\n'); 
fprintf(fid_out, '<TABLE width="100%%" cellpadding="0" cellspacing="1" border="0" bgcolor="#000000">\n');
fprintf(fid_out, '<FONT SIZE="1">\n');

for xn = 1:xn_max
    fprintf(fid_out, '<tr><td bgcolor="#ffffff"><b>n</b></td>');
    for jj = 1:max_n_on_row    
        n8 = jj + max_n_on_row*(xn - 1);
        fprintf(fid_out, '<td bgcolor="#FFFFFF">%.3f</td>', n8/8);
    end
    fprintf(fid_out, '<td bgcolor="#ffffff"><b>n</b></td></tr>\n');
    
    fprintf(fid_out, '<tr><td bgcolor="#ffffff"><b>1/beta</b></td>');
    for jj = 1:max_n_on_row
        n8 = jj + max_n_on_row*(xn - 1);
        fprintf(fid_out, '<td bgcolor="#FFFFFF">%.2f</td>', sum_dB_BPSK(n8));
    end
    fprintf(fid_out, '<td bgcolor="#ffffff"><b>1/beta</b></td></tr>\n');
    
    if xn ~= xn_max
        fprintf(fid_out, '<tr><td bgcolor="#ffffff"><b>&nbsp;</b></td>');
        for jj = 1:max_n_on_row
            fprintf(fid_out, '<td bgcolor="#FFFFFF"><b>&nbsp;</b></td>');
        end
        fprintf(fid_out, '<td bgcolor="#ffffff"><b>&nbsp;</b></td></tr>\n');
    end

end
fprintf(fid_out, '</FONT>');
fprintf(fid_out, '</TABLE>\n');
fprintf(fid_out, '</td></tr></table><br>\n<br>\n');
%     fprintf(fid_out, '&nbsp;<br>&nbsp;<br>&nbsp;<br>');


fprintf(fid_out, '</body> </html> \n');

if fid_out ~= 1
    fclose(fid_out);
end
    