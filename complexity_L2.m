%/**
% Данный скрипт:
% 1) либо рассчитывает базовый показатель сложности - только по
% полосе и многопиковости без учета всяких спецблоков и т.п., создает 
% соответствующую wiki-таблицу
% 2) либо формирует html-документ с рисунками и подписями
%
%@param Signal_Type задает BoCsin, BoCcos или BPSK
%@param Table_Type задает wiki или html
%*/

clear 
close all
clc

n8max = 80;
m8max = 80;
farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты

path_to_stuff = '/complexity_stuff/L2'; 
% Путь должен существовать, а в нем каталоги png и fig
% Туда будем класть картинки, туда же сдедует положить html
path_to_results = [pwd '/results/complexity_L2'];
path_to_ro = [pwd '/ro'];

% С помощью данной секции можно создать пустые массивы
% Complexity_BoCsin = nan(m8max, n8max, fmax);
% Complexity_BoCcos = nan(m8max, n8max, fmax);
% Complexity_BPSK = nan(n8max, fmax);
% save([path_to_results '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
% save([path_to_results '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
% save([path_to_results '/Complexity_BPSK.mat'], 'Complexity_BPSK');
% return
load([path_to_results '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
load([path_to_results '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
load([path_to_results '/Complexity_BPSK.mat'], 'Complexity_BPSK');

Signals_L1; % Параметры приемлимых сигналов

% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin; 2 - BOCcos; 3 - BPSK.

% Тип формируемого:
% wiki - таблички для wiki
% html - html-файл с рисунками АКФ и подписями
T_wiki = 1; T_html = 2;
Table_Type = 0; % 1 - wiki, 2 - html

load([path_to_ro '/Td.mat']);

if Table_Type == T_wiki
    fprintf('{| class="wikitable sortable" border="1" \n');
elseif Table_Type == T_html
    fprintf('<!DOCTYPE html> \n');
    fprintf('<html><head><title>!DOCTYPE</title> \n');
    fprintf('<meta charset="utf-8"> </head> <body> \n');
end


if Signal_Type == BOCsin
    Nsig = size(BoCsin_Freq_L1_num, 1);
    Sig_Arr = BoCsin_Freq_L1_num;
    if Table_Type == T_wiki
        fprintf('|+ Partial complexity of BoC<sub>sin</sub> L2 signals, scores\n');    
    end
elseif Signal_Type == BOCcos
    Nsig = size(BoCcos_Freq_L1_num, 1);
    Sig_Arr = BoCcos_Freq_L1_num;    
    if Table_Type == T_wiki
        fprintf('|+ Partial complexity of BoC<sub>cos</sub> L2 signals, scores\n');    
    end  
elseif Signal_Type == BPSK
    Nsig = size(BPSK_Freq_L1_num, 1);
    Sig_Arr = zeros(Nsig, 3);
    Sig_Arr(:, 2:3) = BPSK_Freq_L1_num;
    if Table_Type == T_wiki
        fprintf('|+ Partial complexity of BPSK L2 signals, scores\n');    
    end
end


if Table_Type == T_wiki
    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!Band \n');
    fprintf('!Second peak \n');
    fprintf('!Sum');
end


for n8 = 1:n8max
    for m8 = 1:m8max
        if m8 < n8
            if Signal_Type ~= BPSK
                continue;
            end
        end
        n = n8 / 8;
        m = m8 / 8 * (Signal_Type ~= BPSK);

        if Table_Type == T_wiki
            fprintf('\n|- align="center"\n');
            if Signal_Type == BOCsin 
                fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            elseif Signal_Type == BOCcos
                fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            elseif Signal_Type == BPSK
                fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
            end
        end
            
        % Число балов сложности за полосу сигнала
        SignalBand = 2*(n + m);
        Score_SignalBand = SignalBand * 3;
        if Table_Type == T_wiki
            fprintf('|| %.1f (%.0f<math>f_b</math>) ', round(Score_SignalBand*10)/10, round(SignalBand*10)/10 );
        end

        % Открываем файл АКФ указанного сигнала
        ro_our = get_ro(m, n, Signal_Type, path_to_ro);
        N_ro = length(ro_our);
        abs_ro = abs(ro_our);
        diff_abs_ro = diff(abs_ro);

        [a, b] = max(abs_ro);
        sign_old = -1;
        max2 = 0;
        for j = (b+1):(N_ro-1)
            if (sign_old > 0)&&(diff_abs_ro(j) < 0)
                max2 = abs_ro(j);
                break;
            end
            sign_old = diff_abs_ro(j);
        end
        if max2 > 0.01
            Score_Peak = max2*12;
        else
            Score_Peak = 0;
        end
        
%         if n == round(n)
%             Score_n = 0;
%         elseif 2*n == round(2*n) 
%             Score_n = 1;
%         elseif 4*n == round(4*n)    
%             Score_n = 4;
%         elseif 8*n == round(n)        
%             Score_n = 8;
%         else
%             Score_n = 16;
%         end

       for j = 1:50
        if (n*1023*j) == round(n*1023*j)
            break;
        end
       end
       Score_n = j;
       
       for j = 1:50
        if (m*1023*j) == round(m*1023*j)
            break;
        end
       end
       Score_m = j;       
        
       Score_mn = max([Score_m Score_n])*3;
%         if m/n == round(m/n)
%             Score_mn = 0;
%         elseif 2*(m/n) == round(2*m/n) 
%             Score_mn = 2;
%         elseif 4*(m/n) == round(4*m/n)    
%             Score_mn = 4;
%         elseif 8*(m/n) == round(*8m/n)        
%             Score_mn = 8;
%         else
%             Score_mn = 16;
%         end

%         Score_m = 0;
%         if m == round(m)
% %             Score_m = 0;
%         elseif 2*m == round(m) 
% %             Score_m = 1;
%         elseif 4*m == round(m)    
% %             Score_m = 4;
%         elseif 8*m == round(m)        
% %             Score_m = 8;
%         else
% %             Score_m = 16;
%         end
        
        if Table_Type == T_wiki
            fprintf('|| %.1f (%.0f%%) ', round(Score_Peak*10)/10, round(max2*100*(max2>0.01)) );
            fprintf('|| %.1f ', round((Sum_Score)*10)/10);
        end
        for freq_index = 1:fmax
            Score_5 = 10 * ((freq_index ~= 1)&&(freq_index ~= 6)&&(freq_index ~= 11)&&(freq_index ~= 16));
            Score_Center = 3*abs(farr(freq_index) - (max(farr) + min(farr))/2);
%             Sum_Score = Score_SignalBand + Score_Peak + Score_n + Score_m
%             + Score_5;
            Sum_Score = Score_SignalBand + Score_Peak + Score_mn + Score_5 + Score_Center;
            if Signal_Type == BOCsin
                Complexity_BoCsin(m8, n8, freq_index) = Sum_Score;
            elseif Signal_Type == BOCcos
                Complexity_BoCcos(m8, n8, freq_index) = Sum_Score;
            elseif Signal_Type == BPSK
                Complexity_BPSK(n8, freq_index) = Sum_Score;
            end
        end 
        if Signal_Type == BPSK
            break;
        end
    end
    n
    save([path_to_results '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
    save([path_to_results '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
    save([path_to_results '/Complexity_BPSK.mat'], 'Complexity_BPSK');
end
    
%     hF = 0;
% 
%     N_co = (length(ro_our) - 1) / 2;
%     hF = figure(hF + 1);
%     tau = (-N_co:1:N_co)*Td*1e6;
%     plot(tau, ro_our);
%     xlabel('\tau, \mu{s}', 'FontSize', 14)
%     ylabel('\rho(\tau)', 'FontSize', 14)
%     grid on;
%     set(gca, 'FontSize', 14)
%     if Signal_Type == BOCsin
%         title(sprintf('{\\rho}({\\tau}) for BOC_{sin}(%s, %s)', sprintf('%.3f', m), sprintf('%.3f', n)), 'FontSize', 14);
%     elseif Signal_Type == BOCcos
%         title(sprintf('{\\rho}({\\tau}) for BOC_{cos}(%s, %s)', sprintf('%.3f', m), sprintf('%.3f', n)), 'FontSize', 14);
%     elseif Signal_Type == BPSK
%         title(sprintf('{\\rho}({\\tau}) for BPSK(%s)', sprintf('%.3f', n)), 'FontSize', 14);
%     end
%     drawnow
%     if (Signal_Type == BOCsin) 
%         if Table_Type == T_html
%             fprintf('<img src="%s" width="530" /><br> \n', ['png/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
%             fprintf('<b>Рисунок 1.7.zzsin%.0f</b> -  <i>Автокорреляционная функция<br>сигнала BoC<sub>sin</sub>(%.3f, %.3f)</i><br><br> \n', i, m, n);
%         end
%         saveas(hF, [pwd path_to_stuff '/png/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
%         saveas(hF, [pwd path_to_stuff '/fig/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').fig']);
%     elseif (Signal_Type == BOCcos)
%         if Table_Type == T_html
%             fprintf('<img src="%s" width="530" /><br> \n', ['png/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
%             fprintf('<b>Рисунок 1.7.zzcos%.0f</b> -  <i>Автокорреляционная функция<br>сигнала BoC<sub>cos</sub>(%.3f, %.3f)</i><br><br> \n', i, m, n);
%         end
%         saveas(hF, [pwd path_to_stuff '/png/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
%         saveas(hF, [pwd path_to_stuff '/fig/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').fig']);
%     elseif (Signal_Type == BPSK)
%         if Table_Type == T_html
%             fprintf('<img src="%s" width="530" /><br> \n', ['png/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', n) ').png']);
%             fprintf('<b>Рисунок 1.7.zzBPSK%.0f</b> -  <i>Автокорреляционная функция<br>сигнала BPSK(%.3f)</i><br><br> \n', i, n);
%         end
%         saveas(hF, [pwd path_to_stuff '/png/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
%         saveas(hF, [pwd path_to_stuff '/fig/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').fig']);
%     end
% % end
% if Table_Type == T_wiki
%     fprintf('\n|} \n');
% elseif Table_Type == T_html
%     fprintf('</body> </html> \n');
% end

save([path_to_results '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
save([path_to_results '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
save([path_to_results '/Complexity_BPSK.mat'], 'Complexity_BPSK');
