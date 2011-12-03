%/**
% Данный скрипт рассчитывает коэффициенты спектрального разделения между
% сигналами полного выборочного множества и сигналами системы GPS при том 
% условии, что фронтенд настроен на прием наших сигналов (полоса фронтенда 
% равна полосе сигнала по нулям в теории; на практике в данном скрипте
% фильтр максимально широкий для выделенного нам частотного диапазона, но
% при этом симметричный относительно несущей сигнала). Результат вычисления 
% коэффициентов частотного разделения сохраняется в массивах 
% InterSysJam_*** с индексами (8*m, 8*n, freq_ind) для BoC-сигналов или 
% (8*n, freq_ind) для BPSK-сигналов, которые сохраняются в соответствующих 
% mat-файлах. 
%@param Signal_Type задает BoCsin, BoCcos или BPSK
%@param m8max, n8max - максимальное значение индексов 8*m, 8*n. Пробегаем
%по n и m  с шагом 0.125
%@param farr - массив номированных несущих частот
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/intersystem_L2'];

n8max = 80;
m8max = 80;
farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты

% Signals_L2; % Параметры приемлимых сигналов

% С помощью данной секции можно создать пустые массивы
% InterSysJam_BoCsin_L2_BoC_10_5 = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L2_BoC_10_5 = nan(m8max, n8max, fmax);
% InterSysJam_BPSK_L2_BoC_10_5 = nan(n8max, fmax);
% save([path_to_results '/InterSysJam_BoCsin_L2_BoC_10_5.mat'], 'InterSysJam_BoCsin_L2_BoC_10_5');
% save([path_to_results '/InterSysJam_BoCcos_L2_BoC_10_5.mat'], 'InterSysJam_BoCcos_L2_BoC_10_5');
% save([path_to_results '/InterSysJam_BPSK_L2_BoC_10_5.mat'], 'InterSysJam_BPSK_L2_BoC_10_5');
% 
% InterSysJam_BoCsin_L2_BoC_0_1 = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L2_BoC_0_1 = nan(m8max, n8max, fmax);
% InterSysJam_BPSK_L2_BoC_0_1 = nan(n8max, fmax);
% save([path_to_results '/InterSysJam_BoCsin_L2_BoC_0_1.mat'], 'InterSysJam_BoCsin_L2_BoC_0_1');
% save([path_to_results '/InterSysJam_BoCcos_L2_BoC_0_1.mat'], 'InterSysJam_BoCcos_L2_BoC_0_1');
% save([path_to_results '/InterSysJam_BPSK_L2_BoC_0_1.mat'], 'InterSysJam_BPSK_L2_BoC_0_1');
% 
% InterSysJam_BoCsin_L2_BoC_0_10 = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L2_BoC_0_10 = nan(m8max, n8max, fmax);
% InterSysJam_BPSK_L2_BoC_0_10 = nan(n8max, fmax);
% save([path_to_results '/InterSysJam_BoCsin_L2_BoC_0_10.mat'], 'InterSysJam_BoCsin_L2_BoC_0_10');
% save([path_to_results '/InterSysJam_BoCcos_L2_BoC_0_10.mat'], 'InterSysJam_BoCcos_L2_BoC_0_10');
% save([path_to_results '/InterSysJam_BPSK_L2_BoC_0_10.mat'], 'InterSysJam_BPSK_L2_BoC_0_10');


load([path_to_results '/InterSysJam_BoCsin_L2_BoC_10_5.mat']);
load([path_to_results '/InterSysJam_BoCcos_L2_BoC_10_5.mat']);
load([path_to_results '/InterSysJam_BPSK_L2_BoC_10_5.mat']);
load([path_to_results '/InterSysJam_BoCsin_L2_BoC_0_1.mat']);
load([path_to_results '/InterSysJam_BoCcos_L2_BoC_0_1.mat']);
load([path_to_results '/InterSysJam_BPSK_L2_BoC_0_1.mat']);
load([path_to_results '/InterSysJam_BoCsin_L2_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BoCcos_L2_BoC_0_10.mat']);
load([path_to_results '/InterSysJam_BPSK_L2_BoC_0_10.mat']);

% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

load([pwd '/ro/Td.mat']);


% АКФ сигнала GPS BoC(10, 5)
ro_GPS_BoC_10_5 = get_ro(10, 5, BOCsin, path_to_ro );

% АКФ сигнала GPS BPSK(1)
ro_GPS_BoC_0_1 = get_ro(0, 1, BPSK, path_to_ro );

% АКФ сигнала GPS BPSK(10)
ro_GPS_BoC_0_10 = get_ro(0, 10, BPSK, path_to_ro );

f_intro = 1227.6e6; % Частота мешающих сигналов

N_ro_old = 0;
f_index_old = 0;
n_plus_m_old = 0;
Hd_old = 0;
% Этот блок используется, если хотим пробежать только по списку сигналов
% for i = 1:size(BPSK_Freq_L2_num, 1)
% 
% % n = BoCcos_Freq_L2_num(i, 2); n8 = n*8;
% % m = BoCcos_Freq_L2_num(i, 1); m8 = m*8;
% % freq = BoCcos_Freq_L2_num(i, 3);
% n = BPSK_Freq_L2_num(i, 1); n8 = n*8;
% m = 0; m8 = m*8;
% freq = BPSK_Freq_L2_num(i, 2);
% freq_index = freq - 1558 + 1;
% for f_index = freq_index
%     for n8 = n8
%         for m8 = m8 
        
for f_index = 1:fmax
    for n8 = 1:80
        for m8 = 1:80

            if m8 < n8
                if Signal_Type ~= BPSK
                    continue;
                end
            end

            if Signal_Type == BOCsin 
                if (~isnan(InterSysJam_BoCsin_L2_BoC_10_5(m8, n8, f_index))) && ...
                        (~isnan(InterSysJam_BoCsin_L2_BoC_0_1(m8, n8, f_index))) && ...
                        (~isnan(InterSysJam_BoCsin_L2_BoC_0_10(m8, n8, f_index)))
                    continue;
                end
            elseif Signal_Type == BOCcos
                if (~isnan(InterSysJam_BoCcos_L2_BoC_10_5(m8, n8, f_index))) && ...
                        (~isnan(InterSysJam_BoCcos_L2_BoC_0_1(m8, n8, f_index))) && ...
                        (~isnan(InterSysJam_BoCcos_L2_BoC_0_10(m8, n8, f_index)))
                    continue;
                end
            elseif Signal_Type == BPSK
                if (~isnan(InterSysJam_BPSK_L2_BoC_10_5(n8, f_index))) && ...
                        (~isnan(InterSysJam_BPSK_L2_BoC_0_1(n8, f_index))) && ...
                        (~isnan(InterSysJam_BPSK_L2_BoC_0_10(n8, f_index)))
                    continue;
                end                
            end

            if ((m8+n8)/8 > f_index + 2) || ((m8+n8)/8 > (18-f_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end

            m = m8/8 * (Signal_Type ~= BPSK);
            n = n8/8;

            % Открываем файл АКФ указанного сигнала
            ro_our = get_ro(m, n, Signal_Type, path_to_ro);
            if ro_our == 0
                continue
            end
            N_ro = length(ro_our);
            N_ro_dop = N_ro*4 + 1;
            ro_our_dop = zeros(1, N_ro_dop);
            ro_our_dop(1:N_ro) = ro_our;

            N_ro_1 = length(ro_GPS_BoC_0_1); % Число точек при n >= 1
            offset = (N_ro_1 - 1) / 2;
            offset = (N_ro - 1) / 2 - offset;
            
            % Фильтрация 
            if n_plus_m_old ~= ceil(n+m)
                load([pwd '/filters/Hd_Band_' sprintf('%.0f', ceil(n+m)) '.mat']);        
                n_plus_m_old = ceil(n+m);
            end
            ro_our_f = filter(Hd, ro_our_dop);
        
            if (N_ro_old ~= N_ro) % Если вдруг поменялось
                % Дополнение нулями
                ro_GPS_BoC_10_5_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_10_5_dop( (1:length(ro_GPS_BoC_10_5)) + offset ) = ro_GPS_BoC_10_5;
                ro_GPS_BoC_0_1_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_0_1_dop( (1:length(ro_GPS_BoC_0_1)) + offset ) = ro_GPS_BoC_0_1;
                ro_GPS_BoC_0_10_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_0_10_dop( (1:length(ro_GPS_BoC_0_10)) + offset ) = ro_GPS_BoC_0_10;
            end
        
            if (Hd_old ~= Hd)||(N_ro_old ~= N_ro)||(f_index_old ~= f_index)
                % Конвертер частоты
                f_sig = 1.023e6*farr(f_index); % Частота нашего сигнала        
                dFreq = f_intro - f_sig;
                cos_df = cos(2*pi*dFreq*( (1:N_ro_1) - ((N_ro_1 - 1)/2 + 1) )*Td); % С исходным темпом
                cos_df_dop = exp(1i*2*pi*dFreq*( (1:N_ro_dop) - ((N_ro - 1)/2 + 1) )*Td); % Дополненный справа

                % Да, при таком подходе расчет долог. Но зато прозрачно и кода
                % меньше. Если бы не менялось N_ro, можно было бы один раз это
                % посчитать и использовать
                ro_GPS_BoC_10_5_dop_f = filter(Hd, cos_df_dop.*ro_GPS_BoC_10_5_dop);
                ro_GPS_BoC_0_1_dop_f = filter(Hd, cos_df_dop.*ro_GPS_BoC_0_1_dop);
                ro_GPS_BoC_0_10_dop_f = filter(Hd, cos_df_dop.*ro_GPS_BoC_0_10_dop);
                    
            end
            Hd_old = Hd;
            N_ro_old = N_ro;
            f_index_old = f_index;

            % Расчет интегралов - коэффициентов спектрального разделения
            Inte_GPS_BoC_10_5  = 10*log10(abs(ro_our_f*(ro_GPS_BoC_10_5_dop_f)'*Td));
            Inte_GPS_BoC_0_1  = 10*log10(abs(ro_our_f*(ro_GPS_BoC_0_1_dop_f)'*Td));
            Inte_GPS_BoC_0_10  = 10*log10(abs(ro_our_f*(ro_GPS_BoC_0_10_dop_f)'*Td));
            
            if Signal_Type == BOCsin 
                InterSysJam_BoCsin_L2_BoC_10_5(m8, n8, f_index) = Inte_GPS_BoC_10_5;
                InterSysJam_BoCsin_L2_BoC_0_1(m8, n8, f_index) = Inte_GPS_BoC_0_1;
                InterSysJam_BoCsin_L2_BoC_0_10(m8, n8, f_index) = Inte_GPS_BoC_0_10;
            elseif Signal_Type == BOCcos
                InterSysJam_BoCcos_L2_BoC_10_5(m8, n8, f_index) = Inte_GPS_BoC_10_5;
                InterSysJam_BoCcos_L2_BoC_0_1(m8, n8, f_index) = Inte_GPS_BoC_0_1;
                InterSysJam_BoCcos_L2_BoC_0_10(m8, n8, f_index) = Inte_GPS_BoC_0_10;
            elseif Signal_Type == BPSK
                InterSysJam_BPSK_L2_BoC_10_5(n8, f_index) = Inte_GPS_BoC_10_5;
                InterSysJam_BPSK_L2_BoC_0_1(n8, f_index) = Inte_GPS_BoC_0_1;
                InterSysJam_BPSK_L2_BoC_0_10(n8, f_index) = Inte_GPS_BoC_0_10;
            end
            
            fprintf('Intersystem Jamm BoC(%.3f, %.3f) at %.0f \n \t with GPS BoC(10,5) = %.2f dB\n \t with GPS BPSK(1) = %.2f dB\n \t with GPS BPSK(10) = %.2f dB\n', m, n, farr(f_index), ...
                Inte_GPS_BoC_10_5, Inte_GPS_BoC_0_1, Inte_GPS_BoC_0_10);
            
            % Для BPSK по m пробегать не надо
            if (Signal_Type == BPSK)
                break;
            end
        end
        
        if Signal_Type == BOCsin 
            save([path_to_results '/InterSysJam_BoCsin_L2_BoC_10_5.mat'], 'InterSysJam_BoCsin_L2_BoC_10_5');
            save([path_to_results '/InterSysJam_BoCsin_L2_BoC_0_1.mat'], 'InterSysJam_BoCsin_L2_BoC_0_1');
            save([path_to_results '/InterSysJam_BoCsin_L2_BoC_0_10.mat'], 'InterSysJam_BoCsin_L2_BoC_0_10');
        elseif Signal_Type == BOCcos
            save([path_to_results '/InterSysJam_BoCcos_L2_BoC_10_5.mat'], 'InterSysJam_BoCcos_L2_BoC_10_5');
            save([path_to_results '/InterSysJam_BoCcos_L2_BoC_0_1.mat'], 'InterSysJam_BoCcos_L2_BoC_0_1');
            save([path_to_results '/InterSysJam_BoCcos_L2_BoC_0_10.mat'], 'InterSysJam_BoCcos_L2_BoC_0_10');    
        elseif Signal_Type == BPSK
            save([path_to_results '/InterSysJam_BPSK_L2_BoC_10_5.mat'], 'InterSysJam_BPSK_L2_BoC_10_5');
            save([path_to_results '/InterSysJam_BPSK_L2_BoC_0_1.mat'], 'InterSysJam_BPSK_L2_BoC_0_1');
            save([path_to_results '/InterSysJam_BPSK_L2_BoC_0_10.mat'], 'InterSysJam_BPSK_L2_BoC_0_10');               
        end        
    end
end
% end
hF = 0;


for i = 1:fmax
    hF = figure(hF+1);
    if (Signal_Type == BOCsin)
        pcolor((1:80)/8, (1:80)/8, InterSysJam_BoCsin_L2_BoC_0_1(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (Signal_Type == BOCcos)
        pcolor((1:80)/8, (1:80)/8, InterSysJam_BoCcos_L2_BoC_0_1(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (Signal_Type == BPSK)
        plot((1:80)/8, InterSysJam_BPSK_L2_BoC_0_1(1:80, i));        
        xlabel('n')
        ylabel('k_cd')
    end
    title(sprintf('Normalized freq = %.0f', farr(i)));
end

hF = figure(hF + 1);
plot(1:N_ro, ro_our, 1:N_ro_1, ro_GPS_BoC_10_5, 1:N_ro_1, ro_GPS_BoC_10_5.*cos_df, ...
     1:N_ro_dop, ro_our_f, 1:N_ro_dop, ro_GPS_BoC_10_5_dop_f);
title('All')

ff = (-(N_ro/2 - 1):1:N_ro/2)/(Td*1e6)/N_ro; % Ось частот для недополненной АКФ
ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td*1e6)/N_ro_dop;  % Ось частот для дополненной АКФ

hF = figure(hF + 1);
plot(ff_dop, (abs(fftshift(fft(ro_our_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_GPS_BoC_10_5_dop)))), ...
     ff_dop, abs((fftshift(fft(ro_GPS_BoC_10_5_dop.*cos_df_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_our_f)))), ...
     ff_dop, (abs(fftshift(fft(ro_GPS_BoC_10_5_dop_f))))   )
xlabel('MHz')