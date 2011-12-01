%/**
% Данный скрипт рассчитывает коэффициенты спектрального разделения между
% сигналами полного выборочного множества и сигналами ГЛОНАСС FDMA при том 
% условии, что фронтенд настроен на прием и тех, и других сигналов (полоса
% фронтенда равна полосе совокупности сигнала по нулям в теории; на 
% практике в данном скрипте фильтр просто отсутствует). Результат 
% вычисления коэффициентов частотного разделения сохраняется в массивах 
% InterSysJam_*** с индексами (8*m, 8*n, freq_ind, lit) для BoC-сигналов 
% или (8*n, freq_ind, lit) для BPSK-сигналов, которые сохраняются в 
% соответствующих mat-файлах. Есть ещё файлы с суффиксом _mean, в которых
% хранится усредненное по литерам значение. 
%@param Signal_Type задает BoCsin, BoCcos или BPSK
%@param m8max, n8max - максимальное значение индексов 8*m, 8*n. Пробегаем
%по n и m  с шагом 0.125
%@param farr - массив номированных несущих частот
%*/

clear 
close all
clc

n8max = 80;
m8max = 80;
farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
lit_arr = -7:6; lit_size = length(lit_arr);

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/intersystem_L3'];

% Empty file creation
% InterSysJam_BoCsin_L3_GloST = nan(m8max, n8max, fmax, lit_size);
% InterSysJam_BoCcos_L3_GloST = nan(m8max, n8max, fmax, lit_size);
% save([path_to_results '/InterSysJam_BoCsin_L3_GloST.mat'], 'InterSysJam_BoCsin_L3_GloST');
% save([path_to_results '/InterSysJam_BoCcos_L3_GloST.mat'], 'InterSysJam_BoCcos_L3_GloST');
% InterSysJam_BoCsin_L3_GloVT = nan(m8max, n8max, fmax, lit_size);
% InterSysJam_BoCcos_L3_GloVT = nan(m8max, n8max, fmax, lit_size);
% save([path_to_results '/InterSysJam_BoCsin_L3_GloVT.mat'], 'InterSysJam_BoCsin_L3_GloVT');
% save([path_to_results '/InterSysJam_BoCcos_L3_GloVT.mat'], 'InterSysJam_BoCcos_L3_GloVT');
% InterSysJam_BPSK_L3_GloST = nan(n8max, fmax, lit_size);
% InterSysJam_BPSK_L3_GloVT = nan(n8max, fmax, lit_size);
% save([path_to_results '/InterSysJam_BPSK_L3_GloST.mat'], 'InterSysJam_BPSK_L3_GloST');
% save([path_to_results '/InterSysJam_BPSK_L3_GloVT.mat'], 'InterSysJam_BPSK_L3_GloVT');
% 
% InterSysJam_BoCsin_L3_GloST_mean = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L3_GloST_mean = nan(m8max, n8max, fmax);
% save([path_to_results '/InterSysJam_BoCsin_L3_GloST_mean.mat'], 'InterSysJam_BoCsin_L3_GloST_mean');
% save([path_to_results '/InterSysJam_BoCcos_L3_GloST_mean.mat'], 'InterSysJam_BoCcos_L3_GloST_mean');
% InterSysJam_BoCsin_L3_GloVT_mean = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L3_GloVT_mean = nan(m8max, n8max, fmax);
% save([path_to_results '/InterSysJam_BoCsin_L3_GloVT_mean.mat'], 'InterSysJam_BoCsin_L3_GloVT_mean');
% save([path_to_results '/InterSysJam_BoCcos_L3_GloVT_mean.mat'], 'InterSysJam_BoCcos_L3_GloVT_mean');
% InterSysJam_BPSK_L3_GloST_mean = nan(n8max, fmax);
% InterSysJam_BPSK_L3_GloVT_mean = nan(n8max, fmax);
% save([path_to_results '/InterSysJam_BPSK_L3_GloST_mean.mat'], 'InterSysJam_BPSK_L3_GloST_mean');
% save([path_to_results '/InterSysJam_BPSK_L3_GloVT_mean.mat'], 'InterSysJam_BPSK_L3_GloVT_mean');


% Load files
load([path_to_results '/InterSysJam_BoCsin_L3_GloST.mat'], 'InterSysJam_BoCsin_L3_GloST');
load([path_to_results '/InterSysJam_BoCcos_L3_GloST.mat'], 'InterSysJam_BoCcos_L3_GloST');
load([path_to_results '/InterSysJam_BoCsin_L3_GloVT.mat'], 'InterSysJam_BoCsin_L3_GloVT');
load([path_to_results '/InterSysJam_BoCcos_L3_GloVT.mat'], 'InterSysJam_BoCcos_L3_GloVT');
load([path_to_results '/InterSysJam_BPSK_L3_GloST.mat'], 'InterSysJam_BPSK_L3_GloST');
load([path_to_results '/InterSysJam_BPSK_L3_GloVT.mat'], 'InterSysJam_BPSK_L3_GloVT');
load([path_to_results '/InterSysJam_BoCsin_L3_GloST_mean.mat'], 'InterSysJam_BoCsin_L3_GloST_mean');
load([path_to_results '/InterSysJam_BoCcos_L3_GloST_mean.mat'], 'InterSysJam_BoCcos_L3_GloST_mean');
load([path_to_results '/InterSysJam_BoCsin_L3_GloVT_mean.mat'], 'InterSysJam_BoCsin_L3_GloVT_mean');
load([path_to_results '/InterSysJam_BoCcos_L3_GloVT_mean.mat'], 'InterSysJam_BoCcos_L3_GloVT_mean');
load([path_to_results '/InterSysJam_BPSK_L3_GloST_mean.mat'], 'InterSysJam_BPSK_L3_GloST_mean');
load([path_to_results '/InterSysJam_BPSK_L3_GloVT_mean.mat'], 'InterSysJam_BPSK_L3_GloVT_mean');


% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 2; % 1 - BOCsin; 2 - BOCcos; 3 - BPSK.

load([pwd '/ro/Td.mat']);

% АКФ сигнала GLO BPSK(0.5)
ro_GLO_ST = get_ro(0, 0.5, BPSK, path_to_ro);

% АКФ сигнала GLO BPSK(5)
ro_GLO_VT = get_ro(0, 5, BPSK, path_to_ro);

N_ro_1 = length(ro_GLO_VT); % Число точек при n >= 1
N_ro_05 = length(ro_GLO_ST); % Число точек при n = 0.5
offset_x1 = (N_ro_1 - 1) / 2;
offset_x8 = offset_x1*8;
offset_x2 = offset_x1*2;

offset_for_05_old = 1e11;
offset_for_1_old = 1e11;
N_ro_dop_old = 0;
for f_index = 1:fmax
    for n8 = 1:80
        for m8 = 1:80     

            do_out = 1;
            
            Inte_Glo_ST_nodB  = zeros(1, lit_size);
            Inte_Glo_VT_nodB  = zeros(1, lit_size);

            if m8 < n8
                if Signal_Type ~= BPSK
                    do_out = 0;
                    continue;
                end
            end
            
            for lit_index = 1:lit_size
                % Если уже посчитано, то идем дальше
                if Signal_Type == BOCsin 
                    if (~isnan(InterSysJam_BoCsin_L3_GloST(m8, n8, f_index, lit_index))) && ...
                            (~isnan(InterSysJam_BoCsin_L3_GloVT(m8, n8, f_index, lit_index)))
                        do_out = 0;
                        break;
                    end
                elseif Signal_Type == BOCcos
                    if (~isnan(InterSysJam_BoCcos_L3_GloST(m8, n8, f_index, lit_index))) && ...
                            (~isnan(InterSysJam_BoCcos_L3_GloVT(m8, n8, f_index, lit_index)))
                        do_out = 0;
                        break;
                    end
                elseif Signal_Type == BPSK
                    if (~isnan(InterSysJam_BPSK_L3_GloST(n8, f_index, lit_index))) && ...
                            (~isnan(InterSysJam_BPSK_L3_GloVT(n8, f_index, lit_index)))
                        do_out = 0;
                        break;
                    end                
                end

                m = m8/8 * (Signal_Type ~= BPSK);
                n = n8/8;

                % Открываем файл АКФ указанного сигнала
                if lit_index == 1
                    ro_our = get_ro(m, n, Signal_Type, path_to_ro);
                    if ro_our == 0
                        do_out = 0;
                        break; % Если файла нет
                    end
                    N_ro = length(ro_our);
                    offset_for_05 = (N_ro - 1) / 2 - offset_x2;
                    offset_for_1 = (N_ro - 1) / 2 - offset_x1;
                    N_ro_dop = max([N_ro 2*N_ro_1]) *4 + 1;
                    ro_our_dop = zeros(1, N_ro_dop);
                    if offset_for_05 <= 0
                        ro_our_dop((1:N_ro) - offset_for_05) = ro_our;
                        offset_for_1 = offset_for_1 - offset_for_05;
                        offset_for_05 = 0;
                    else
                        ro_our_dop((1:N_ro)) = ro_our;
                    end                    

                    if (N_ro_dop_old ~= N_ro_dop)||...
                            (offset_for_05_old ~= offset_for_05)||...
                            (offset_for_1_old ~= offset_for_1_old)% Если вдруг что поменялось

                        % Дополнение нулями и совмещаем максимумы
                        ro_GLO_ST_dop = zeros(1, N_ro_dop);
                            ro_GLO_ST_dop( (1:N_ro_05) + offset_for_05 ) = ro_GLO_ST;
                        ro_GLO_VT_dop = zeros(1, N_ro_dop);
                            ro_GLO_VT_dop( (1:N_ro_1) + offset_for_1 ) = ro_GLO_VT;

                    end
                    N_ro_dop_old = N_ro_dop;
                    offset_for_05_old = offset_for_05;
                    offset_for_1_old = offset_for_1;
                end
                
                % Конвертер частоты
                f_sig = 1.023e6*farr(f_index); % Частота нашего сигнала     
                f_intro = 1246e6 + lit_arr(lit_index)*0.4375e6;
                dFreq = f_intro - f_sig;
%                 cos_df_dop = cos(2*pi*dFreq*( (1:N_ro_dop) - ((N_ro - 1)/2 + 1) )*Td); % Дополненный справа
                cos_df_dop = exp(1i*2*pi*dFreq*( (1:N_ro_dop) - ((max([N_ro N_ro_05]) - 1)/2 + 1) )*Td); % Дополненный справа
                
                ro_GLO_ST_dop_f = cos_df_dop.*ro_GLO_ST_dop;
                ro_GLO_VT_dop_f = cos_df_dop.*ro_GLO_VT_dop;                
                ro_our_f =  ro_our_dop; 

                % Расчет интегралов - коэффициентов спектрального разделения
        %         Inte = ro_our*(ro_GPS_BoC_1_1.*cos_df)'*Td;
                Inte_Glo_ST_nodB(lit_index)  = abs(ro_our_f*(ro_GLO_ST_dop_f)'*Td);
                Inte_Glo_VT_nodB(lit_index)  = abs(ro_our_f*(ro_GLO_VT_dop_f)'*Td);
                Inte_Glo_ST  = 10*log10(Inte_Glo_ST_nodB(lit_index));
                Inte_Glo_VT  = 10*log10(Inte_Glo_VT_nodB(lit_index));

                if Signal_Type == BOCsin 
                    InterSysJam_BoCsin_L3_GloST(m8, n8, f_index, lit_index) = Inte_Glo_ST;
                    InterSysJam_BoCsin_L3_GloVT(m8, n8, f_index, lit_index) = Inte_Glo_VT;
                elseif Signal_Type == BOCcos
                    InterSysJam_BoCcos_L3_GloST(m8, n8, f_index, lit_index) = Inte_Glo_ST;
                    InterSysJam_BoCcos_L3_GloVT(m8, n8, f_index, lit_index) = Inte_Glo_VT;
                elseif Signal_Type == BPSK
                    InterSysJam_BPSK_L3_GloST(n8, f_index, lit_index) = Inte_Glo_ST;
                    InterSysJam_BPSK_L3_GloVT(n8, f_index, lit_index) = Inte_Glo_VT;
                end

            end

            if do_out 
                Inte_Glo_ST_mean = 10*log10(mean(Inte_Glo_ST_nodB));
                Inte_Glo_VT_mean = 10*log10(mean(Inte_Glo_VT_nodB));
                fprintf('Intersystem Jamm BoC(%.3f, %.3f) at %.0f \n \t with GloST = %.2f dB\n \t with GloVT = %.2f dB\n', m, n, farr(f_index), ...
                    Inte_Glo_ST_mean, Inte_Glo_VT_mean);

                if Signal_Type == BOCsin 
                    InterSysJam_BoCsin_L3_GloST_mean(m8, n8, f_index) = Inte_Glo_ST_mean;
                    InterSysJam_BoCsin_L3_GloVT_mean(m8, n8, f_index) = Inte_Glo_VT_mean;
                elseif Signal_Type == BOCcos
                    InterSysJam_BoCcos_L3_GloST_mean(m8, n8, f_index) = Inte_Glo_ST_mean;
                    InterSysJam_BoCcos_L3_GloVT_mean(m8, n8, f_index) = Inte_Glo_VT_mean;
                elseif Signal_Type == BPSK
                    InterSysJam_BPSK_L3_GloST_mean(n8, f_index) = Inte_Glo_ST_mean;
                    InterSysJam_BPSK_L3_GloVT_mean(n8, f_index) = Inte_Glo_VT_mean;
                end        
            end
            % Для BPSK по m пробегать не надо
            if (Signal_Type == BPSK)
                break;
            end
        end
        if Signal_Type == BOCsin 
            save([path_to_results '/InterSysJam_BoCsin_L3_GloST.mat'], 'InterSysJam_BoCsin_L3_GloST');
            save([path_to_results '/InterSysJam_BoCsin_L3_GloVT.mat'], 'InterSysJam_BoCsin_L3_GloVT');
            save([path_to_results '/InterSysJam_BoCsin_L3_GloST_mean.mat'], 'InterSysJam_BoCsin_L3_GloST_mean');
            save([path_to_results '/InterSysJam_BoCsin_L3_GloVT_mean.mat'], 'InterSysJam_BoCsin_L3_GloVT_mean');
        elseif Signal_Type == BOCcos
            save([path_to_results '/InterSysJam_BoCcos_L3_GloST.mat'], 'InterSysJam_BoCcos_L3_GloST');
            save([path_to_results '/InterSysJam_BoCcos_L3_GloVT.mat'], 'InterSysJam_BoCcos_L3_GloVT');
            save([path_to_results '/InterSysJam_BoCcos_L3_GloST_mean.mat'], 'InterSysJam_BoCcos_L3_GloST_mean');
            save([path_to_results '/InterSysJam_BoCcos_L3_GloVT_mean.mat'], 'InterSysJam_BoCcos_L3_GloVT_mean');
        elseif Signal_Type == BPSK
            save([path_to_results '/InterSysJam_BPSK_L3_GloST.mat'], 'InterSysJam_BPSK_L3_GloST');
            save([path_to_results '/InterSysJam_BPSK_L3_GloVT.mat'], 'InterSysJam_BPSK_L3_GloVT');
            save([path_to_results '/InterSysJam_BPSK_L3_GloST_mean.mat'], 'InterSysJam_BPSK_L3_GloST_mean');
            save([path_to_results '/InterSysJam_BPSK_L3_GloVT_mean.mat'], 'InterSysJam_BPSK_L3_GloVT_mean');           
        end                
    end
end
hF = 0;


for i = 1:fmax
    hF = figure(hF+1);
    if (Signal_Type == BOCsin)
        pcolor((1:80)/8, (1:80)/8, InterSysJam_BoCsin_L3_GloST_mean(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (Signal_Type == BOCcos)
        pcolor((1:80)/8, (1:80)/8, InterSysJam_BoCcos_L3_GloST_mean(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (Signal_Type == BPSK)
        plot((1:80)/8, InterSysJam_BPSK_L3_GloST_mean(1:80, i));        
        xlabel('n')
        ylabel('k_cd')
    end
    title(sprintf('Normalized freq = %.0f', farr(i)));
end

hF = figure(hF + 1);
plot(1:N_ro, ro_our, 1:N_ro_1, ro_GLO_VT, 1:N_ro_dop, ro_GLO_VT_dop, ...
     1:N_ro_dop, ro_our_f, 1:N_ro_dop, ro_GLO_VT_dop_f);
title('All')

ff = (-(N_ro/2 - 1):1:N_ro/2)/(Td*1e6)/N_ro; % Ось частот для недополненной АКФ
ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td*1e6)/N_ro_dop;  % Ось частот для дополненной АКФ

hF = figure(hF + 1);
plot(ff_dop, (abs(fftshift(fft(ro_our_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_GLO_VT_dop)))), ...
     ff_dop, abs((fftshift(fft(ro_GLO_VT_dop.*cos_df_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_our_f)))), ...
     ff_dop, (abs(fftshift(fft(ro_GLO_VT_dop_f))))   )
xlabel('MHz')