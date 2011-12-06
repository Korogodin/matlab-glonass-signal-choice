clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/back_intersystem_L3'];

n8max = 80;
m8max = 80;
farr = 1164:1184; fmax = length(farr); % Нормированный центральные частоты
Signal_Band = 20;

% BackInterSysJam_BoCsin_L3_BoC_0_10 = nan(m8max, n8max, fmax);
% BackInterSysJam_BoCcos_L3_BoC_0_10 = nan(m8max, n8max, fmax);
% BackInterSysJam_BPSK_L3_BoC_0_10 = nan(n8max, fmax);
% save([path_to_results '/BackInterSysJam_BoCsin_L3_BoC_0_10.mat'], 'BackInterSysJam_BoCsin_L3_BoC_0_10');
% save([path_to_results '/BackInterSysJam_BoCcos_L3_BoC_0_10.mat'], 'BackInterSysJam_BoCcos_L3_BoC_0_10');
% save([path_to_results '/BackInterSysJam_BPSK_L3_BoC_0_10.mat'], 'BackInterSysJam_BPSK_L3_BoC_0_10');
% return

load([path_to_results '/BackInterSysJam_BoCsin_L3_BoC_0_10.mat']);
load([path_to_results '/BackInterSysJam_BoCcos_L3_BoC_0_10.mat']);
load([path_to_results '/BackInterSysJam_BPSK_L3_BoC_0_10.mat']);

% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

load([pwd '/ro/Td.mat']);

% АКФ сигнала GPS BPSK(10)
ro_GPS_BoC_0_10 = get_ro(0, 10, BPSK, path_to_ro );

f_intro = 1176.45e6; % Частота рассматриваемых сигналов (теперь полезных)

load([pwd '/filters/Hd_for_BoC_0_10.mat']);
Hd_0_10 = Hd;

N_ro_old = 0;
f_index_old = 0;

for f_index = 1:fmax
    for n8 = 1:80
        for m8 = 1:80  
         
            if m8 < n8
                if Signal_Type ~= BPSK
                    continue;
                end
            end
            
            m = m8/8 * (Signal_Type ~= BPSK);
            n = n8/8;

            
            if ((m+n) > f_index + 2) || ((m+n) > (Signal_Band-f_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            % Если уже посчитано, то идем дальше
            if Signal_Type == BOCsin 
                if (~isnan(BackInterSysJam_BoCsin_L3_BoC_0_10(m8, n8, f_index))) 
                    continue;
                end
            elseif Signal_Type == BOCcos
                if (~isnan(BackInterSysJam_BoCcos_L3_BoC_0_10(m8, n8, f_index))) 
                    continue;
                end
            elseif Signal_Type == BPSK
                if (~isnan(BackInterSysJam_BPSK_L3_BoC_0_10(n8, f_index)))
                    continue;
                end                
            end

            % Открываем файл АКФ указанного сигнала
            ro_our = get_ro(m, n, Signal_Type, path_to_ro);
            if ro_our == 0
                continue
            end
            N_ro = length(ro_our);
            N_ro_dop = N_ro*4 + 1;
            ro_our_dop = zeros(1, N_ro_dop);
            ro_our_dop(1:N_ro) = ro_our;
            
            N_ro_1 = length(ro_GPS_BoC_0_10); % Число точек при n >= 1
            offset = (N_ro_1 - 1) / 2;
            offset = (N_ro - 1) / 2 - offset;
            if (N_ro_old ~= N_ro) || (f_index_old ~= f_index) % Если вдруг поменялось
                % Дополнение нулями
      
                ro_GPS_BoC_0_10_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_0_10_dop( (1:length(ro_GPS_BoC_0_10)) + offset ) = ro_GPS_BoC_0_10;

                % Конвертер частоты
                f_sig = 1.023e6*farr(f_index); % Частота нашего сигнала        
                dFreq = f_intro - f_sig;
                cos_df = cos(2*pi*dFreq*( (1:N_ro_1) - ((N_ro_1 - 1)/2 + 1) )*Td); % С исходным темпом
                cos_df_dop = exp(1i*2*pi*dFreq*( (1:N_ro_dop) - ((N_ro - 1)/2 + 1) )*Td); % Дополненный справа

                % Фильтрация (теперь для каждого полезного сигнала свой фильтр)
           
                ro_GPS_BoC_0_10_dop_f = filter(Hd_0_10, ro_GPS_BoC_0_10_dop);                
            end
            
            ro_our_f_0_10 = filter(Hd_0_10, cos_df_dop.*ro_our_dop);
            
            N_ro_old = N_ro;
            f_index_old = f_index;

            % Расчет интегралов - коэффициентов спектрального разделения
            Inte_GPS_BoC_0_10  = 10*log10(abs(ro_our_f_0_10*(ro_GPS_BoC_0_10_dop_f)'*Td));
            
            if Signal_Type == BOCsin 
                BackInterSysJam_BoCsin_L3_BoC_0_10(m8, n8, f_index) = Inte_GPS_BoC_0_10;
            elseif Signal_Type == BOCcos
                BackInterSysJam_BoCcos_L3_BoC_0_10(m8, n8, f_index) = Inte_GPS_BoC_0_10;
            elseif Signal_Type == BPSK
                BackInterSysJam_BPSK_L3_BoC_0_10(n8, f_index) = Inte_GPS_BoC_0_10;
            end
            
            fprintf('BackIntersystem Jamm BoC(%.3f, %.3f) at %.0f with GPS BPSK(10) = %.2f dB\n', m, n, farr(f_index), ...
                Inte_GPS_BoC_0_10);
            
            % Для BPSK по m пробегать не надо
            if (Signal_Type == BPSK)
                break;
            end
        end
        if Signal_Type == BOCsin 
            save([path_to_results '/BackInterSysJam_BoCsin_L3_BoC_0_10.mat'], 'BackInterSysJam_BoCsin_L3_BoC_0_10');
        elseif Signal_Type == BOCcos
            save([path_to_results '/BackInterSysJam_BoCcos_L3_BoC_0_10.mat'], 'BackInterSysJam_BoCcos_L3_BoC_0_10');    
        elseif Signal_Type == BPSK
            save([path_to_results '/BackInterSysJam_BPSK_L3_BoC_0_10.mat'], 'BackInterSysJam_BPSK_L3_BoC_0_10');               
        end
    end
end
% end

hF = 0;

% Сигнал BPSK(10) есть во всех частотных диапазонах, удобно для унификации
for i = 1:fmax
    hF = figure(hF+1);
    if (Signal_Type == BOCsin)
        pcolor((1:80)/8, (1:80)/8, BackInterSysJam_BoCsin_L3_BoC_0_10(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (Signal_Type == BOCcos)
        pcolor((1:80)/8, (1:80)/8, BackInterSysJam_BoCcos_L3_BoC_0_10(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (Signal_Type == BPSK)
        plot((1:80)/8, BackInterSysJam_BPSK_L3_BoC_0_10(1:80, i));        
        xlabel('n')
        ylabel('k_cd')
    end
    title(sprintf('Normalized freq = %.0f', farr(i)));
end

hF = figure(hF + 1);
plot(1:N_ro, ro_our, 1:N_ro_1, ro_GPS_BoC_0_10, 1:N_ro_dop, cos_df_dop.*ro_our_dop, ...
     1:N_ro_dop, ro_GPS_BoC_0_10_dop_f, 1:N_ro_dop, ro_our_f_0_10);
title('All')

ff = (-(N_ro/2 - 1):1:N_ro/2)/(Td*1e6)/N_ro; % Ось частот для недополненной АКФ
ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td*1e6)/N_ro_dop;  % Ось частот для дополненной АКФ

hF = figure(hF + 1);
plot(ff_dop, (abs(fftshift(fft(ro_our_dop.*cos_df_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_GPS_BoC_0_10_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_our_f_0_10)))), ...
     ff_dop, (abs(fftshift(fft(ro_GPS_BoC_0_10_dop_f))))   )
xlabel('MHz')