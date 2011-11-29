clear 
close all
clc

n8max = 80;
m8max = 80;
farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты

BoCsincos_L1; % Параметры приемлимых сигналов

% BackInterSysJam_BoCsin_L1_BoC_1_1 = nan(m8max, n8max, fmax);
% BackInterSysJam_BoCcos_L1_BoC_1_1 = nan(m8max, n8max, fmax);
% BackInterSysJam_BPSK_L1_BoC_1_1 = nan(n8max, fmax);
% save('results/back/BackInterSysJam_BoCsin_L1_BoC_1_1.mat', 'BackInterSysJam_BoCsin_L1_BoC_1_1');
% save('results/back/BackInterSysJam_BoCcos_L1_BoC_1_1.mat', 'BackInterSysJam_BoCcos_L1_BoC_1_1');
% save('results/back/BackInterSysJam_BPSK_L1_BoC_1_1.mat', 'BackInterSysJam_BPSK_L1_BoC_1_1');
% 
% BackInterSysJam_BoCsin_L1_BoC_6_1 = nan(m8max, n8max, fmax);
% BackInterSysJam_BoCcos_L1_BoC_6_1 = nan(m8max, n8max, fmax);
% BackInterSysJam_BPSK_L1_BoC_6_1 = nan(n8max, fmax);
% save('results/back/BackInterSysJam_BoCsin_L1_BoC_6_1.mat', 'BackInterSysJam_BoCsin_L1_BoC_6_1');
% save('results/back/BackInterSysJam_BoCcos_L1_BoC_6_1.mat', 'BackInterSysJam_BoCcos_L1_BoC_6_1');
% save('results/back/BackInterSysJam_BPSK_L1_BoC_6_1.mat', 'BackInterSysJam_BPSK_L1_BoC_6_1');
% 
% BackInterSysJam_BoCsin_L1_BoC_10_5 = nan(m8max, n8max, fmax);
% BackInterSysJam_BoCcos_L1_BoC_10_5 = nan(m8max, n8max, fmax);
% BackInterSysJam_BPSK_L1_BoC_10_5 = nan(n8max, fmax);
% save('results/back/BackInterSysJam_BoCsin_L1_BoC_10_5.mat', 'BackInterSysJam_BoCsin_L1_BoC_10_5');
% save('results/back/BackInterSysJam_BoCcos_L1_BoC_10_5.mat', 'BackInterSysJam_BoCcos_L1_BoC_10_5');
% save('results/back/BackInterSysJam_BPSK_L1_BoC_10_5.mat', 'BackInterSysJam_BPSK_L1_BoC_10_5');
% 
% BackInterSysJam_BoCsin_L1_BoC_0_1 = nan(m8max, n8max, fmax);
% BackInterSysJam_BoCcos_L1_BoC_0_1 = nan(m8max, n8max, fmax);
% BackInterSysJam_BPSK_L1_BoC_0_1 = nan(n8max, fmax);
% save('results/back/BackInterSysJam_BoCsin_L1_BoC_0_1.mat', 'BackInterSysJam_BoCsin_L1_BoC_0_1');
% save('results/back/BackInterSysJam_BoCcos_L1_BoC_0_1.mat', 'BackInterSysJam_BoCcos_L1_BoC_0_1');
% save('results/back/BackInterSysJam_BPSK_L1_BoC_0_1.mat', 'BackInterSysJam_BPSK_L1_BoC_0_1');
% 
% BackInterSysJam_BoCsin_L1_BoC_0_10 = nan(m8max, n8max, fmax);
% BackInterSysJam_BoCcos_L1_BoC_0_10 = nan(m8max, n8max, fmax);
% BackInterSysJam_BPSK_L1_BoC_0_10 = nan(n8max, fmax);
% save('results/back/BackInterSysJam_BoCsin_L1_BoC_0_10.mat', 'BackInterSysJam_BoCsin_L1_BoC_0_10');
% save('results/back/BackInterSysJam_BoCcos_L1_BoC_0_10.mat', 'BackInterSysJam_BoCcos_L1_BoC_0_10');
% save('results/back/BackInterSysJam_BPSK_L1_BoC_0_10.mat', 'BackInterSysJam_BPSK_L1_BoC_0_10');
% return

load('results/back/BackInterSysJam_BoCsin_L1_BoC_1_1.mat'); 
load('results/back/BackInterSysJam_BoCcos_L1_BoC_1_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_1_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_6_1.mat'); 
load('results/back/BackInterSysJam_BoCcos_L1_BoC_6_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_6_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_6_1.mat'); 
load('results/back/BackInterSysJam_BoCcos_L1_BoC_6_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_6_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_10_5.mat');
load('results/back/BackInterSysJam_BoCcos_L1_BoC_10_5.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_10_5.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_0_1.mat');
load('results/back/BackInterSysJam_BoCcos_L1_BoC_0_1.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_0_1.mat');
load('results/back/BackInterSysJam_BoCsin_L1_BoC_0_10.mat');
load('results/back/BackInterSysJam_BoCcos_L1_BoC_0_10.mat');
load('results/back/BackInterSysJam_BPSK_L1_BoC_0_10.mat');

% Параметры нашего сигнала
BOC_Type = 2; % 1 - sin, 2 - cos, 3 - BPSK
load([pwd '/ro/Td.mat']);

% АКФ сигнала GPS BoC(1, 1)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 1) ', ' sprintf('%.3f', 1) ').mat'])
ro_GPS_BoC_1_1 = ro; 

% АКФ сигнала GPS BoC(6, 1)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 6) ', ' sprintf('%.3f', 1) ').mat'])
ro_GPS_BoC_6_1 = ro; 

% АКФ сигнала GPS BoC(10, 5)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 10) ', ' sprintf('%.3f', 5) ').mat'])
ro_GPS_BoC_10_5 = ro; 

% АКФ сигнала GPS BPSK(1)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', 1) ').mat'])
ro_GPS_BoC_0_1 = ro; 

% АКФ сигнала GPS BPSK(10)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', 10) ').mat'])
ro_GPS_BoC_0_10 = ro; 

f_intro = 1575.42e6; % Частота рассматриваемых сигналов (теперь полезных)

load([pwd '/filters/Hd_for_BoC_6_1.mat']);
Hd_6_1 = Hd;

load([pwd '/filters/Hd_for_BoC_10_5.mat']);
Hd_10_5 = Hd;

load([pwd '/filters/Hd_for_BoC_0_1.mat']);
Hd_0_1 = Hd;

load([pwd '/filters/Hd_for_BoC_0_10.mat']);
Hd_0_10 = Hd;

N_ro_old = 0;
f_index_old = 0;

% for i = 1:size(BPSK_Freq_L1_num, 1)

% n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
% m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
% freq = BoCcos_Freq_L1_num(i, 3);
% n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
% m = 0; m8 = m*8;
% freq = BPSK_Freq_L1_num(i, 2);
% freq_index = freq - 1558 + 1;
% for f_index = freq_index
%     for n8 = n8
%         for m8 = m8 
for f_index = 1:fmax
    for n8 = 1:80
        for m8 = 1:80  
% for f_index = 8
%     for n8 = 20
%         for m8 = 40               
         
            if m8 < n8
                continue;
            end
            if ((m8+n8)/8 > f_index + 2) || ((m8+n8)/8 > (16-f_index) +2) % Если этот сигнал не влазиет в полосу
                continue;
            end
            % Если уже посчитано, то идем дальше
            if BOC_Type == 1 
%                 if n8 >= 8
                    if (~isnan(BackInterSysJam_BoCsin_L1_BoC_1_1(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCsin_L1_BoC_6_1(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, f_index)))
                        continue;
                    end
%                 end
            elseif BOC_Type == 2
%                 if n8 >= 8                
                    if (~isnan(BackInterSysJam_BoCcos_L1_BoC_1_1(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCcos_L1_BoC_6_1(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, f_index)))
                        continue;
                    end
%                 end
            elseif BOC_Type == 3
%                 if n8 >= 8
                    if (~isnan(BackInterSysJam_BPSK_L1_BoC_1_1(n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BPSK_L1_BoC_6_1(n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BPSK_L1_BoC_10_5(n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BPSK_L1_BoC_0_1(n8, f_index))) && ...
                            (~isnan(BackInterSysJam_BPSK_L1_BoC_0_10(n8, f_index)))
                        continue;
                    end                
%                 end
            end

            m = m8/8;
            n = n8/8;

            % Открываем файл АКФ указанного сигнала
            try
                if BOC_Type == 1
                    load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').mat'])
                elseif BOC_Type == 2
                    load([pwd '/ro/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').mat'])
                elseif BOC_Type == 3
                    load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', n) ').mat'])
                end
            catch exception
                continue; % Если файла нет
            end
            ro_our = ro;
            N_ro = length(ro);
            N_ro_dop = N_ro*4 + 1;
            ro_our_dop = zeros(1, N_ro_dop);
            ro_our_dop(1:N_ro) = ro_our;
            
            N_ro_1 = length(ro_GPS_BoC_1_1); % Число точек при n >= 1
            offset = (N_ro_1 - 1) / 2;
            offset = (N_ro - 1) / 2 - offset;
            if (N_ro_old ~= N_ro) || (f_index_old ~= f_index) % Если вдруг поменялось
                % Дополнение нулями
                ro_GPS_BoC_1_1_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_1_1_dop( (1:length(ro_GPS_BoC_1_1)) + offset ) = ro_GPS_BoC_1_1;
                ro_GPS_BoC_6_1_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_6_1_dop( (1:length(ro_GPS_BoC_6_1)) + offset ) = ro_GPS_BoC_6_1;                    
                ro_GPS_BoC_10_5_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_10_5_dop( (1:length(ro_GPS_BoC_10_5)) + offset ) = ro_GPS_BoC_10_5;
                ro_GPS_BoC_0_1_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_0_1_dop( (1:length(ro_GPS_BoC_0_1)) + offset ) = ro_GPS_BoC_0_1;
                ro_GPS_BoC_0_10_dop = zeros(1, N_ro_dop);
                    ro_GPS_BoC_0_10_dop( (1:length(ro_GPS_BoC_0_10)) + offset ) = ro_GPS_BoC_0_10;

                % Конвертер частоты
                f_sig = 1.023e6*farr(f_index); % Частота нашего сигнала        
                dFreq = f_intro - f_sig;
                cos_df = cos(2*pi*dFreq*( (1:N_ro_1) - ((N_ro_1 - 1)/2 + 1) )*Td); % С исходным темпом
%                 cos_df_dop = cos(2*pi*dFreq*( (1:N_ro_dop) - ((N_ro - 1)/2 + 1) )*Td); % Дополненный справа
                cos_df_dop = exp(1i*2*pi*dFreq*( (1:N_ro_dop) - ((N_ro - 1)/2 + 1) )*Td); % Дополненный справа

                % Фильтрация (теперь для каждого полезного сигнала свой фильтр)
                ro_GPS_BoC_1_1_dop_f = filter(Hd_6_1, ro_GPS_BoC_1_1_dop);
                ro_GPS_BoC_6_1_dop_f = filter(Hd_6_1, ro_GPS_BoC_6_1_dop);                
                ro_GPS_BoC_10_5_dop_f = filter(Hd_10_5, ro_GPS_BoC_10_5_dop);
                ro_GPS_BoC_0_1_dop_f = filter(Hd_0_1, ro_GPS_BoC_0_1_dop);
                ro_GPS_BoC_0_10_dop_f = filter(Hd_0_10, ro_GPS_BoC_0_10_dop);                
            end
            
            ro_our_f_6_1 = filter(Hd_6_1, cos_df_dop.*ro_our_dop);
            ro_our_f_10_5 = filter(Hd_10_5, cos_df_dop.*ro_our_dop);
            ro_our_f_0_1 = filter(Hd_0_1, cos_df_dop.*ro_our_dop);
            ro_our_f_0_10 = filter(Hd_0_10, cos_df_dop.*ro_our_dop);
            
            N_ro_old = N_ro;
            f_index_old = f_index;

            % Расчет интегралов - коэффициентов спектрального разделения
            Inte_GPS_BoC_1_1  = 10*log10(abs(ro_our_f_6_1*(ro_GPS_BoC_1_1_dop_f)'*Td));
            Inte_GPS_BoC_6_1  = 10*log10(abs(ro_our_f_6_1*(ro_GPS_BoC_6_1_dop_f)'*Td));
            Inte_GPS_BoC_10_5  = 10*log10(abs(ro_our_f_10_5*(ro_GPS_BoC_10_5_dop_f)'*Td));
            Inte_GPS_BoC_0_1  = 10*log10(abs(ro_our_f_0_1*(ro_GPS_BoC_0_1_dop_f)'*Td));
            Inte_GPS_BoC_0_10  = 10*log10(abs(ro_our_f_0_10*(ro_GPS_BoC_0_10_dop_f)'*Td));
            
            if BOC_Type == 1 
                BackInterSysJam_BoCsin_L1_BoC_1_1(m8, n8, f_index) = Inte_GPS_BoC_1_1;
                BackInterSysJam_BoCsin_L1_BoC_6_1(m8, n8, f_index) = Inte_GPS_BoC_6_1;
                BackInterSysJam_BoCsin_L1_BoC_10_5(m8, n8, f_index) = Inte_GPS_BoC_10_5;
                BackInterSysJam_BoCsin_L1_BoC_0_1(m8, n8, f_index) = Inte_GPS_BoC_0_1;
                BackInterSysJam_BoCsin_L1_BoC_0_10(m8, n8, f_index) = Inte_GPS_BoC_0_10;
            elseif BOC_Type == 2
                BackInterSysJam_BoCcos_L1_BoC_1_1(m8, n8, f_index) = Inte_GPS_BoC_1_1;                
                BackInterSysJam_BoCcos_L1_BoC_6_1(m8, n8, f_index) = Inte_GPS_BoC_6_1;
                BackInterSysJam_BoCcos_L1_BoC_10_5(m8, n8, f_index) = Inte_GPS_BoC_10_5;
                BackInterSysJam_BoCcos_L1_BoC_0_1(m8, n8, f_index) = Inte_GPS_BoC_0_1;
                BackInterSysJam_BoCcos_L1_BoC_0_10(m8, n8, f_index) = Inte_GPS_BoC_0_10;
            elseif BOC_Type == 3
                BackInterSysJam_BPSK_L1_BoC_1_1(n8, f_index) = Inte_GPS_BoC_1_1;                
                BackInterSysJam_BPSK_L1_BoC_6_1(n8, f_index) = Inte_GPS_BoC_6_1; 
                BackInterSysJam_BPSK_L1_BoC_10_5(n8, f_index) = Inte_GPS_BoC_10_5;
                BackInterSysJam_BPSK_L1_BoC_0_1(n8, f_index) = Inte_GPS_BoC_0_1;
                BackInterSysJam_BPSK_L1_BoC_0_10(n8, f_index) = Inte_GPS_BoC_0_10;
            end
            
            fprintf('BackIntersystem Jamm BoC(%.3f, %.3f) at %.0f \n \t with GPS BoC(1,1) = %.2f dB\n \t with GPS BoC(6,1) = %.2f dB\n \t with GPS BoC(10,5) = %.2f dB\n \t with GPS BPSK(1) = %.2f dB\n \t with GPS BPSK(10) = %.2f dB\n', m, n, farr(f_index), ...
                Inte_GPS_BoC_1_1, Inte_GPS_BoC_6_1, Inte_GPS_BoC_10_5, Inte_GPS_BoC_0_1, Inte_GPS_BoC_0_10);
        end
        if BOC_Type == 1 
            save('results/back/BackInterSysJam_BoCsin_L1_BoC_1_1.mat', 'BackInterSysJam_BoCsin_L1_BoC_1_1');
            save('results/back/BackInterSysJam_BoCsin_L1_BoC_6_1.mat', 'BackInterSysJam_BoCsin_L1_BoC_6_1');
            save('results/back/BackInterSysJam_BoCsin_L1_BoC_10_5.mat', 'BackInterSysJam_BoCsin_L1_BoC_10_5');
            save('results/back/BackInterSysJam_BoCsin_L1_BoC_0_1.mat', 'BackInterSysJam_BoCsin_L1_BoC_0_1');
            save('results/back/BackInterSysJam_BoCsin_L1_BoC_0_10.mat', 'BackInterSysJam_BoCsin_L1_BoC_0_10');
        elseif BOC_Type == 2
            save('results/back/BackInterSysJam_BoCcos_L1_BoC_1_1.mat', 'BackInterSysJam_BoCcos_L1_BoC_1_1');
            save('results/back/BackInterSysJam_BoCcos_L1_BoC_6_1.mat', 'BackInterSysJam_BoCcos_L1_BoC_6_1');
            save('results/back/BackInterSysJam_BoCcos_L1_BoC_10_5.mat', 'BackInterSysJam_BoCcos_L1_BoC_10_5');
            save('results/back/BackInterSysJam_BoCcos_L1_BoC_0_1.mat', 'BackInterSysJam_BoCcos_L1_BoC_0_1');
            save('results/back/BackInterSysJam_BoCcos_L1_BoC_0_10.mat', 'BackInterSysJam_BoCcos_L1_BoC_0_10');    
        elseif BOC_Type == 3
            save('results/back/BackInterSysJam_BPSK_L1_BoC_1_1.mat', 'BackInterSysJam_BPSK_L1_BoC_1_1');
            save('results/back/BackInterSysJam_BPSK_L1_BoC_6_1.mat', 'BackInterSysJam_BPSK_L1_BoC_6_1');
            save('results/back/BackInterSysJam_BPSK_L1_BoC_10_5.mat', 'BackInterSysJam_BPSK_L1_BoC_10_5');
            save('results/back/BackInterSysJam_BPSK_L1_BoC_0_1.mat', 'BackInterSysJam_BPSK_L1_BoC_0_1');
            save('results/back/BackInterSysJam_BPSK_L1_BoC_0_10.mat', 'BackInterSysJam_BPSK_L1_BoC_0_10');               
        end        
    end
end
% end

hF = 0;


for i = 1:fmax
    hF = figure(hF+1);
    if (BOC_Type == 1)
        pcolor((1:80)/8, (1:80)/8, BackInterSysJam_BoCsin_L1_BoC_1_1(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (BOC_Type == 2)
        pcolor((1:80)/8, (1:80)/8, BackInterSysJam_BoCcos_L1_BoC_1_1(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (BOC_Type == 3)
        plot((1:80)/8, BackInterSysJam_BPSK_L1_BoC_1_1(1:80, i));        
        xlabel('n')
        ylabel('k_cd')
    end
    title(sprintf('Normalized freq = %.0f', farr(i)));
end


hF = figure(hF + 1);
plot(1:N_ro, ro_our, 1:N_ro_1, ro_GPS_BoC_10_5, 1:N_ro_dop, cos_df_dop.*ro_our_dop, ...
     1:N_ro_dop, ro_GPS_BoC_10_5_dop_f, 1:N_ro_dop, ro_our_f_10_5);
title('All')

ff = (-(N_ro/2 - 1):1:N_ro/2)/(Td*1e6)/N_ro; % Ось частот для недополненной АКФ
ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td*1e6)/N_ro_dop;  % Ось частот для дополненной АКФ

hF = figure(hF + 1);
plot(ff_dop, (abs(fftshift(fft(ro_our_dop.*cos_df_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_GPS_BoC_10_5_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_our_f_10_5)))), ...
     ff_dop, (abs(fftshift(fft(ro_GPS_BoC_10_5_dop_f))))   )
xlabel('MHz')