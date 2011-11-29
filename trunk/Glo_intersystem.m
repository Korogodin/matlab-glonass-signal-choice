clear 
close all
clc

n8max = 80;
m8max = 80;
farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
lit_arr = -7:6; lit_size = length(lit_arr);

% Empty file creation
% InterSysJam_BoCsin_L1_GloST = nan(m8max, n8max, fmax, lit_size);
% InterSysJam_BoCcos_L1_GloST = nan(m8max, n8max, fmax, lit_size);
% save('results/InterSysJam_BoCsin_L1_GloST.mat', 'InterSysJam_BoCsin_L1_GloST');
% save('results/InterSysJam_BoCcos_L1_GloST.mat', 'InterSysJam_BoCcos_L1_GloST');
% InterSysJam_BoCsin_L1_GloVT = nan(m8max, n8max, fmax, lit_size);
% InterSysJam_BoCcos_L1_GloVT = nan(m8max, n8max, fmax, lit_size);
% save('results/InterSysJam_BoCsin_L1_GloVT.mat', 'InterSysJam_BoCsin_L1_GloVT');
% save('results/InterSysJam_BoCcos_L1_GloVT.mat', 'InterSysJam_BoCcos_L1_GloVT');
% InterSysJam_BPSK_L1_GloST = nan(n8max, fmax, lit_size);
% InterSysJam_BPSK_L1_GloVT = nan(n8max, fmax, lit_size);
% save('results/InterSysJam_BPSK_L1_GloST.mat', 'InterSysJam_BPSK_L1_GloST');
% save('results/InterSysJam_BPSK_L1_GloVT.mat', 'InterSysJam_BPSK_L1_GloVT');
% 
% InterSysJam_BoCsin_L1_GloST_mean = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L1_GloST_mean = nan(m8max, n8max, fmax);
% save('results/InterSysJam_BoCsin_L1_GloST_mean.mat', 'InterSysJam_BoCsin_L1_GloST_mean');
% save('results/InterSysJam_BoCcos_L1_GloST_mean.mat', 'InterSysJam_BoCcos_L1_GloST_mean');
% InterSysJam_BoCsin_L1_GloVT_mean = nan(m8max, n8max, fmax);
% InterSysJam_BoCcos_L1_GloVT_mean = nan(m8max, n8max, fmax);
% save('results/InterSysJam_BoCsin_L1_GloVT_mean.mat', 'InterSysJam_BoCsin_L1_GloVT_mean');
% save('results/InterSysJam_BoCcos_L1_GloVT_mean.mat', 'InterSysJam_BoCcos_L1_GloVT_mean');
% InterSysJam_BPSK_L1_GloST_mean = nan(n8max, fmax);
% InterSysJam_BPSK_L1_GloVT_mean = nan(n8max, fmax);
% save('results/InterSysJam_BPSK_L1_GloST_mean.mat', 'InterSysJam_BPSK_L1_GloST_mean');
% save('results/InterSysJam_BPSK_L1_GloVT_mean.mat', 'InterSysJam_BPSK_L1_GloVT_mean');


% Load files
load('results/InterSysJam_BoCsin_L1_GloST.mat', 'InterSysJam_BoCsin_L1_GloST');
load('results/InterSysJam_BoCcos_L1_GloST.mat', 'InterSysJam_BoCcos_L1_GloST');
load('results/InterSysJam_BoCsin_L1_GloVT.mat', 'InterSysJam_BoCsin_L1_GloVT');
load('results/InterSysJam_BoCcos_L1_GloVT.mat', 'InterSysJam_BoCcos_L1_GloVT');
load('results/InterSysJam_BPSK_L1_GloST.mat', 'InterSysJam_BPSK_L1_GloST');
load('results/InterSysJam_BPSK_L1_GloVT.mat', 'InterSysJam_BPSK_L1_GloVT');
load('results/InterSysJam_BoCsin_L1_GloST_mean.mat', 'InterSysJam_BoCsin_L1_GloST_mean');
load('results/InterSysJam_BoCcos_L1_GloST_mean.mat', 'InterSysJam_BoCcos_L1_GloST_mean');
load('results/InterSysJam_BoCsin_L1_GloVT_mean.mat', 'InterSysJam_BoCsin_L1_GloVT_mean');
load('results/InterSysJam_BoCcos_L1_GloVT_mean.mat', 'InterSysJam_BoCcos_L1_GloVT_mean');
load('results/InterSysJam_BPSK_L1_GloST_mean.mat', 'InterSysJam_BPSK_L1_GloST_mean');
load('results/InterSysJam_BPSK_L1_GloVT_mean.mat', 'InterSysJam_BPSK_L1_GloVT_mean');


% Параметры нашего сигнала
BOC_Type = 1; % 1 - sin, 2 - cos, 3 - BPSK
load([pwd '/ro/Td.mat']);

% АКФ сигнала GLO BPSK(0.5)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', 0.5) ').mat'])
ro_GLO_ST = ro; 

% АКФ сигнала GLO BPSK(5)
load([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', 5) ').mat'])
ro_GLO_VT = ro; 

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
%             for m8 = 0
            Inte_Glo_ST_nodB  = zeros(1, lit_size);
            Inte_Glo_VT_nodB  = zeros(1, lit_size);
            for lit_index = 1:lit_size
         
                do_out = 1;
                if m8 < n8
                    do_out = 0;
                    break;
                end
                
                % Если уже посчитано, то идем дальше
                if BOC_Type == 1 
                    if (~isnan(InterSysJam_BoCsin_L1_GloST(m8, n8, f_index, lit_index))) && ...
                            (~isnan(InterSysJam_BoCsin_L1_GloVT(m8, n8, f_index, lit_index)))
                        do_out = 0;
                        break;
                    end
                elseif BOC_Type == 2
                    if (~isnan(InterSysJam_BoCcos_L1_GloST(m8, n8, f_index, lit_index))) && ...
                            (~isnan(InterSysJam_BoCcos_L1_GloVT(m8, n8, f_index, lit_index)))
                        do_out = 0;
                        break;
                    end
                elseif BOC_Type == 3
                    if (~isnan(InterSysJam_BPSK_L1_BoC_GloST(n8, f_index, lit_index))) && ...
                            (~isnan(InterSysJam_BPSK_L1_GloVT(n8, f_index, lit_index)))
                        do_out = 0;
                        break;
                    end                
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
                    do_out = 0;
                    break; % Если файла нет
                end
                ro_our = ro;
                N_ro = length(ro);
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
                
                % Конвертер частоты
                f_sig = 1.023e6*farr(f_index); % Частота нашего сигнала     
                f_intro = 1602e6 + lit_arr(lit_index)*0.5625e6;
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

                if BOC_Type == 1 
                    InterSysJam_BoCsin_L1_GloST(m8, n8, f_index, lit_index) = Inte_Glo_ST;
                    InterSysJam_BoCsin_L1_GloVT(m8, n8, f_index, lit_index) = Inte_Glo_VT;
                elseif BOC_Type == 2
                    InterSysJam_BoCcos_L1_GloST(m8, n8, f_index, lit_index) = Inte_Glo_ST;
                    InterSysJam_BoCcos_L1_GloVT(m8, n8, f_index, lit_index) = Inte_Glo_VT;
                elseif BOC_Type == 3
                    InterSysJam_BPSK_L1_GloST(n8, f_index, lit_index) = Inte_Glo_ST;
                    InterSysJam_BPSK_L1_GloVT(n8, f_index, lit_index) = Inte_Glo_VT;
                end

            end

            if do_out 
                Inte_Glo_ST_mean = 10*log10(mean(Inte_Glo_ST_nodB));
                Inte_Glo_VT_mean = 10*log10(mean(Inte_Glo_VT_nodB));
                fprintf('Intersystem Jamm BoC(%.3f, %.3f) at %.0f \n \t with GloST = %.2f dB\n \t with GloVT = %.2f dB\n', m, n, farr(f_index), ...
                    Inte_Glo_ST_mean, Inte_Glo_VT_mean);

                if BOC_Type == 1 
                    InterSysJam_BoCsin_L1_GloST_mean(m8, n8, f_index) = Inte_Glo_ST_mean;
                    InterSysJam_BoCsin_L1_GloVT_mean(m8, n8, f_index) = Inte_Glo_VT_mean;
                elseif BOC_Type == 2
                    InterSysJam_BoCcos_L1_GloST_mean(m8, n8, f_index) = Inte_Glo_ST_mean;
                    InterSysJam_BoCcos_L1_GloVT_mean(m8, n8, f_index) = Inte_Glo_VT_mean;
                elseif BOC_Type == 3
                    InterSysJam_BPSK_L1_GloST_mean(n8, f_index) = Inte_Glo_ST_mean;
                    InterSysJam_BPSK_L1_GloVT_mean(n8, f_index) = Inte_Glo_VT_mean;
                end        
            end
        end
        if BOC_Type == 1 
            save('results/InterSysJam_BoCsin_L1_GloST.mat', 'InterSysJam_BoCsin_L1_GloST');
            save('results/InterSysJam_BoCsin_L1_GloVT.mat', 'InterSysJam_BoCsin_L1_GloVT');
            save('results/InterSysJam_BoCsin_L1_GloST_mean.mat', 'InterSysJam_BoCsin_L1_GloST_mean');
            save('results/InterSysJam_BoCsin_L1_GloVT_mean.mat', 'InterSysJam_BoCsin_L1_GloVT_mean');
        elseif BOC_Type == 2
            save('results/InterSysJam_BoCcos_L1_GloST.mat', 'InterSysJam_BoCcos_L1_GloST');
            save('results/InterSysJam_BoCcos_L1_GloVT.mat', 'InterSysJam_BoCcos_L1_GloVT');
            save('results/InterSysJam_BoCcos_L1_GloST_mean.mat', 'InterSysJam_BoCcos_L1_GloST_mean');
            save('results/InterSysJam_BoCcos_L1_GloVT_mean.mat', 'InterSysJam_BoCcos_L1_GloVT_mean');
        elseif BOC_Type == 3
            save('results/InterSysJam_BPSK_L1_GloST.mat', 'InterSysJam_BPSK_L1_GloST');
            save('results/InterSysJam_BPSK_L1_GloVT.mat', 'InterSysJam_BPSK_L1_GloVT');
            save('results/InterSysJam_BPSK_L1_GloST_mean.mat', 'InterSysJam_BPSK_L1_GloST_mean');
            save('results/InterSysJam_BPSK_L1_GloVT_mean.mat', 'InterSysJam_BPSK_L1_GloVT_mean');           
        end                
    end
end
hF = 0;


for i = 1:fmax
    hF = figure(hF+1);
    if (BOC_Type == 1)
        pcolor((1:80)/8, (1:80)/8, InterSysJam_BoCsin_L1_GloST_mean(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (BOC_Type == 2)
        pcolor((1:80)/8, (1:80)/8, InterSysJam_BoCcos_L1_GloST_mean(1:80,1:80, i));
        xlabel('n')
        ylabel('m')
    elseif (BOC_Type == 3)
        plot((1:80)/8, InterSysJam_BPSK_L1_GloST_mean(1:80, i));        
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