%/**
% Данный скрипт предназначен для расчета коэффициента спектрального
% разделения между сигналами полного выборочного множества и полосовым
% шумом
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/pomexoyst'];

n8max = 80;
m8max = 80;

% Signals_L1; % Параметры приемлимых сигналов

% Pomexoyst_BoCsin = nan(m8max, n8max);
% Pomexoyst_BoCcos = nan(m8max, n8max);
% Pomexoyst_BPSK = nan(n8max);
% save([path_to_results '/Pomexoyst_BoCsin.mat'], 'Pomexoyst_BoCsin');
% save([path_to_results '/Pomexoyst_BoCcos.mat'], 'Pomexoyst_BoCcos');
% save([path_to_results '/Pomexoyst_BPSK.mat'], 'Pomexoyst_BPSK');

load([path_to_results '/Pomexoyst_BoCsin.mat'], 'Pomexoyst_BoCsin');
load([path_to_results '/Pomexoyst_BoCcos.mat'], 'Pomexoyst_BoCcos');
load([path_to_results '/Pomexoyst_BPSK.mat'], 'Pomexoyst_BPSK');

% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

load([pwd '/ro/Td.mat']);

for n8 = 1:80
    for m8 = 1:80  
         
        if m8 < n8
            continue;
        end
        
        % Если уже посчитано, то идем дальше
        if Signal_Type == BOCsin 
                if (~isnan(Pomexoyst_BoCsin(m8, n8))) 
                    continue;
                end
        elseif Signal_Type == BOCcos
                if (~isnan(Pomexoyst_BoCcos(m8, n8)))
                    continue;
                end
        elseif Signal_Type == BPSK
                if (~isnan(Pomexoyst_BPSK(n8))) 
                    continue;
                end                
        end

        m = m8/8 * (Signal_Type ~= BPSK);
        n = n8/8;

        % Открываем файл АКФ указанного сигнала
        ro_our = get_ro(m, n, Signal_Type, path_to_ro);
        if ro_our == 0 % Если нет такого
            continue
        end
        N_ro = length(ro_our);
        N_ro_dop = N_ro*4 + 1;
        ro_our_dop = zeros(1, N_ro_dop);
        ro_pomex_dop = zeros(1, N_ro_dop);
        ro_our_dop(1:N_ro) = ro_our;

        offset = (N_ro - 1) / 2;
        delta_f = (m+n)*2*1.023e6; % Полоса сигнала
        tau_p = ((1:N_ro) - offset - 1)*Td;
        ro_pomex = sinc(pi * delta_f * tau_p /pi);
        ro_pomex_dop(1:N_ro) = ro_pomex;

%         ro_our_f_6_1 = filter(Hd_6_1, cos_df_dop.*ro_our_dop);
        ro_pomex_dop_f = ro_pomex_dop;
        ro_our_dop_f = ro_our_dop;

        % Расчет интегралов - коэффициентов спектрального разделения
        k_cd  = 10*log10(abs(ro_our_dop_f*(ro_pomex_dop_f)'*Td));

        if Signal_Type == BOCsin 
            Pomexoyst_BoCsin(m8, n8) = k_cd;
        elseif Signal_Type == BOCcos
            Pomexoyst_BoCcos(m8, n8) = k_cd;
        elseif Signal_Type == BPSK
            Pomexoyst_BPSK(n8) = k_cd;
        end

        fprintf('k_cd fo BoC(%.3f, %.3f)  = %.2f dB\n', m, n, k_cd);

        % Для BPSK по m пробегать не надо
        if (Signal_Type == BPSK)
            break;
        end
    end
    if Signal_Type == BOCsin 
        save([path_to_results '/Pomexoyst_BoCsin.mat'], 'Pomexoyst_BoCsin');
    elseif Signal_Type == BOCcos
        save([path_to_results '/Pomexoyst_BoCcos.mat'], 'Pomexoyst_BoCcos');
    elseif Signal_Type == BPSK
        save([path_to_results '/Pomexoyst_BPSK.mat'], 'Pomexoyst_BPSK');
    end
end

hF = 0;

hF = figure(hF + 1);
if Signal_Type == BOCsin
    mesh((1:80)/8, (1:80)/8, Pomexoyst_BoCsin)
    xlabel('n')
    ylabel('m')
    zlabel('k_{cd}=k_{sd}, dB')    
elseif Signal_Type == BOCcos
    mesh((1:80)/8, (1:80)/8, Pomexoyst_BoCcos)
    xlabel('n')
    ylabel('m')
    zlabel('k_{cd}=k_{sd}, dB')
elseif Signal_Type == BPSK
    plot(1:80, Pomexoyst_BPSK)
    xlabel('n')
    ylabel('k_{cd}=k_{sd}, dB')
    grid on
end


if Signal_Type == BOCsin
    hF = figure(hF + 1);
    pcolor((1:80)/8, (1:80)/8, Pomexoyst_BoCsin)
    xlabel('n')
    ylabel('m')
    zlabel('k_{cd}=k_{sd}, dB')    
elseif Signal_Type == BOCcos
    hF = figure(hF + 1);
    pcolor((1:80)/8, (1:80)/8, Pomexoyst_BoCcos)
    xlabel('n')
    ylabel('m')
    zlabel('k_{cd}=k_{sd}, dB')
end


hF = figure(hF + 1);
plot(1:N_ro_dop, ro_our_dop, 1:N_ro_dop, ro_pomex_dop);
title('All')

ff = (-(N_ro/2 - 1):1:N_ro/2)/(Td*1e6)/N_ro; % Ось частот для недополненной АКФ
ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td*1e6)/N_ro_dop;  % Ось частот для дополненной АКФ

hF = figure(hF + 1);
plot(ff_dop, (abs(fftshift(fft(ro_our_dop)))), ...
     ff_dop, (abs(fftshift(fft(ro_pomex_dop)))) )
xlabel('MHz')