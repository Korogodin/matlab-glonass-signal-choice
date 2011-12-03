%/**
% Данный скрипт предназначен для расчета показателя beta, определяющего 
% потенциальную точность оценки задержки
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/beta'];

n8max = 80;
m8max = 80;

% Signals_L1; % Параметры приемлимых сигналов

% Beta_BoCsin = nan(m8max, n8max);
% Beta_BoCcos = nan(m8max, n8max);
% Beta_BPSK = nan(n8max);
% save([path_to_results '/Beta_BoCsin.mat'], 'Beta_BoCsin');
% save([path_to_results '/Beta_BoCcos.mat'], 'Beta_BoCcos');
% save([path_to_results '/Beta_BPSK.mat'], 'Beta_BPSK');

load([path_to_results '/Beta_BoCsin.mat'], 'Beta_BoCsin');
load([path_to_results '/Beta_BoCcos.mat'], 'Beta_BoCcos');
load([path_to_results '/Beta_BPSK.mat'], 'Beta_BPSK');

% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

load([pwd '/ro/Td.mat']);

for n8 = 1:80
    for m8 = 1:80  
         
        if m8 < n8
            continue;
        end
        
        % Если уже посчитано, то идем дальше
        if Signal_Type == BOCsin 
                if (~isnan(Beta_BoCsin(m8, n8))) 
                    continue;
                end
        elseif Signal_Type == BOCcos
                if (~isnan(Beta_BoCcos(m8, n8)))
                    continue;
                end
        elseif Signal_Type == BPSK
                if (~isnan(Beta_BPSK(n8))) 
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
        N_ro_dop = N_ro*50 + 1;
        ro_our_dop = zeros(1, N_ro_dop);
        ro_our_dop( ( (N_ro_dop - 1)/2 - (N_ro - 1)/2 ):( (N_ro_dop - 1)/2 + (N_ro - 1)/2 ) ) = ...
            ro_our;
       
        Bahd_f = (m+n)*2*1.023e6; % Полоса сигнала
        df = 1/(Td)/N_ro_dop;
        N_int = ceil(Bahd_f / df);
        
        fft_ro_our_dop = abs(fftshift(fft(ro_our_dop))) / N_ro_dop;
        norm_power = sum(fft_ro_our_dop*df);
        fft_ro_our_dop = fft_ro_our_dop / norm_power;
        ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td)/N_ro_dop; % Ось частот для недополненной АКФ
        ind_int = ( (N_ro_dop - 1)/2 - fix((N_int - 1)/2) ):( (N_ro_dop - 1)/2 + fix((N_int - 1)/2) );        
        Window = (2*pi*ff_dop(ind_int)).^2;
        
        Beta2 = fft_ro_our_dop(ind_int) * Window' * df;

        if Signal_Type == BOCsin 
            Beta_BoCsin(m8, n8) = sqrt(Beta2);
        elseif Signal_Type == BOCcos
            Beta_BoCcos(m8, n8) = sqrt(Beta2);
        elseif Signal_Type == BPSK
            Beta_BPSK(n8) = sqrt(Beta2);
        end

        fprintf('Beta fo BoC(%.3f, %.3f)  = %.2f MHz\n', m, n, sqrt(Beta2)/1e6);

        % Для BPSK по m пробегать не надо
        if (Signal_Type == BPSK)
            break;
        end
    end
    if Signal_Type == BOCsin 
        save([path_to_results '/Beta_BoCsin.mat'], 'Beta_BoCsin');
    elseif Signal_Type == BOCcos
        save([path_to_results '/Beta_BoCcos.mat'], 'Beta_BoCcos');
    elseif Signal_Type == BPSK
        save([path_to_results '/Beta_BPSK.mat'], 'Beta_BPSK');
    end
end

hF = 4;

hF = figure(hF + 1);
if Signal_Type == BOCsin
    mesh((1:80)/8, (1:80)/8, Beta_BoCsin/1e6)
    xlabel('n')
    ylabel('m')
    zlabel('\beta, MHz')    
elseif Signal_Type == BOCcos
    mesh((1:80)/8, (1:80)/8, Beta_BoCcos/1e6)
    xlabel('n')
    ylabel('m')
    zlabel('\beta, MHz')    
elseif Signal_Type == BPSK
    plot(1:80, Beta_BPSK/1e6)
    xlabel('n')
    ylabel('\beta, MHz')    
    grid on
end


if Signal_Type == BOCsin
    hF = figure(hF + 1);
    pcolor((1:80)/8, (1:80)/8, Beta_BoCsin/1e6)
    xlabel('n')
    ylabel('m')
    zlabel('\beta, MHz')    
elseif Signal_Type == BOCcos
    hF = figure(hF + 1);
    pcolor((1:80)/8, (1:80)/8, Beta_BoCcos/1e6)
    xlabel('n')
    ylabel('m')
    zlabel('\beta, MHz')    
end


hF = figure(hF + 1);
plot(ff_dop/1e6, fft_ro_our_dop/max(fft_ro_our_dop), ff_dop(ind_int) / 1e6, Window / max(Window))
title('Normalized to 1')
xlabel('f, MHz')
ylabel('fft')
