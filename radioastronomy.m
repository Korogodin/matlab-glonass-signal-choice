%/**
% Данный скрипт предназначен для расчета максимального значения
% спектральной плотности мощности в диапазоне РА
%*/

clear 
close all
clc

path_to_ro = [pwd '/ro'];
path_to_results = [pwd '/results/radioastronomy'];

n8max = 80;
m8max = 80;

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
farr_Hz = farr * 1.023e6;

% Radioastronomy_BoCsin = nan(m8max, n8max, fmax);
% Radioastronomy_BoCcos = nan(m8max, n8max, fmax);
% Radioastronomy_BPSK = nan(n8max, fmax);
% save([path_to_results '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
% save([path_to_results '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
% save([path_to_results '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');
% return

load([path_to_results '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
load([path_to_results '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
load([path_to_results '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');

% Параметры нашего сигнала
BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

load([pwd '/ro/Td.mat']);

RA_1 = 1610.6e6; % RA bandspace
RA_2 = 1613.8e6;

for n8 = 1:80
    for m8 = 1:80  
         
        if m8 < n8
            continue;
        end
        
        % Если уже всё посчитано, то идем дальше
        if Signal_Type == BOCsin 
                if (~max(isnan(Radioastronomy_BoCsin(m8, n8, :)))) 
                    continue;
                end
        elseif Signal_Type == BOCcos
                if (~max(isnan(Radioastronomy_BoCcos(m8, n8, :))))
                    continue;
                end
        elseif Signal_Type == BPSK
                if (~max(isnan(Radioastronomy_BPSK(n8, :)))) 
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
        
        df = 1/(Td)/N_ro_dop;
        fft_ro_our_dop = abs(fftshift(fft(ro_our_dop))) / N_ro_dop;
        norm_power = sum(fft_ro_our_dop*df);
        fft_ro_our_dop = fft_ro_our_dop / norm_power;
        ff_dop = (-(N_ro_dop/2 - 1):1:N_ro_dop/2)/(Td)/N_ro_dop; % Ось частот для недополненной АКФ
        
        for freq_index = 1:fmax
            max_SP = -1;
            for jj = (N_ro_dop - 1)/2:N_ro_dop
                if ff_dop(jj) < (RA_1 - farr_Hz(freq_index))
                    continue
                end
                if ff_dop(jj) > (RA_2 - farr_Hz(freq_index))
                    break;
                end
                if max_SP < fft_ro_our_dop(jj)
                    max_SP = fft_ro_our_dop(jj);
                end
            end
            max_SP = 10*log10(max_SP);
            if Signal_Type == BOCsin 
                Radioastronomy_BoCsin(m8, n8, freq_index) = max_SP;
            elseif Signal_Type == BOCcos
                Radioastronomy_BoCcos(m8, n8, freq_index) = max_SP;
            elseif Signal_Type == BPSK
                Radioastronomy_BPSK(n8, freq_index) = max_SP;
            end

            fprintf('Radioastronomy fo BoC(%.3f, %.3f) at %.0f = %.2f dB\n', m, n, farr(freq_index), max_SP);
        end

        % Для BPSK по m пробегать не надо
        if (Signal_Type == BPSK)
            break;
        end
    end
    if Signal_Type == BOCsin 
        save([path_to_results '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
    elseif Signal_Type == BOCcos
        save([path_to_results '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
    elseif Signal_Type == BPSK
        save([path_to_results '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');
    end
end

hF = 0;

% hF = figure(hF + 1);
% for jj = 1:fmax
%     if Signal_Type == BOCsin
%         mesh((1:80)/8, (1:80)/8, Radioastronomy_BoCsin)
%         xlabel('n')
%         ylabel('m')
%         zlabel('Max power density in RA span, dB')    
%     elseif Signal_Type == BOCcos
%         mesh((1:80)/8, (1:80)/8, Radioastronomy_BoCcos)
%         xlabel('n')
%         ylabel('m')
%         zlabel('Max power density in RA span, dB')    
%     end
% end

for jj = 1:fmax
    if Signal_Type == BOCsin
        hF = figure(hF + 1);
        pcolor((1:80)/8, (1:80)/8, Radioastronomy_BoCsin(:, :, jj))
        xlabel('n')
        ylabel('m')
        zlabel('Max power density in RA span, dB')    
    elseif Signal_Type == BOCcos
        hF = figure(hF + 1);
        pcolor((1:80)/8, (1:80)/8, Radioastronomy_BoCcos(:, :, jj))
        xlabel('n')
        ylabel('m')
        zlabel('Max power density in RA span, dB')  
    elseif Signal_Type == BPSK
        hF = figure(hF + 1);
        plot(1:80, Radioastronomy_BPSK(:, jj))
        xlabel('n')
        ylabel('Max power density in RA span, dB')    
        grid on
    end
    title(sprintf('Normalized freq = %.0f', farr(jj)));
end


hF = figure(hF + 1);
min_lev = min(10*log10(fft_ro_our_dop));
max_lev = max(10*log10(fft_ro_our_dop));
plot(ff_dop/1e6, 10*log10(fft_ro_our_dop), ...
    ([RA_1 RA_1] - farr_Hz(freq_index))/1e6, [min_lev max_lev], 'r', ...
    ([RA_2 RA_2] - farr_Hz(freq_index))/1e6, [min_lev max_lev], 'r');
xlabel('f, MHz')
ylabel('fft, dB')
