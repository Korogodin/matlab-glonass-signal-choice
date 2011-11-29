%/**
% Данный скрипт рассчитывает автокорреляционную функцию низкочастотного 
% эквивалента (скорее даже огибающей) BoCsin, BoCcos или BPSK сигнала, а
% затем рассчитывает коэффициент кодового разделения k_cd. АКФ запоминается
% в mat-файле с соответствующим именем, её график в png и fig файлах. 
% Результат вычисления коэффициента кодового разделения сохраняется в
% массивах InSysJam_*** с индексами (8*m, 8*n) для BoC-сигналов или (8*n) 
% для BPSK-сигналов, которые сохраняются в соответствующих mat-файлах. 
%
%@param Signal_Type задает BoCsin, BoCcos или BPSK
%@param m8max, n8max - максимальное значение индексов 8*m, 8*n. Пробегаем
%по n и m  с шагом 0.125
%*/

clear
close all
clc

n8max = 80;
m8max = 80;

% С помощью данной секции можно создать 
% InSysJam_BoCsin = nan(m8max, n8max);
% InSysJam_BoCcos = nan(m8max, n8max);
% InSysJam_BPSK = nan(m8max, n8max);

load([pwd '/results/intrasystem/InSysJam_BoCsin.mat']);
load([pwd '/results/intrasystem/InSysJam_BoCcos.mat']);
load([pwd '/results/intrasystem/InSysJam_BPSK.mat']);

n_min = 1;
fc = 1023;
N_ms = 8;
N_1ms = 2*496*fc*n_min;
N = N_1ms*N_ms;
Td = 1 / N_1ms / 1000;

hF = 0;
% BOC(m, n); m >= n
% BPSK(n)
Signal_Type = 	2; % 1 - BOCsin; 2 - BOCcos; 3 - BPSK;

Signals_L1;
for n8 = 1:1:(8*10)
    for m8 = 1:1:(8*10) 

% for n8 = 1:1:80
%     for m8 = 0

% for ii = 1:size(BoCcos_Freq_L1_num, 1)        
%     for jj = 1
%     m8 = BoCcos_Freq_L1_num(ii, 1)*8;
%     n8 = BoCcos_Freq_L1_num(ii, 2)*8;
       

        if m8 < n8
            if Signal_Type ~= 3
                continue;
            end
        end
%         if ((m8+n8)/8 > 12)  % Если этот сигнал не влазиет в полосу
%             continue;
%         end

        % Если уже посчитано, то идем дальше
        if Signal_Type == 1 
            if ~isnan(InSysJam_BoCsin(m8, n8))
                continue;
            end
        elseif Signal_Type == 2
            if ~isnan(InSysJam_BoCcos(m8, n8))
                continue;
            end
        elseif Signal_Type == 3 
            if ~isnan(InSysJam_BPSK(n8))
                continue;
            end
        end
            
        m = m8 / 8;
        n = n8 / 8;
        
        t = (1:N)/N_1ms;
        if (Signal_Type == 1)
            g_boc = sign(sin(2*pi*m*fc*t) + eps);
        elseif Signal_Type == 2
            g_boc = sign(cos(2*pi*m*fc*t) + eps);
        elseif (Signal_Type == 3)
            g_boc = ones(1, length(t));            
        end
        g_bpsk = sign(randn(1, fc*n*N_ms));

        g_boc_N = resize_array(g_boc, N);
        g_bpsk_N = resize_array(g_bpsk, N);

        g_N = g_boc_N.*g_bpsk_N;

        N_co = fix(N_1ms / (1023 * n_min)); % Точек на чип для BPSK(n_min)
        N_main = ceil(N_1ms / (1023 * n)); % Точек на чип для выбранного n
        if N_main > N_co
            N_co = N_main;
        end
        ro_half = zeros(1, N_co); % Половинка от корреляционной ф-ии 
        ro = zeros(1, 2*N_co + 1); % сдвиг на N_co <-> 1 чипу bpsk <-> ro = 0;
        for i = 1:N_co
            if i < N_main
                for j = 1:N
                    if (j + i) <= N
                        ro_half(i) = ro_half(i) + g_N(j) * g_N(j + i);
                    else
                        ro_half(i) = ro_half(i) + g_N(j) * g_N(j + i - N);
                    end
                end
            else
                ro_half(i) = 0;
            end
            if ~mod(i, round(N_co/10))
                fprintf('Progress m = %.3f n = %.3f: %.0f%%\n', m, n, i/N_co*100);
            end
        end

        for i = 1:N_co
            ro(i) = ro_half(N_co - i + 1);
            ro(2*N_co + 2 - i) = ro(i);
        end
            ro(N_co + 1) = g_N*g_N';
            ro = ro / N;

        hF = 0;
        % hF = figure(hF + 1);
        % plot(t, g_N);
        % xlabel('t, ms')
        % ylabel('PRN')
        % 
        % hF = figure(hF + 1);
        % ff = (-(N/2 - 1):1:N/2)/(Td*1e6)/N;
        % plot(ff, 20*log10(abs(fftshift(fft(g_N / N)))));
        % xlabel('freq, MHz')
        % ylabel('Code FFT, dB')
        % 
        hF = figure(hF + 1);
        tau = (-N_co:1:N_co)*Td*1e6;
        plot(tau, ro);
        xlabel('\tau, \mu{s}')
        ylabel('\rho(\tau)')
        grid on;
        drawnow
        if (Signal_Type == 1) || (Signal_Type == 3)
            saveas(hF, [pwd '/ro/png/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
            saveas(hF, [pwd '/ro/fig/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').fig']);
        elseif Signal_Type == 2
            saveas(hF, [pwd '/ro/png/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').png']);
            saveas(hF, [pwd '/ro/fig/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').fig']);
        end
        % 
        % hF = figure(hF + 1);
        % ff2 = (-N_co:1:N_co)/(Td*1e6)/N_co;
        % plot(ff2, 10*log10(abs(fftshift(fft(ro)))));
        % xlabel('freq, MHz')
        % ylabel('ro FFT, dB')

        if Signal_Type == 3
            InSysJam_BPSK(n8) = 10*log10(ro*ro'*Td);
            fprintf('k_cd для BPSK(%.3f) составляет %.1f дБ \n', n, round(10 * InSysJam_BPSK(n8)) /10);
            save([pwd '/ro/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', n) ').mat'], 'ro');
            save([pwd '/results/intrasystem/InSysJam_BPSK.mat'], 'InSysJam_BPSK');           
        elseif Signal_Type == 1
            InSysJam_BoCsin(m8, n8) = 10*log10(ro*ro'*Td);
            fprintf('k_cd для BoCsin(%.3f, %.3f) составляет %.1f дБ \n', m, n, round(10 * InSysJam_BoCsin(m8, n8)) /10);
            save([pwd '/ro/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').mat'], 'ro');
            save([pwd '/results/intrasystem/InSysJam_BoCsin.mat'], 'InSysJam_BoCsin');
        elseif Signal_Type == 2
            InSysJam_BoCcos(m8, n8) = 10*log10(ro*ro'*Td);
            fprintf('k_cd для BoCcos(%.3f, %.3f) составляет %.1f дБ \n', m, n, round(10 * InSysJam_BoCcos(m8, n8)) /10);
            save([pwd '/ro/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').mat'], 'ro');
            save([pwd '/results/intrasystem/InSysJam_BoCcos.mat'], 'InSysJam_BoCcos');
        end

    end
end

if Signal_Type == 3
    hF = figure(hF + 1);
    index_bpsk = 8:8:80;
    index_n_bpsk = index_bpsk / 8;
    plot(index_n_bpsk, InSysJam_BPSK(index_bpsk))
    xlabel('n')
    ylabel('k_{cd}=k_{sd}, dB')
    grid on
elseif Signal_Type == 1
    Inn = InSysJam_BoCsin;
    Inn_ur = nan(80, 80);
    for i = 1:80
        for j = 1:80
            Inn_ur(i, j) = Inn(i, j);
        end
    end
    hF = figure(hF + 1);
    mesh((1:80)/8, (1:80)/8, Inn_ur)
    xlabel('n')
    ylabel('m')
    zlabel('k_{cd}=k_{sd}, dB')    
elseif Signal_Type == 2
    Inn = InSysJam_BoCcos;
    for i = 1:80
        for j = 1:80
            Inn_ur(i, j) = Inn(i, j);
        end
    end
    hF = figure(hF + 1);
    mesh((1:80)/8, (1:80)/8, Inn_ur)
    xlabel('n')
    ylabel('m')
    zlabel('k_{cd}=k_{sd}, dB')
end


if (Signal_Type == 1)
    hF = figure(hF+1);
    pcolor((1:80)/8, (1:80)/8, InSysJam_BoCsin(1:80,1:80));
    xlabel('n')
    ylabel('m')    
    title('k_{cd} = k_{sd}, dB')
    colorbar
    set(hF, 'Position', [1 1 620 485]);
elseif (Signal_Type == 2)
    hF = figure(hF+1);
    pcolor((1:80)/8, (1:80)/8, InSysJam_BoCcos(1:80,1:80));
    xlabel('n')
    ylabel('m')    
    title('k_{cd} = k_{sd}, dB')
    colorbar
    set(hF, 'Position', [1 1 620 485]);    
end

if (Signal_Type == 1)
    hF = figure(hF+1);
    plot((1:79)/8, InSysJam_BoCsin(80, 1:79));
    grid on 
    xlabel('n')
    ylabel('k_{cd} = k_{sd}, dB')    
elseif (Signal_Type == 2)
    hF = figure(hF+1);
    plot((1:80)/8, diag(InSysJam_BoCcos));
    grid on
    xlabel('n')
    ylabel('k_{cd} = k_{sd}, dB')    
elseif (Signal_Type == 3)
    hF = figure(hF+1);
    plot((1:80)/8, (InSysJam_BPSK));
    grid on
    xlabel('n')
    ylabel('k_{cd} = k_{sd}, dB')    
end

