%/**
% Скрипт формирования усеченного множества сигналов
%*/

clear 
close all
clc

path_to_results_RA = [pwd '/results/radioastronomy'];
load([path_to_results_RA '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
load([path_to_results_RA '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
load([path_to_results_RA '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');

path_to_results_pomexoyst = [pwd '/results/pomexoyst'];
load([path_to_results_pomexoyst '/Pomexoyst_BoCsin.mat'], 'Pomexoyst_BoCsin');
load([path_to_results_pomexoyst '/Pomexoyst_BoCcos.mat'], 'Pomexoyst_BoCcos');
load([path_to_results_pomexoyst '/Pomexoyst_BPSK.mat'], 'Pomexoyst_BPSK');

path_to_results_beta = [pwd '/results/beta'];
load([path_to_results_beta '/Beta_BoCsin.mat'], 'Beta_BoCsin');
load([path_to_results_beta '/Beta_BoCcos.mat'], 'Beta_BoCcos');
load([path_to_results_beta '/Beta_BPSK.mat'], 'Beta_BPSK');

path_to_results_intra = [pwd '/results/intrasystem'];
load([path_to_results_intra '/InSysJam_BoCsin.mat'], 'InSysJam_BoCsin');
load([path_to_results_intra '/InSysJam_BoCcos.mat'], 'InSysJam_BoCcos');
load([path_to_results_intra '/InSysJam_BPSK.mat'], 'InSysJam_BPSK');


path_to_results_compl = [pwd '/results/complexity_L1'];
load([path_to_results_compl '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
load([path_to_results_compl '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
load([path_to_results_compl '/Complexity_BPSK.mat'], 'Complexity_BPSK');


BOCsin = 1; BOCcos = 2; BPSK = 3;
Signal_Type = 3; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

OpenSignal = 1; CloseSignal = 2; % Сигнал ОД или СД
Access_Type = 2;

if Access_Type == OpenSignal
    Access_str = 'open';
    Threshold_RA = -84; % Порог по РА
    Threshold_Pomex = Pomexoyst_BPSK(0.5*8);
    Threshold_Accur = 1/Beta_BPSK(0.5*8);
    Threshold_Intra = InSysJam_BPSK(0.5*8);
    Threshold_Complexity = Complexity_BoCsin(1*8, 1*8, 3);
elseif Access_Type == CloseSignal
    Access_str = 'secure';
    Threshold_RA = -90; % Порог по РА
    Threshold_Pomex = Pomexoyst_BPSK(5*8);
    Threshold_Accur = 1/Beta_BPSK(5*8);
    Threshold_Intra = InSysJam_BPSK(5*8);
    Threshold_Complexity = Complexity_BoCsin(10*8, 5*8, 3);
end

    fprintf('{| class="wikitable sortable" border="1" \n');
if Signal_Type == BOCsin   
    fprintf('|+ BOC<sub>sin</sub>(m, n) %s signals after clipping \n', Access_str);
elseif Signal_Type == BOCcos
    fprintf('|+ BOC<sub>cos</sub>(m, n) %s signals after clipping \n', Access_str);
elseif Signal_Type == BPSK
    fprintf('|+ BPSK(n) %s signals after clipping \n', Access_str);
end    

    fprintf('|- align="center"\n');
    fprintf('!Type \n');
    if Signal_Type ~= BPSK
        fprintf('!''''m''''\n');
    end
    fprintf('!''''n''''\n');
    fprintf('!''''f''''<sub>n</sub>\n');
    fprintf('!''''PSD''''<sub>max</sub>, dB \n');
    fprintf('!''''k''''<sub>cd,pomex</sub>, dB \n');
    fprintf('!1/<math>\\beta</math>, nm \n')
    fprintf('!''''k''''<sub>cd,intra</sub>, dB \n')
    fprintf('!Complexity, scores')
    
counter = 0;    
for n8 = 1:80
    for m8 = 1:80
        if m8 < n8
            if Signal_Type ~= BPSK
                continue;
            end
        end
        n = n8 / 8;
        m = m8 / 8 * (Signal_Type ~= BPSK);
        for freq_index = 1:fmax
            if Signal_Type == BOCsin
                if Radioastronomy_BoCsin(m8, n8, freq_index) <= Threshold_RA
                    if Pomexoyst_BoCsin(m8, n8) <= Threshold_Pomex
                        if (1/Beta_BoCsin(m8, n8)) <= Threshold_Accur
                            if InSysJam_BoCsin(m8, n8) <= Threshold_Intra
                                if Complexity_BoCsin(m8, n8, freq_index) <= Threshold_Complexity
                                    if ((m+n + farr(freq_index)) <= (farr(end) + 1)) && ((farr(freq_index) - (m+n)) >= (farr(1) - 1))
                                        fprintf('\n|- align="center"\n');
                                        fprintf('| BOC<sub>sin</sub>');
                                        fprintf('|| %.3f ', m8/8);
                                        fprintf('|| %.3f ', n8/8);
                                        fprintf('|| %4.0f ', farr(freq_index));
                                        fprintf('|| %3.1f ', Radioastronomy_BoCsin(m8, n8, freq_index));
                                        fprintf('|| %.1f ', Pomexoyst_BoCsin(m8, n8));
                                        fprintf('|| %.2f ', 1/Beta_BoCsin(m8, n8)*1e9);
                                        fprintf('|| %.1f ', InSysJam_BoCsin(m8, n8));
                                        fprintf('|| %.1f ', Complexity_BoCsin(m8, n8, freq_index));
                                        counter = counter + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            elseif Signal_Type == BOCcos
                if Radioastronomy_BoCcos(m8, n8, freq_index) <= Threshold_RA
                    if Pomexoyst_BoCcos(m8, n8) <= Threshold_Pomex
                        if (1/Beta_BoCcos(m8, n8)) <= Threshold_Accur
                            if InSysJam_BoCcos(m8, n8) <= Threshold_Intra
                                if Complexity_BoCcos(m8, n8, freq_index) <= Threshold_Complexity
                                    if ((m+n + farr(freq_index)) <= (farr(end) + 1)) && ((farr(freq_index) - (m+n)) >= (farr(1) - 1))
                                        fprintf('\n|- align="center"\n');
                                        fprintf('| BOC<sub>cos</sub>');
                                        fprintf('|| %.3f ', m8/8);
                                        fprintf('|| %.3f ', n8/8);
                                        fprintf('|| %4.0f ', farr(freq_index));
                                        fprintf('|| %3.1f ', Radioastronomy_BoCcos(m8, n8, freq_index));
                                        fprintf('|| %.1f ', Pomexoyst_BoCcos(m8, n8));
                                        fprintf('|| %.2f ', 1/Beta_BoCcos(m8, n8)*1e9);
                                        fprintf('|| %.1f ', InSysJam_BoCcos(m8, n8));
                                        fprintf('|| %.1f ', Complexity_BoCcos(m8, n8, freq_index));
                                        counter = counter + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            elseif Signal_Type == BPSK
                if Radioastronomy_BPSK(n8, freq_index) <= Threshold_RA
                    if Pomexoyst_BPSK(n8) <= Threshold_Pomex
                        if (1/Beta_BPSK(n8)) <= Threshold_Accur
                            if InSysJam_BPSK(n8) <= Threshold_Intra
                                if Complexity_BPSK(n8, freq_index) <= Threshold_Complexity
                                    if ((m+n + farr(freq_index)) <= (farr(end) + 1)) && ((farr(freq_index) - (m+n)) >= (farr(1) - 1))
                                        fprintf('\n|- align="center"\n');
                                        fprintf('| BPSK');
                                        fprintf('|| %.3f ', n8/8);
                                        fprintf('|| %4.0f ', farr(freq_index));
                                        fprintf('|| %3.1f ', Radioastronomy_BPSK(n8, freq_index));
                                        fprintf('|| %.1f ', Pomexoyst_BPSK(n8));
                                        fprintf('|| %.2f ', 1/Beta_BPSK(n8)*1e9);
                                        fprintf('|| %.1f ', InSysJam_BPSK(n8));
                                        fprintf('|| %.1f ', Complexity_BPSK(n8, freq_index));
                                        counter = counter + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        if Signal_Type == BPSK
            break;
        end
    end
end
fprintf('\n|}\n');