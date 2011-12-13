close all
clc
clear
%/**
% Скрипт формирования усеченного множества сигналов
%*/

% path_to_results_RA = [pwd '/results/radioastronomy'];
% load([path_to_results_RA '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
% load([path_to_results_RA '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
% load([path_to_results_RA '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');

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


path_to_results_compl = [pwd '/results/complexity_L2'];
load([path_to_results_compl '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
load([path_to_results_compl '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
load([path_to_results_compl '/Complexity_BPSK.mat'], 'Complexity_BPSK');

path_to_Evgeniy = [pwd '/results/Evgeniy'];
load([path_to_Evgeniy '/AcqBOC.mat'], 'AcqBOC');
load([path_to_Evgeniy '/AcqBPSK.mat'], 'AcqBPSK');

path_to_Evgeniy = [pwd '/results/Evgeniy'];
load([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
load([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
load([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');

path_to_results_inter_C = [pwd '/results/intersystem_L2/Common'];
load([path_to_results_inter_C '/InterSysJam_BoCsin_L2.mat'], 'InterSysJam_BoCsin_L2');
load([path_to_results_inter_C '/InterSysJam_BoCcos_L2.mat'], 'InterSysJam_BoCcos_L2');
load([path_to_results_inter_C '/InterSysJam_BPSK_L2.mat'], 'InterSysJam_BPSK_L2');

path_to_results_back_inter_C = [pwd '/results/back_intersystem_L2/Common'];
load([path_to_results_back_inter_C '/BackInterSysJam_BoCsin_L2.mat'], 'BackInterSysJam_BoCsin_L2');
load([path_to_results_back_inter_C '/BackInterSysJam_BoCcos_L2.mat'], 'BackInterSysJam_BoCcos_L2');
load([path_to_results_back_inter_C '/BackInterSysJam_BPSK_L2.mat'], 'BackInterSysJam_BPSK_L2');


BOCsin = 1; BOCcos = 2; BPSK = 3;
% Signal_Type = 1; % 1 - BOCsin, 2 - BOCcos, 3 - BPSK

farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
m8max = 80; n8max = 80;

OpenSignal = 1; CloseSignal = 2; % Сигнал ОД или СД
% Access_Type = 2;


counter_all = 0;
for Access_Type = 2:2
    
    fprintf('{| class="wikitable sortable" border="1" \n');
    
    if Access_Type == OpenSignal
        fprintf('|+ L2 open signals after clipping \n');
        Threshold_Pomex = Pomexoyst_BPSK(0.5*8);
        Threshold_Accur = 1/Beta_BPSK(0.5*8);
        Threshold_Intra = InSysJam_BPSK(0.5*8);
        Threshold_Complexity = Complexity_BoCsin(1*8, 1*8, 8);
    %     Threshold_Search = max(Best_Search);
    elseif Access_Type == CloseSignal
        fprintf('|+ L2 special signals after clipping \n');
        Threshold_Pomex = Pomexoyst_BPSK(5*8);
        Threshold_Accur = 1/Beta_BPSK(5*8);
        Threshold_Intra = InSysJam_BPSK(5*8);
        Threshold_Complexity = Complexity_BoCsin(10*8, 5*8, 8);
    %     Threshold_Search = max(Best_Search);
    end           
    
    WideTable = 1;

    fprintf('|- align="center"\n');
    fprintf('! \n');
    fprintf('! \n');    
    fprintf('!  \n');
    fprintf('!  \n');
    fprintf('! Antijamming \n');
    fprintf('! Accurancy \n');
    fprintf('! Intra \n');
    if WideTable
        fprintf('! Inter to \n');
        fprintf('! Inter from \n');
        fprintf('! MultiPath \n');        
        fprintf('! Search \n');
    end
    fprintf('!Complexity\n');
    
    fprintf('|- align="center"\n');
    fprintf('!Type \n');
    fprintf('!''''m''''\n');    
    fprintf('!''''n''''\n');
    fprintf('!''''f''''<sub>n</sub>\n');
    fprintf('!''''k''''<sub>cd,pomex</sub>, dB \n');
    fprintf('!1/<math>\\beta</math>, nm \n');
    fprintf('!''''k''''<sub>cd,intra</sub>, dB \n');
    if WideTable
        fprintf('!''''k''''<sub>sd, to</sub>, dB \n');
        fprintf('!''''k''''<sub>sd, from</sub>, dB \n');
        fprintf('!''''MP''''<sub>err</sub>, m \n');        
        fprintf('!''''T''''<sub>acq</sub>, s \n');
    end
    fprintf('!scores');
    
    for Signal_Type = 3:3


        if Access_Type == OpenSignal
            max_HalfBand = 1 + 1 + 100;
            if Signal_Type == BOCsin
                fid_BOCsin = fopen('stuff/complexity_stuff/L2OC_BOCsin.txt', 'w');
            elseif Signal_Type == BOCcos
                fid_BOCcos = fopen('stuff/complexity_stuff/L2OC_BOCcos.txt', 'w');
            elseif Signal_Type == BPSK
                fid_BPSK = fopen('stuff/complexity_stuff/L2OC_BPSK.txt', 'w');
            end
        elseif Access_Type == CloseSignal
            max_HalfBand = 5 + 2.5;
            if Signal_Type == BOCsin
                fid_BOCsin = fopen('stuff/complexity_stuff/L2SC_BOCsin.txt', 'w');
            elseif Signal_Type == BOCcos
                fid_BOCcos = fopen('stuff/complexity_stuff/L2SC_BOCcos.txt', 'w');
            elseif Signal_Type == BPSK
                fid_BPSK = fopen('stuff/complexity_stuff/L2SC_BPSK.txt', 'w');
            end
        end

        Signals_BOCsin_L2OC = [NaN NaN NaN];
        Signals_BOCcos_L2OC = [NaN NaN NaN];
        Signals_BPSK_L2OC = [NaN NaN];
        Signals_BOCsin_L2SC = [NaN NaN NaN];
        Signals_BOCcos_L2SC = [NaN NaN NaN];
        Signals_BPSK_L2SC = [NaN NaN];


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
                        if Pomexoyst_BoCsin(m8, n8) <= Threshold_Pomex
                            if (1/Beta_BoCsin(m8, n8)) <= Threshold_Accur
                                if InSysJam_BoCsin(m8, n8) <= Threshold_Intra
                                    if Complexity_BoCsin(m8, n8, freq_index) <= Threshold_Complexity
                                        if ((m+n + farr(freq_index)) <= (farr(end) + 1)) && ((farr(freq_index) - (m+n)) >= (farr(1) - 1))...
                                                && ((m+n) <= max_HalfBand)
                                            if (Access_Type == OpenSignal)||(n == 2.5)||(n == 2.625)%||(n == 2.25)
                                                fprintf('\n|- align="center"\n');
                                                fprintf('| BOC<sub>sin</sub>');
                                                fprintf('|| %.3f ', m8/8);
                                                fprintf('|| %.3f ', n8/8);
                                                fprintf('|| %4.0f ', farr(freq_index));
                                                fprintf('|| %.1f ', Pomexoyst_BoCsin(m8, n8));
                                                fprintf('|| %.2f ', 1/Beta_BoCsin(m8, n8)*1e9);
                                                fprintf('|| %.1f ', InSysJam_BoCsin(m8, n8));
                                                if WideTable
                                                    fprintf('|| %.1f ', BackInterSysJam_BoCsin_L2(m8, n8, freq_index));
                                                    fprintf('|| %.1f ', InterSysJam_BoCsin_L2(m8, n8, freq_index));
                                                    fprintf('|| %.1f ', MP_BOCsin(m8, n8));
                                                    fprintf('|| %.0f ', AcqBOC(n8));
                                                end
                                                fprintf('|| %.1f ', Complexity_BoCsin(m8, n8, freq_index));
                                                counter = counter + 1;
                                                fprintf(fid_BOCsin, '%.3f %.3f \n', m8/8, n8/8);
                                                if Access_Type == OpenSignal
                                                    Signals_BOCsin_L2OC(counter, 1) = m;
                                                    Signals_BOCsin_L2OC(counter, 2) = n;
                                                    Signals_BOCsin_L2OC(counter, 3) = farr(freq_index);
                                                elseif Access_Type == CloseSignal
                                                    Signals_BOCsin_L2SC(counter, 1) = m;
                                                    Signals_BOCsin_L2SC(counter, 2) = n;
                                                    Signals_BOCsin_L2SC(counter, 3) = farr(freq_index);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    elseif Signal_Type == BOCcos
                        if Pomexoyst_BoCcos(m8, n8) <= Threshold_Pomex
                            if (1/Beta_BoCcos(m8, n8)) <= Threshold_Accur
                                if InSysJam_BoCcos(m8, n8) <= Threshold_Intra
                                    if Complexity_BoCcos(m8, n8, freq_index) <= Threshold_Complexity
                                        if ((m+n + farr(freq_index)) <= (farr(end) + 1)) && ((farr(freq_index) - (m+n)) >= (farr(1) - 1))...
                                                && ((m+n) <= max_HalfBand)
                                            if (Access_Type == OpenSignal)||(n == 2.5)||(n == 2.625)%||(n == 2.25)
                                                fprintf('\n|- align="center"\n');
                                                fprintf('| BOC<sub>cos</sub>');
                                                fprintf('|| %.3f ', m8/8);
                                                fprintf('|| %.3f ', n8/8);
                                                fprintf('|| %4.0f ', farr(freq_index));
                                                fprintf('|| %.1f ', Pomexoyst_BoCcos(m8, n8));
                                                fprintf('|| %.2f ', 1/Beta_BoCcos(m8, n8)*1e9);
                                                fprintf('|| %.1f ', InSysJam_BoCcos(m8, n8));
                                                if WideTable
                                                    fprintf('|| %.1f ', BackInterSysJam_BoCcos_L2(m8, n8, freq_index));
                                                    fprintf('|| %.1f ', InterSysJam_BoCcos_L2(m8, n8, freq_index));
                                                    fprintf('|| %.1f ', MP_BOCcos(m8, n8));
                                                    fprintf('|| %.0f ', AcqBOC(n8));
                                                end                                                
                                                fprintf('|| %.1f ', Complexity_BoCcos(m8, n8, freq_index));
                                                counter = counter + 1;
                                                fprintf(fid_BOCcos, '%.3f %.3f \n', m8/8, n8/8);
                                                if Access_Type == OpenSignal
                                                    Signals_BOCcos_L2OC(counter, 1) = m;
                                                    Signals_BOCcos_L2OC(counter, 2) = n;
                                                    Signals_BOCcos_L2OC(counter, 3) = farr(freq_index);
                                                elseif Access_Type == CloseSignal
                                                    Signals_BOCcos_L2SC(counter, 1) = m;
                                                    Signals_BOCcos_L2SC(counter, 2) = n;
                                                    Signals_BOCcos_L2SC(counter, 3) = farr(freq_index);
                                                end                                        
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    elseif Signal_Type == BPSK
                        if Pomexoyst_BPSK(n8) <= Threshold_Pomex
                            if (1/Beta_BPSK(n8)) <= Threshold_Accur
                                if InSysJam_BPSK(n8) <= Threshold_Intra
                                    if Complexity_BPSK(n8, freq_index) <= Threshold_Complexity
                                        if ((m+n + farr(freq_index)) <= (farr(end) + 1)) && ((farr(freq_index) - (m+n)) >= (farr(1) - 1))...
                                                && ((n) <= max_HalfBand)
                                            if (Access_Type == OpenSignal)||(n == 2.5)||(n == 2.625)%||(n == 2.25)
                                                fprintf('\n|- align="center"\n');
                                                fprintf('| BPSK');
                                                fprintf('|| - ');
                                                fprintf('|| %.3f ', n8/8);
                                                fprintf('|| %4.0f ', farr(freq_index));
                                                fprintf('|| %.1f ', Pomexoyst_BPSK(n8));
                                                fprintf('|| %.2f ', 1/Beta_BPSK(n8)*1e9);
                                                fprintf('|| %.1f ', InSysJam_BPSK(n8));
                                                if WideTable
                                                    fprintf('|| %.1f ', BackInterSysJam_BPSK_L2(n8, freq_index));
                                                    fprintf('|| %.1f ', InterSysJam_BPSK_L2(n8, freq_index));
                                                    fprintf('|| %.1f ', MP_BPSK(n8));
                                                    fprintf('|| %.0f ', AcqBPSK(n8));
                                                end                                                
                                                fprintf('|| %.1f ', Complexity_BPSK(n8, freq_index));
                                                counter = counter + 1;
                                                fprintf(fid_BPSK, '%.3f \n', n8/8);
                                                if Access_Type == OpenSignal
                                                    Signals_BPSK_L2OC(counter, 1) = m;
                                                    Signals_BPSK_L2OC(counter, 2) = n;
                                                    Signals_BPSK_L2OC(counter, 3) = farr(freq_index);
                                                elseif Access_Type == CloseSignal
                                                    Signals_BPSK_L2SC(counter, 1) = m;
                                                    Signals_BPSK_L2SC(counter, 2) = n;
                                                    Signals_BPSK_L2SC(counter, 3) = farr(freq_index);
                                                end                                        
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
        
        if Signal_Type == BPSK
            fclose(fid_BPSK);
        elseif Signal_Type == BOCsin
            fclose(fid_BOCsin);
        elseif Signal_Type == BOCcos
            fclose(fid_BOCcos);
        end

        if Signal_Type == BOCsin
            if Access_Type == OpenSignal
                save('results/clipping/Signals_BOCsin_L2OC.mat', 'Signals_BOCsin_L2OC');
            elseif Access_Type == CloseSignal
                save('results/clipping/Signals_BOCsin_L2SC.mat', 'Signals_BOCsin_L2SC');
            end
        elseif Signal_Type == BOCcos
            if Access_Type == OpenSignal
                save('results/clipping/Signals_BOCcos_L2OC.mat', 'Signals_BOCcos_L2OC');
            elseif Access_Type == CloseSignal
                save('results/clipping/Signals_BOCcos_L2SC.mat', 'Signals_BOCcos_L2SC');
            end
        elseif Signal_Type == BPSK
            if Access_Type == OpenSignal
                save('results/clipping/Signals_BPSK_L2OC.mat', 'Signals_BPSK_L2OC')
            elseif Access_Type == CloseSignal
                save('results/clipping/Signals_BPSK_L2SC.mat', 'Signals_BPSK_L2SC')
            end
        end

        counter_all = counter_all + counter;    
    end 
    fprintf('\n|}\n');
end

fprintf('Counter all насчитал %.0f сигналов\n', counter_all);