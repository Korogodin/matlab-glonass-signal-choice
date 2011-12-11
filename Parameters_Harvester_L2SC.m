%/**
% Вспомогательный скрипт для проведения многокритериальной оптимизации
%*/

clear 
close all
clc

path_to_results_dir = [pwd '/results'];
path_to_prop = [pwd '/results/prop_L2'];

% Усеченное множество сигналов L2SC
load([path_to_results_dir '/clipping/Signals_BOCsin_L2SC.mat'], 'Signals_BOCsin_L2SC');
load([path_to_results_dir '/clipping/Signals_BOCcos_L2SC.mat'], 'Signals_BOCcos_L2SC');
load([path_to_results_dir '/clipping/Signals_BPSK_L2SC.mat'], 'Signals_BPSK_L2SC')
size_BOCsin = size(Signals_BOCsin_L2SC, 1)*(~isnan(Signals_BOCsin_L2SC(1,1)));
size_BOCcos = size(Signals_BOCcos_L2SC, 1)*(~isnan(Signals_BOCcos_L2SC(1,1)));
size_BPSK = size(Signals_BPSK_L2SC, 1)*(~isnan(Signals_BPSK_L2SC(1,1)));

n8max = 80;
m8max = 80;
farr = 1210:1228; fmax = length(farr); % Нормированный центральные частоты
qual_max = 12; % Сколько всего может быть ПК (с запасом взял)

ntmin = 0; % Интервал, к которому приводятся параметры
ntmax = 1;

% С помощью данной секции можно создать пустые массивы
Qual_BoCsin_L2SC = nan(m8max, n8max, fmax, qual_max);
save([path_to_prop '/Qual_BoCsin_L2SC.mat'], 'Qual_BoCsin_L2SC');
Qual_BoCcos_L2SC = nan(m8max, n8max, fmax, qual_max);
save([path_to_prop '/Qual_BoCcos_L2SC.mat'], 'Qual_BoCcos_L2SC');
Qual_BPSK_L2SC = nan(n8max, fmax, qual_max);
save([path_to_prop '/Qual_BPSK_L2SC.mat'], 'Qual_BPSK_L2SC');

constants;

load([path_to_prop '/Qual_BoCsin_L2SC.mat'], 'Qual_BoCsin_L2SC');
load([path_to_prop '/Qual_BoCcos_L2SC.mat'], 'Qual_BoCcos_L2SC');
load([path_to_prop '/Qual_BPSK_L2SC.mat'], 'Qual_BPSK_L2SC');

BOCsin = 1; BOCcos = 2; BPSK = 3;

% %%%%%%%%%%%%%%%%%%%% Search %%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_Evgeniy = [pwd '/results/Evgeniy'];
load([path_to_Evgeniy '/AcqBOC.mat'], 'AcqBOC');
load([path_to_Evgeniy '/AcqBPSK.mat'], 'AcqBPSK');

min_Search = 0;

max_Search = - 99999;
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    if max_Search < AcqBOC(n8)
        max_Search = AcqBOC(n8);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        if max_Search < AcqBOC(n8)
            max_Search = AcqBOC(n8);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                if max_Search < AcqBPSK(n8)
                    max_Search = AcqBPSK(n8);
                end
            end
        end

for i = 1:size_BOCsin
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    resrec = recalc_threshold(AcqBOC(n8), ntmin, ntmax, min_Search, max_Search);
    for freq_index = 1:fmax
        Qual_BoCsin_L2SC(m8, n8, freq_index, Search) = resrec;
    end
end
    for i = 1:size_BOCcos
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        resrec = recalc_threshold(AcqBOC(n8), ntmin, ntmax, min_Search, max_Search);
        for freq_index = 1:fmax
            Qual_BoCcos_L2SC(m8, n8, freq_index, Search) = resrec;
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                resrec = recalc_threshold(AcqBPSK(n8), ntmin, ntmax, min_Search, max_Search);
                for freq_index = 1:fmax
                    Qual_BPSK_L2SC(n8, freq_index, Search) = resrec;
                end
            end
        end
        

% % % %%%%%%%%%%%%%%%%%% MultiPath %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m8max = 80; n8max = 80;
% path_to_Evgeniy = [pwd '/results/Evgeniy'];
% load([path_to_Evgeniy '/BOCsin_Evg_L2.mat'], 'BOCsin_Evg_L2');
% load([path_to_Evgeniy '/BOCcos_Evg_L2.mat'], 'BOCcos_Evg_L2');
% 
% MP_BOCsin = nan(m8max, n8max);
% MP_BOCcos = nan(m8max, n8max);
% MP_BPSK = nan(1, n8max);
% 
% for i = 1:size(BOCsin_Evg_L2, 1)
%     m = BOCsin_Evg_L2(i, 1); m8 = m*8;
%     n = BOCsin_Evg_L2(i, 2); n8 = n*8;
%     MP_BOCsin(m8, n8) = BOCsin_Evg_L2(i, 4);
% end
% for i = 1:size(BOCcos_Evg_L2, 1)
%     m = BOCcos_Evg_L2(i, 1); m8 = m*8;
%     n = BOCcos_Evg_L2(i, 2); n8 = n*8;
%     MP_BOCcos(m8, n8) = BOCcos_Evg_L2(i, 4);
% end
% 
% save([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
% save([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
% save([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');


% %%%%%%%%%%%%%%%%%%%% MultiPath %%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_Evgeniy = [pwd '/results/Evgeniy'];
load([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
load([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
load([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');

MP_BOCsin = zeros(m8max, n8max);
MP_BOCcos = zeros(m8max, n8max);
MP_BPSK = zeros(1, n8max);

min_MP = 0;

max_MP = - 99999;
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    if max_MP < MP_BOCsin(m8, n8)
        max_MP = MP_BOCsin(m8, n8);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        if max_MP < MP_BOCcos(m8, n8)
            max_MP = MP_BOCcos(m8, n8);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                if max_MP < MP_BPSK(n8)
                    max_MP = MP_BPSK(n8);
                end
            end
        end

for i = 1:size_BOCsin
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    resrec = recalc_threshold(MP_BOCsin(m8, n8), ntmin, ntmax, min_MP, max_MP);
    for freq_index = 1:fmax
        Qual_BoCsin_L2SC(m8, n8, freq_index, MP) = resrec;
    end
end
    for i = 1:size_BOCcos
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        resrec = recalc_threshold(MP_BOCcos(m8, n8), ntmin, ntmax, min_MP, max_MP);
        for freq_index = 1:fmax
            Qual_BoCcos_L2SC(m8, n8, freq_index, MP) = resrec;
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                resrec = recalc_threshold(MP_BPSK(n8), ntmin, ntmax, min_MP, max_MP);
                for freq_index = 1:fmax
                    Qual_BPSK_L2SC(n8, freq_index, MP) = resrec;
                end
            end
        end
        

% %%%%%%%%%%%%%%%%%%%%%%%% Intrasystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_intra = [pwd '/results/intrasystem'];
load([path_to_results_intra '/InSysJam_BoCsin.mat']);
load([path_to_results_intra '/InSysJam_BoCcos.mat']);
load([path_to_results_intra '/InSysJam_BPSK.mat']);

max_Intra = - 99999;
min_Intra = + 99999;
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    if max_Intra < InSysJam_BoCsin(m8, n8)
        max_Intra = InSysJam_BoCsin(m8, n8);
    end
    if min_Intra > InSysJam_BoCsin(m8, n8)
        min_Intra = InSysJam_BoCsin(m8, n8);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        if max_Intra < InSysJam_BoCcos(m8, n8)
            max_Intra = InSysJam_BoCcos(m8, n8);
        end
        if min_Intra > InSysJam_BoCcos(m8, n8)
            min_Intra = InSysJam_BoCcos(m8, n8);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = 0; m8 = m*8;    
                if max_Intra < InSysJam_BPSK(n8)
                    max_Intra = InSysJam_BPSK(n8);
                end
                if min_Intra > InSysJam_BPSK(n8)
                    min_Intra = InSysJam_BPSK(n8);
                end
            end
        end
max_Intra = InSysJam_BPSK(1);

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    resrec = recalc_threshold(InSysJam_BoCsin(m8, n8), ntmin, ntmax, min_Intra, max_Intra);
    for freq_index = 1:fmax
        Qual_BoCsin_L2SC(m8, n8, freq_index, IntraS) = resrec;
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        resrec = recalc_threshold(InSysJam_BoCcos(m8, n8), ntmin, ntmax, min_Intra, max_Intra);
        for freq_index = 1:fmax
            Qual_BoCcos_L2SC(m8, n8, freq_index, IntraS) = resrec;
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                resrec = recalc_threshold(InSysJam_BPSK(n8), ntmin, ntmax, min_Intra, max_Intra);
                for freq_index = 1:fmax
                    Qual_BPSK_L2SC(n8, freq_index, IntraS) = resrec;
                end
            end
        end

    
% %%%%%%%%%%%%%%%%%%%%%%%% Complexity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_results_comp = [pwd '/results/complexity_L2'];

load([path_to_results_comp '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
load([path_to_results_comp '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
load([path_to_results_comp '/Complexity_BPSK.mat'], 'Complexity_BPSK');

max_Comp = - 99999;
min_Comp = 0;

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    if max_Comp < Complexity_BoCsin(m8, n8, freq_index)
        max_Comp = Complexity_BoCsin(m8, n8, freq_index);
    end
    if min_Comp > Complexity_BoCsin(m8, n8, freq_index)
        min_Comp = Complexity_BoCsin(m8, n8, freq_index);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;    
        if max_Comp < Complexity_BoCcos(m8, n8, freq_index)
            max_Comp = Complexity_BoCcos(m8, n8, freq_index);
        end
        if min_Comp > Complexity_BoCcos(m8, n8, freq_index)
            min_Comp = Complexity_BoCcos(m8, n8, freq_index);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = 0; m8 = m*8;    
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;    
                if max_Comp < Complexity_BPSK(n8, freq_index)
                    max_Comp = Complexity_BPSK(n8, freq_index);
                end
                if min_Comp > Complexity_BPSK(n8, freq_index)
                    min_Comp = Complexity_BPSK(n8, freq_index);
                end
            end
        end

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    resrec = recalc_threshold(Complexity_BoCsin(m8, n8, freq_index), ntmin, ntmax, min_Comp, max_Comp);
    Qual_BoCsin_L2SC(m8, n8, freq_index, Compl) = resrec;
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;           
        resrec = recalc_threshold(Complexity_BoCcos(m8, n8, freq_index), ntmin, ntmax, min_Comp, max_Comp);
        Qual_BoCcos_L2SC(m8, n8, freq_index, Compl) = resrec;
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;               
                resrec = recalc_threshold(Complexity_BPSK(n8, freq_index), ntmin, ntmax, min_Comp, max_Comp);
                Qual_BPSK_L2SC(n8, freq_index, Compl) = resrec;
            end
        end


% %%%%%%%%%%%%%%%%%%%%%%%% Intersystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_results_inter_C = [pwd '/results/intersystem_L2/Common'];

load([path_to_results_inter_C '/InterSysJam_BoCsin_L2.mat'], 'InterSysJam_BoCsin_L2');
load([path_to_results_inter_C '/InterSysJam_BoCcos_L2.mat'], 'InterSysJam_BoCcos_L2');
load([path_to_results_inter_C '/InterSysJam_BPSK_L2.mat'], 'InterSysJam_BPSK_L2');

max_Inter = - 99999;
min_Inter = + 99999;

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    if max_Inter < InterSysJam_BoCsin_L2(m8, n8, freq_index)
        max_Inter = InterSysJam_BoCsin_L2(m8, n8, freq_index);
    end
    if min_Inter > InterSysJam_BoCsin_L2(m8, n8, freq_index)
        min_Inter = InterSysJam_BoCsin_L2(m8, n8, freq_index);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;    
        if max_Inter < InterSysJam_BoCcos_L2(m8, n8, freq_index)
            max_Inter = InterSysJam_BoCcos_L2(m8, n8, freq_index);
        end
        if min_Inter > InterSysJam_BoCcos_L2(m8, n8, freq_index)
            min_Inter = InterSysJam_BoCcos_L2(m8, n8, freq_index);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = 0; m8 = m*8;    
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;    
                if max_Inter < InterSysJam_BPSK_L2(n8, freq_index)
                    max_Inter = InterSysJam_BPSK_L2(n8, freq_index);
                end
                if min_Inter > InterSysJam_BPSK_L2(n8, freq_index)
                    min_Inter = InterSysJam_BPSK_L2(n8, freq_index);
                end
            end
        end
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    resrec = recalc_threshold(InterSysJam_BoCsin_L2(m8, n8, freq_index), ntmin, ntmax, min_Inter, max_Inter);
    Qual_BoCsin_L2SC(m8, n8, freq_index, InterNam) = resrec;
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;           
        resrec = recalc_threshold(InterSysJam_BoCcos_L2(m8, n8, freq_index), ntmin, ntmax, min_Inter, max_Inter);
        Qual_BoCcos_L2SC(m8, n8, freq_index, InterNam) = resrec;
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;               
                resrec = recalc_threshold(InterSysJam_BPSK_L2(n8, freq_index), ntmin, ntmax, min_Inter, max_Inter);
                Qual_BPSK_L2SC(n8, freq_index, InterNam) = resrec;
            end        
        end
        


% %%%%%%%%%%%%%%%%%%%%%%%% BackIntersystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_results_inter_C = [pwd '/results/back_intersystem_L2/Common'];

load([path_to_results_inter_C '/BackInterSysJam_BoCsin_L2.mat'], 'BackInterSysJam_BoCsin_L2');
load([path_to_results_inter_C '/BackInterSysJam_BoCcos_L2.mat'], 'BackInterSysJam_BoCcos_L2');
load([path_to_results_inter_C '/BackInterSysJam_BPSK_L2.mat'], 'BackInterSysJam_BPSK_L2');

max_BackInter = - 99999;
min_BackInter = + 99999;

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    if max_BackInter < BackInterSysJam_BoCsin_L2(m8, n8, freq_index)
        max_BackInter = BackInterSysJam_BoCsin_L2(m8, n8, freq_index);
    end
    if min_BackInter > BackInterSysJam_BoCsin_L2(m8, n8, freq_index)
        min_BackInter = BackInterSysJam_BoCsin_L2(m8, n8, freq_index);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;    
        if max_BackInter < BackInterSysJam_BoCcos_L2(m8, n8, freq_index)
            max_BackInter = BackInterSysJam_BoCcos_L2(m8, n8, freq_index);
        end
        if min_BackInter > BackInterSysJam_BoCcos_L2(m8, n8, freq_index)
            min_BackInter = BackInterSysJam_BoCcos_L2(m8, n8, freq_index);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = 0; m8 = m*8;    
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;    
                if max_BackInter < BackInterSysJam_BPSK_L2(n8, freq_index)
                    max_BackInter = BackInterSysJam_BPSK_L2(n8, freq_index);
                end
                if min_BackInter > BackInterSysJam_BPSK_L2(n8, freq_index)
                    min_BackInter = BackInterSysJam_BPSK_L2(n8, freq_index);
                end
            end
        end
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    resrec = recalc_threshold(BackInterSysJam_BoCsin_L2(m8, n8, freq_index), ntmin, ntmax, min_BackInter, max_BackInter);
    Qual_BoCsin_L2SC(m8, n8, freq_index, InterIm) = resrec;
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;           
        resrec = recalc_threshold(BackInterSysJam_BoCcos_L2(m8, n8, freq_index), ntmin, ntmax, min_BackInter, max_BackInter);
        Qual_BoCcos_L2SC(m8, n8, freq_index, InterIm) = resrec;
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;               
                resrec = recalc_threshold(BackInterSysJam_BPSK_L2(n8, freq_index), ntmin, ntmax, min_BackInter, max_BackInter);
                Qual_BPSK_L2SC(n8, freq_index, InterIm) = resrec;
            end   
        end
        
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%% RA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_RA = [pwd '/results/radioastronomy'];
load([path_to_results_RA '/Radioastronomy_BoCsin.mat'], 'Radioastronomy_BoCsin');
load([path_to_results_RA '/Radioastronomy_BoCcos.mat'], 'Radioastronomy_BoCcos');
load([path_to_results_RA '/Radioastronomy_BPSK.mat'], 'Radioastronomy_BPSK');

max_RA = - 99999;
min_RA = + 99999;

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    if max_RA < Radioastronomy_BoCsin(m8, n8, freq_index)
        max_RA = Radioastronomy_BoCsin(m8, n8, freq_index);
    end
    if min_RA > Radioastronomy_BoCsin(m8, n8, freq_index)
        min_RA = Radioastronomy_BoCsin(m8, n8, freq_index);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;    
        if max_RA < Radioastronomy_BoCcos(m8, n8, freq_index)
            max_RA = Radioastronomy_BoCcos(m8, n8, freq_index);
        end
        if min_RA > Radioastronomy_BoCcos(m8, n8, freq_index)
            min_RA = Radioastronomy_BoCcos(m8, n8, freq_index);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = Signals_BPSK_L2SC(i, 1); m8 = m*8; % = 0
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;    
                if max_RA < Radioastronomy_BPSK(n8, freq_index)
                    max_RA = Radioastronomy_PSK(n8, freq_index);
                end
                if min_RA > Radioastronomy_BPSK(n8, freq_index)
                    min_RA = Radioastronomy_BPSK(n8, freq_index);
                end
            end    
        end

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;  
    resrec = recalc_threshold(Radioastronomy_BoCsin(m8, n8, freq_index), ntmin, ntmax, min_RA, max_RA);
    Qual_BoCsin_L2SC(m8, n8, freq_index, RA) = resrec;
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;  
        resrec = recalc_threshold(Radioastronomy_BoCcos(m8, n8, freq_index), ntmin, ntmax, min_RA, max_RA);
        Qual_BoCcos_L2SC(m8, n8, freq_index, RA) = resrec;
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = Signals_BPSK_L2SC(i, 1); m8 = m*8; % = 0
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;  
                resrec = recalc_threshold(Radioastronomy_BPSK(n8, freq_index), ntmin, ntmax, min_RA, max_RA);
                Qual_BPSK_L2SC(n8, freq_index, RA) = resrec;
            end
        end
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%% Accur %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_accur = [pwd '/results/beta'];

load([path_to_results_accur '/Beta_BoCsin.mat'], 'Beta_BoCsin');
load([path_to_results_accur '/Beta_BoCcos.mat'], 'Beta_BoCcos');
load([path_to_results_accur '/Beta_BPSK.mat'], 'Beta_BPSK');

max_Accur = - 99999;
% min_Accur = + 99999;
min_Accur = 0;

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    if max_Accur < (1 / Beta_BoCsin(m8, n8))
        max_Accur = 1 / Beta_BoCsin(m8, n8);
    end
    if min_Accur > (1 / Beta_BoCsin(m8, n8))
        min_Accur = 1 / Beta_BoCsin(m8, n8);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        if max_Accur < (1 / Beta_BoCcos(m8, n8))
            max_Accur = 1 / Beta_BoCcos(m8, n8);
        end
        if min_Accur > (1 / Beta_BoCcos(m8, n8))
            min_Accur = 1 / Beta_BoCcos(m8, n8);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))    
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = 0; m8 = m*8;    
                if max_Accur < (1 / Beta_BPSK(n8))
                    max_Accur = 1 / Beta_BPSK(n8);
                end
                if min_Accur > (1 / Beta_BPSK(n8))
                    min_Accur = 1 / Beta_BPSK(n8);
                end
            end
        end
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    resrec = recalc_threshold(1 / Beta_BoCsin(m8, n8), ntmin, ntmax, min_Accur, max_Accur);
    for freq_index = 1:fmax
        Qual_BoCsin_L2SC(m8, n8, freq_index, Accur) = resrec;
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        resrec = recalc_threshold(1 / Beta_BoCcos(m8, n8), ntmin, ntmax, min_Accur, max_Accur);
        for freq_index = 1:fmax
            Qual_BoCcos_L2SC(m8, n8, freq_index, Accur) = resrec;
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))    
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                resrec = recalc_threshold(1 / Beta_BPSK(n8), ntmin, ntmax, min_Accur, max_Accur);
                for freq_index = 1:fmax
                    Qual_BPSK_L2SC(n8, freq_index, Accur) = resrec;
                end
            end   
        end


% %%%%%%%%%%%%%%%%%%%%%%%%%%% Pomex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_pomex = [pwd '/results/pomexoyst'];

load([path_to_results_pomex '/Pomexoyst_BoCsin.mat'], 'Pomexoyst_BoCsin');
load([path_to_results_pomex '/Pomexoyst_BoCcos.mat'], 'Pomexoyst_BoCcos');
load([path_to_results_pomex '/Pomexoyst_BPSK.mat'], 'Pomexoyst_BPSK');

% Зачем я это делал оО
% Pomexoyst_BoCsin = 10.^(Pomexoyst_BoCsin/10);
% Pomexoyst_BoCcos = 10.^(Pomexoyst_BoCcos/10);
% Pomexoyst_BPSK = 10.^(Pomexoyst_BPSK/10);

max_Pomex = - 99999;
min_Pomex = + 99999;

for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    if max_Pomex < Pomexoyst_BoCsin(m8, n8)
        max_Pomex = Pomexoyst_BoCsin(m8, n8);
    end
    if min_Pomex > Pomexoyst_BoCsin(m8, n8)
        min_Pomex = Pomexoyst_BoCsin(m8, n8);
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        if max_Pomex < Pomexoyst_BoCcos(m8, n8)
            max_Pomex = Pomexoyst_BoCcos(m8, n8);
        end
        if min_Pomex > Pomexoyst_BoCcos(m8, n8)
            min_Pomex = Pomexoyst_BoCcos(m8, n8);
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                m = 0; m8 = m*8;    
                if max_Pomex < Pomexoyst_BPSK(n8)
                    max_Pomex = Pomexoyst_BPSK(n8);
                end
                if min_Pomex > Pomexoyst_BPSK(n8)
                    min_Pomex = Pomexoyst_BPSK(n8);
                end
            end
        end
for i = 1:size_BOCsin
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
    resrec = recalc_threshold(Pomexoyst_BoCsin(m8, n8), ntmin, ntmax, min_Pomex, max_Pomex);
    for freq_index = 1:fmax
        Qual_BoCsin_L2SC(m8, n8, freq_index, Pomex) = resrec;
    end
end
    for i = 1:size_BOCcos
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        resrec = recalc_threshold(Pomexoyst_BoCcos(m8, n8), ntmin, ntmax, min_Pomex, max_Pomex);
        for freq_index = 1:fmax
            Qual_BoCcos_L2SC(m8, n8, freq_index, Pomex) = resrec;
        end
    end
        if ~isnan(Signals_BPSK_L2SC(1,1))    
            for i = 1:size_BPSK
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                resrec = recalc_threshold(Pomexoyst_BPSK(n8), ntmin, ntmax, min_Pomex, max_Pomex);
                for freq_index = 1:fmax
                    Qual_BPSK_L2SC(n8, freq_index, Pomex) = resrec;
                end
            end   
        end
        
        
Table_Type = 2; % 1 - All; 2 - Top
% %%%%%%%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SC
C_Compl_SC = 6;
C_Search_SC = 6;
C_Accur_SC = 6;
% C_RA_SC = 2;
C_Pomex_SC = 8;
C_MP_SC = 2;
C_InterNam_SC = 2;
C_InterIm_SC = 1;
C_IntraS_SC = 10;   
        
TopSize = 30;
All_Signals_Size = size_BOCsin + size_BOCcos + size_BPSK;
All_Signals_Norm = nan(All_Signals_Size, 1);
All_Signals_index = 0;

    fprintf('{| class="wikitable sortable" border="1" \n');
    if Table_Type == 1
        fprintf('|+ Multi-objective optimization for L2SC signal\n');
    elseif Table_Type == 2
        if TopSize < All_Signals_Size
            fprintf('|+ Multi-objective optimization for L2SC signal (Top %.0f)\n', TopSize);
        else
            TopSize = All_Signals_Size;
            fprintf('|+ Multi-objective optimization for L2SC signal\n');
        end
    end
    
    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!Antijamming<br>w=%.2f \n', C_Pomex_SC);
%     fprintf('!Astronomy<br>w=%.2f  \n', C_RA_SC);
    fprintf('!Accurancy<br>w=%.2f  \n', C_Accur_SC);
    fprintf('!MultiPath<br>w=%.2f  \n', C_MP_SC);
    fprintf('!Inter from<br>w=%.2f \n', C_InterNam_SC);
    fprintf('!Inter to<br>w=%.2f \n', C_InterIm_SC);
    fprintf('!Search<br>w=%.2f\n', C_Search_SC);
    fprintf('!Complexity<br>w=%.2f\n', C_Compl_SC);
    fprintf('!Intra<br>w=%.2f\n', C_IntraS_SC);
    fprintf('!Norm');    
    

for i = 1:size_BOCsin
    norm_s = 0;
    n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
    m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;

    freq = Signals_BOCsin_L2SC(i, 3);
    freq_index = freq - farr(1) + 1;    
    
    if Table_Type == 1
        fprintf('\n|- align="center"\n');
        fprintf('|BoC<sub>sin</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Pomex));
%         fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, RA));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Accur));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, MP));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, InterNam));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, InterIm));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Search));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Compl));
        fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, IntraS));
    end
    
    norm_s = norm_s + (C_Pomex_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, Pomex))^2;
%     norm_s = norm_s + (C_RA_SC*Qual_BoCsin_L2SC(m8, n8, freq_index,
%     RA))^2;
    norm_s = norm_s + (C_Accur_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, Accur))^2;
    norm_s = norm_s + (C_MP_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, MP))^2;
    norm_s = norm_s + (C_InterNam_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, InterNam))^2;
    norm_s = norm_s + (C_InterIm_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, InterIm))^2;
    norm_s = norm_s + (C_Search_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, Search))^2;
    norm_s = norm_s + (C_Compl_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, Compl))^2;
    norm_s = norm_s + (C_IntraS_SC*Qual_BoCsin_L2SC(m8, n8, freq_index, IntraS))^2;
    
    if Table_Type == 1
        fprintf('|| %.4f ', sqrt(norm_s));
    end
    
    All_Signals_index = All_Signals_index + 1;
    All_Signals_Norm(All_Signals_index) = sqrt(norm_s);

    if (m == 5)&&(n == 2.5)&&(farr(freq_index) == 1220)
        our_inde = All_Signals_index;
%         sqrt(norm_s)
    end
end

    for i = 1:size_BOCcos
        norm_s = 0;
        n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
        m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
        freq = Signals_BOCcos_L2SC(i, 3);
        freq_index = freq - farr(1) + 1;    
        
        if Table_Type == 1
            fprintf('\n|- align="center"\n');
            fprintf('|BoC<sub>cos</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Pomex));
%             fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, RA));
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Accur));
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, MP));
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, InterNam));
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, InterIm));
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Search));
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Compl));        
            fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, IntraS));
        end
        
        norm_s = norm_s + (C_Pomex_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, Pomex))^2;
%         norm_s = norm_s + (C_RA_SC*Qual_BoCcos_L2SC(m8, n8, freq_index,
%         RA))^2;        
        norm_s = norm_s + (C_Accur_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, Accur))^2;
        norm_s = norm_s + (C_MP_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, MP))^2;
        norm_s = norm_s + (C_InterNam_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, InterNam))^2;
        norm_s = norm_s + (C_InterIm_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, InterIm))^2;
        norm_s = norm_s + (C_Search_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, Search))^2;
        norm_s = norm_s + (C_Compl_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, Compl))^2;
        norm_s = norm_s + (C_IntraS_SC*Qual_BoCcos_L2SC(m8, n8, freq_index, IntraS))^2;

        if Table_Type == 1
            fprintf('|| %.4f ', sqrt(norm_s));
        end
        
        All_Signals_index = All_Signals_index + 1;
        All_Signals_Norm(All_Signals_index) = sqrt(norm_s);        

    end
    
        if ~isnan(Signals_BPSK_L2SC(1,1))
            for i = 1:size_BPSK
                norm_s = 0;
                n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                freq = Signals_BPSK_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;    
                
                if Table_Type == 1
                    fprintf('\n|- align="center"\n');
                    fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Pomex));
%                     fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, RA));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Accur));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, MP));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, InterNam));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, InterIm));                
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Search));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Compl));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, IntraS));
                end
                
                norm_s = norm_s + (C_Pomex_SC*Qual_BPSK_L2SC(n8, freq_index, Pomex))^2;
%                 norm_s = norm_s + (C_RA_SC*Qual_BPSK_L2SC(n8, freq_index, RA))^2;
                norm_s = norm_s + (C_Accur_SC*Qual_BPSK_L2SC(n8, freq_index, Accur))^2;
                norm_s = norm_s + (C_MP_SC*Qual_BPSK_L2SC(n8, freq_index, MP))^2;
                norm_s = norm_s + (C_InterNam_SC*Qual_BPSK_L2SC(n8, freq_index, InterNam))^2;
                norm_s = norm_s + (C_InterIm_SC*Qual_BPSK_L2SC(n8, freq_index, InterIm))^2;
                norm_s = norm_s + (C_Search_SC*Qual_BPSK_L2SC(n8, freq_index, Search))^2;
                norm_s = norm_s + (C_Compl_SC*Qual_BPSK_L2SC(n8, freq_index, Compl))^2;
                norm_s = norm_s + (C_IntraS_SC*Qual_BPSK_L2SC(n8, freq_index, IntraS))^2;

                if Table_Type == 1
                    fprintf('|| %.4f ', sqrt(norm_s));
                end
                
                All_Signals_index = All_Signals_index + 1;
                All_Signals_Norm(All_Signals_index) = sqrt(norm_s);   
                
            end 
        end

    
if Table_Type == 2        
    [All_Signals_Norm_a All_Signals_Norm_b] = sort(All_Signals_Norm, 1);
    for j = 1:(TopSize+1)

        if j == (TopSize+1)
            All_Signals_Norm_b(j) = our_inde;
        end
        fprintf('\n|- align="center"\n');

        if All_Signals_Norm_b(j) <= size_BOCsin
            i = All_Signals_Norm_b(j);
            m = Signals_BOCsin_L2SC(i, 1); m8 = m*8;
            n = Signals_BOCsin_L2SC(i, 2); n8 = n*8;
            freq = Signals_BOCsin_L2SC(i, 3);
            freq_index = freq - farr(1) + 1;    
            fprintf('|BOC<sub>sin</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Pomex));
%             fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index,
%             RA));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Accur));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, MP));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, InterNam));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, InterIm));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Search));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, Compl));
            fprintf('|| %.3f ', Qual_BoCsin_L2SC(m8, n8, freq_index, IntraS));
            fprintf('|| %.4f ', All_Signals_Norm( All_Signals_Norm_b(j) ));
        elseif All_Signals_Norm_b(j) <= (size_BOCsin + size_BOCcos)
                i = All_Signals_Norm_b(j) -  size_BOCsin;
                n = Signals_BOCcos_L2SC(i, 2); n8 = n*8;
                m = Signals_BOCcos_L2SC(i, 1); m8 = m*8;
                freq = Signals_BOCcos_L2SC(i, 3);
                freq_index = freq - farr(1) + 1;    
                fprintf('|BoC<sub>cos</sub>(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Pomex));
%                 fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index,
%                 RA));
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Accur));
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, MP));
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, InterNam));
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, InterIm));
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Search));
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, Compl));        
                fprintf('|| %.3f ', Qual_BoCcos_L2SC(m8, n8, freq_index, IntraS));
                fprintf('|| %.4f ', All_Signals_Norm( All_Signals_Norm_b(j) ));
        else
                    i = All_Signals_Norm_b(j) -  size_BOCsin -  size_BOCcos;
                    n = Signals_BPSK_L2SC(i, 2); n8 = n*8;
                    freq = Signals_BPSK_L2SC(i, 3);
                    freq_index = freq - farr(1) + 1;    
                    fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Pomex));
%                     fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index,
%                     RA));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Accur));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, MP));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, InterNam));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, InterIm));                
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Search));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, Compl));
                    fprintf('|| %.3f ', Qual_BPSK_L2SC(n8, freq_index, IntraS));
                    fprintf('|| %.4f ', All_Signals_Norm( All_Signals_Norm_b(j) ));

        end
    end
end

fprintf('\n|} \n');
                    
% %%%%%%%%%%%%%%%%%%%%%%%% Saving of results %%%%%%%%%%%%%%%%%%%%%%%
% save([path_to_prop '/Qual_BoCsin_L2.mat'], 'Qual_BoCsin_L2');
% save([path_to_prop '/Qual_BoCcos_L2.mat'], 'Qual_BoCcos_L2');
% save([path_to_prop '/Qual_BPSK_L2.mat'], 'Qual_BPSK_L2');
