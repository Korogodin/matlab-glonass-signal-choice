%/**
% Вспомогательный скрипт для проведения многокритериальной оптимизации
%*/

clear 
close all
clc

Signal_OC = 0;

path_to_results_dir = [pwd '/results'];
path_to_prop = [pwd '/results/prop_L1'];

n8max = 80;
m8max = 80;
farr = 1558:1573; fmax = length(farr); % Нормированный центральные частоты
qual_max = 12;

ntmin = 0; % Интервал, к которому приводятся параметры
ntmax = 1;

% С помощью данной секции можно создать пустые массивы
Qual_BoCsin_L1 = nan(m8max, n8max, fmax, qual_max);
save([path_to_prop '/Qual_BoCsin_L1.mat'], 'Qual_BoCsin_L1');
Qual_BoCcos_L1 = nan(m8max, n8max, fmax, qual_max);
save([path_to_prop '/Qual_BoCcos_L1.mat'], 'Qual_BoCcos_L1');
Qual_BPSK_L1 = nan(n8max, fmax, qual_max);
save([path_to_prop '/Qual_BPSK_L1.mat'], 'Qual_BPSK_L1');

constants;

load([path_to_prop '/Qual_BoCsin_L1.mat'], 'Qual_BoCsin_L1');
load([path_to_prop '/Qual_BoCcos_L1.mat'], 'Qual_BoCcos_L1');
load([path_to_prop '/Qual_BPSK_L1.mat'], 'Qual_BPSK_L1');

Signals_L1;

BOCsin = 1; BOCcos = 2; BPSK = 3;

% %%%%%%%%%%%%%%%%%%%% Search %%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_Evgeniy = [pwd '/results/Evgeniy'];
load([path_to_Evgeniy '/BOCsin_Evg_L1.mat'], 'BOCsin_Evg_L1');
load([path_to_Evgeniy '/BOCcos_Evg_L1.mat'], 'BOCcos_Evg_L1');
% load([path_to_Evgeniy '/BPSK_Evg_L1.mat'], 'BPSK_Evg_L1');

max_Search = max([max(BOCsin_Evg_L1(:, 3));
                  max(BOCcos_Evg_L1(:, 3)); ]);
%                 max(BPSK_Evg_L1(:, 3))]);

min_Search = min([min(BOCsin_Evg_L1(:, 3));
                  min(BOCcos_Evg_L1(:, 3)); ]);
%                 min(BPSK_Evg_L1(:, 3))]);
min_Search = 0; % Я считаю, что так рациональнее и нагляднее, сохраняется линейность

for i = 1:size(BOCsin_Evg_L1, 1)
    m = BOCsin_Evg_L1(i, 1); m8 = m*8;
    n = BOCsin_Evg_L1(i, 2); n8 = n*8;
%     fprintf('%f, %f \n', m, n);
    resrec = recalc_threshold(BOCsin_Evg_L1(i, 3), ntmin, ntmax, min_Search, max_Search);
    for freq_index = 1:fmax
        Qual_BoCsin_L1(m8, n8, freq_index, Search) = resrec;
    end
%         fprintf('%f\n', resrec);
end
    for i = 1:size(BOCcos_Evg_L1, 1)
        m = BOCcos_Evg_L1(i, 1); m8 = m*8;
        n = BOCcos_Evg_L1(i, 2); n8 = n*8;
        resrec = recalc_threshold(BOCcos_Evg_L1(i, 3), ntmin, ntmax, min_Search, max_Search);
        for freq_index = 1:fmax
            Qual_BoCcos_L1(m8, n8, freq_index, Search) = resrec;
        end
    %         fprintf('%f\n', resrec);
    end
        % for i = 1:size(BOCcos_Evg_L1, 1)
        %     m = BPSK_Evg_L1(i, 1); m8 = m*8*0;
        %     n = BPSK_Evg_L1(i, 2); n8 = n*8;
        %     resrec = recalc_threshold(BPSK_Evg_L1(i, 3), ntmin, ntmax, min_Search, max_Search);
        %     for freq_index = 1:fmax
        %         Qual_BPSK_L1(n8, freq_index, Search) = resrec;
        %     end
        % %         fprintf('%f\n', resrec);
        % end



% %%%%%%%%%%%%%%%%%% MultiPath %%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_Evgeniy = [pwd '/results/Evgeniy'];
load([path_to_Evgeniy '/BOCsin_Evg_L1.mat'], 'BOCsin_Evg_L1');
load([path_to_Evgeniy '/BOCcos_Evg_L1.mat'], 'BOCcos_Evg_L1');
% load([path_to_Evgeniy '/BPSK_Evg_L1.mat'], 'BPSK_Evg_L1');

max_MP = max([max(BOCsin_Evg_L1(:, 4));
              max(BOCcos_Evg_L1(:, 4)); ]);
%             max(BPSK_Evg_L1(:, 4))]);

min_MP = min([min(BOCsin_Evg_L1(:, 4));
              min(BOCcos_Evg_L1(:, 4)); ]);
%             min(BPSK_Evg_L1(:, 4))]);
min_MP = 0; % Я считаю, что так рациональнее и нагляднее, сохраняется линейность

for i = 1:size(BOCsin_Evg_L1, 1)
    m = BOCsin_Evg_L1(i, 1); m8 = m*8;
    n = BOCsin_Evg_L1(i, 2); n8 = n*8;
    resrec = recalc_threshold(BOCsin_Evg_L1(i, 4), ntmin, ntmax, min_MP, max_MP);
    for freq_index = 1:fmax
        Qual_BoCsin_L1(m8, n8, freq_index, MP) = resrec;
    end
%         fprintf('%f\n', resrec);
end
for i = 1:size(BOCcos_Evg_L1, 1)
    m = BOCcos_Evg_L1(i, 1); m8 = m*8;
    n = BOCcos_Evg_L1(i, 2); n8 = n*8;
    resrec = recalc_threshold(BOCcos_Evg_L1(i, 4), ntmin, ntmax, min_MP, max_MP);
    for freq_index = 1:fmax
        Qual_BoCcos_L1(m8, n8, freq_index, MP) = resrec;
    end
%         fprintf('%f\n', resrec);
end
% for i = 1:size(BOCcos_Evg_L1, 1)
%     m = BPSK_Evg_L1(i, 1); m8 = m*8*0;
%     n = BPSK_Evg_L1(i, 2); n8 = n*8;
%     resrec = recalc_threshold(BPSK_Evg_L1(i, 4), ntmin, ntmax, min_MP, max_MP);
%     for freq_index = 1:fmax
%         Qual_BPSK_L1(n8, freq_index, MP) = resrec;
%     end
% %         fprintf('%f\n', resrec);
% end


% %%%%%%%%%%%%%%%%%%%%%%%% Intrasystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_intra = [pwd '/results/intrasystem'];
load([path_to_results_intra '/InSysJam_BoCsin.mat']);
load([path_to_results_intra '/InSysJam_BoCcos.mat']);
load([path_to_results_intra '/InSysJam_BPSK.mat']);

max_Intra = - 99999;
min_Intra = + 99999;
for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    if max_Intra < InSysJam_BoCsin(m8, n8)
        max_Intra = InSysJam_BoCsin(m8, n8);
    end
    if min_Intra > InSysJam_BoCsin(m8, n8)
        min_Intra = InSysJam_BoCsin(m8, n8);
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        if max_Intra < InSysJam_BoCcos(m8, n8)
            max_Intra = InSysJam_BoCcos(m8, n8);
        end
        if min_Intra > InSysJam_BoCcos(m8, n8)
            min_Intra = InSysJam_BoCcos(m8, n8);
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;    
            if max_Intra < InSysJam_BPSK(n8)
                max_Intra = InSysJam_BPSK(n8);
            end
            if min_Intra > InSysJam_BPSK(n8)
                min_Intra = InSysJam_BPSK(n8);
            end
        end

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    resrec = recalc_threshold(InSysJam_BoCsin(m8, n8), ntmin, ntmax, min_Intra, max_Intra);
    for freq_index = 1:fmax
        Qual_BoCsin_L1(m8, n8, freq_index, IntraS) = resrec;
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        resrec = recalc_threshold(InSysJam_BoCcos(m8, n8), ntmin, ntmax, min_Intra, max_Intra);
        for freq_index = 1:fmax
            Qual_BoCcos_L1(m8, n8, freq_index, IntraS) = resrec;
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            resrec = recalc_threshold(InSysJam_BPSK(n8), ntmin, ntmax, min_Intra, max_Intra);
            for freq_index = 1:fmax
                Qual_BPSK_L1(n8, freq_index, IntraS) = resrec;
            end
        end

    
% %%%%%%%%%%%%%%%%%%%%%%%% Complexity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_results_comp = [pwd '/results/complexity_L1'];

load([path_to_results_comp '/Complexity_BoCsin.mat'], 'Complexity_BoCsin');
load([path_to_results_comp '/Complexity_BoCcos.mat'], 'Complexity_BoCcos');
load([path_to_results_comp '/Complexity_BPSK.mat'], 'Complexity_BPSK');

max_Comp = - 99999;
min_Comp = + 99999;

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    if max_Comp < Complexity_BoCsin(m8, n8, freq_index)
        max_Comp = Complexity_BoCsin(m8, n8, freq_index);
    end
    if min_Comp > Complexity_BoCsin(m8, n8, freq_index)
        min_Comp = Complexity_BoCsin(m8, n8, freq_index);
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;    
        if max_Comp < Complexity_BoCcos(m8, n8, freq_index)
            max_Comp = Complexity_BoCcos(m8, n8, freq_index);
        end
        if min_Comp > Complexity_BoCcos(m8, n8, freq_index)
            min_Comp = Complexity_BoCcos(m8, n8, freq_index);
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;    
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;    
            if max_Comp < Complexity_BPSK(n8, freq_index)
                max_Comp = Complexity_BPSK(n8, freq_index);
            end
            if min_Comp > Complexity_BPSK(n8, freq_index)
                min_Comp = Complexity_BPSK(n8, freq_index);
            end
        end

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    resrec = recalc_threshold(Complexity_BoCsin(m8, n8), ntmin, ntmax, min_Comp, max_Comp);
    for freq_index = 1:fmax
        Qual_BoCsin_L1(m8, n8, freq_index, Compl) = resrec;
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;           
        resrec = recalc_threshold(Complexity_BoCcos(m8, n8), ntmin, ntmax, min_Comp, max_Comp);
        for freq_index = 1:fmax
            Qual_BoCcos_L1(m8, n8, freq_index, Compl) = resrec;
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;               
            resrec = recalc_threshold(Complexity_BPSK(n8), ntmin, ntmax, min_Comp, max_Comp);
            for freq_index = 1:fmax
                Qual_BPSK_L1(n8, freq_index, Compl) = resrec;
            end
        end


% %%%%%%%%%%%%%%%%%%%%%%%% Intersystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_results_inter_C = [pwd '/results/intersystem_L1/Common'];

load([path_to_results_inter_C '/InterSysJam_BoCsin_L1.mat'], 'InterSysJam_BoCsin_L1');
load([path_to_results_inter_C '/InterSysJam_BoCcos_L1.mat'], 'InterSysJam_BoCcos_L1');
load([path_to_results_inter_C '/InterSysJam_BPSK_L1.mat'], 'InterSysJam_BPSK_L1');

max_Inter = - 99999;
min_Inter = + 99999;

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    if max_Inter < InterSysJam_BoCsin_L1(m8, n8, freq_index)
        max_Inter = InterSysJam_BoCsin_L1(m8, n8, freq_index);
    end
    if min_Inter > InterSysJam_BoCsin_L1(m8, n8, freq_index)
        min_Inter = InterSysJam_BoCsin_L1(m8, n8, freq_index);
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;    
        if max_Inter < InterSysJam_BoCcos_L1(m8, n8, freq_index)
            max_Inter = InterSysJam_BoCcos_L1(m8, n8, freq_index);
        end
        if min_Inter > InterSysJam_BoCcos_L1(m8, n8, freq_index)
            min_Inter = InterSysJam_BoCcos_L1(m8, n8, freq_index);
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;    
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;    
            if max_Inter < InterSysJam_BPSK_L1(n8, freq_index)
                max_Inter = InterSysJam_BPSK_L1(n8, freq_index);
            end
            if min_Inter > InterSysJam_BPSK_L1(n8, freq_index)
                min_Inter = InterSysJam_BPSK_L1(n8, freq_index);
            end
        end
for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    resrec = recalc_threshold(InterSysJam_BoCsin_L1(m8, n8, freq_index), ntmin, ntmax, min_Inter, max_Inter);
    Qual_BoCsin_L1(m8, n8, freq_index, InterNam) = resrec;
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;           
        resrec = recalc_threshold(InterSysJam_BoCcos_L1(m8, n8, freq_index), ntmin, ntmax, min_Inter, max_Inter);
        Qual_BoCcos_L1(m8, n8, freq_index, InterNam) = resrec;
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;               
            resrec = recalc_threshold(InterSysJam_BPSK_L1(n8, freq_index), ntmin, ntmax, min_Inter, max_Inter);
            Qual_BPSK_L1(n8, freq_index, InterNam) = resrec;
        end        
        


% %%%%%%%%%%%%%%%%%%%%%%%% BackIntersystem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_results_inter_C = [pwd '/results/back_intersystem_L1/Common'];

load([path_to_results_inter_C '/BackInterSysJam_BoCsin_L1.mat'], 'BackInterSysJam_BoCsin_L1');
load([path_to_results_inter_C '/BackInterSysJam_BoCcos_L1.mat'], 'BackInterSysJam_BoCcos_L1');
load([path_to_results_inter_C '/BackInterSysJam_BPSK_L1.mat'], 'BackInterSysJam_BPSK_L1');

max_BackInter = - 99999;
min_BackInter = + 99999;

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    if max_BackInter < BackInterSysJam_BoCsin_L1(m8, n8, freq_index)
        max_BackInter = BackInterSysJam_BoCsin_L1(m8, n8, freq_index);
    end
    if min_BackInter > BackInterSysJam_BoCsin_L1(m8, n8, freq_index)
        min_BackInter = BackInterSysJam_BoCsin_L1(m8, n8, freq_index);
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;    
        if max_BackInter < BackInterSysJam_BoCcos_L1(m8, n8, freq_index)
            max_BackInter = BackInterSysJam_BoCcos_L1(m8, n8, freq_index);
        end
        if min_BackInter > BackInterSysJam_BoCcos_L1(m8, n8, freq_index)
            min_BackInter = BackInterSysJam_BoCcos_L1(m8, n8, freq_index);
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;    
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;    
            if max_BackInter < BackInterSysJam_BPSK_L1(n8, freq_index)
                max_BackInter = BackInterSysJam_BPSK_L1(n8, freq_index);
            end
            if min_BackInter > BackInterSysJam_BPSK_L1(n8, freq_index)
                min_BackInter = BackInterSysJam_BPSK_L1(n8, freq_index);
            end
        end
for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    resrec = recalc_threshold(BackInterSysJam_BoCsin_L1(m8, n8, freq_index), ntmin, ntmax, min_BackInter, max_BackInter);
    Qual_BoCsin_L1(m8, n8, freq_index, InterIm) = resrec;
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;           
        resrec = recalc_threshold(BackInterSysJam_BoCcos_L1(m8, n8, freq_index), ntmin, ntmax, min_BackInter, max_BackInter);
        Qual_BoCcos_L1(m8, n8, freq_index, InterIm) = resrec;
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;               
            resrec = recalc_threshold(BackInterSysJam_BPSK_L1(n8, freq_index), ntmin, ntmax, min_BackInter, max_BackInter);
            Qual_BPSK_L1(n8, freq_index, InterIm) = resrec;
        end   
        
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%% RA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_RA = [pwd '/results/RA_L1'];
load([path_to_RA '/BOCsin_RA_L1.mat'], 'BOCsin_RA_L1');
load([path_to_RA '/BOCcos_RA_L1.mat'], 'BOCcos_RA_L1');
load([path_to_RA '/BPSK_RA_L1.mat'], 'BPSK_RA_L1');

max_RA = max([max(BOCsin_RA_L1(:, 4));
                  max(BOCcos_RA_L1(:, 4));
                  max(BPSK_RA_L1(:, 3))]);

min_RA = min([min(BOCsin_RA_L1(:, 4));
                  min(BOCcos_RA_L1(:, 4));
                  min(BPSK_RA_L1(:, 3))]);

for i = 1:size(BOCsin_RA_L1, 1)
    m = BOCsin_RA_L1(i, 1); m8 = m*8;
    n = BOCsin_RA_L1(i, 2); n8 = n*8;
    freq_val = BOCsin_RA_L1(i, 3);
    freq_index = freq_val - 1558 + 1;
    resrec = recalc_threshold(BOCsin_RA_L1(i, 4), ntmin, ntmax, min_RA, max_RA);
    Qual_BoCsin_L1(m8, n8, freq_index, RA) = resrec;
end
    for i = 1:size(BOCcos_RA_L1, 1)
        m = BOCcos_RA_L1(i, 1); m8 = m*8;
        n = BOCcos_RA_L1(i, 2); n8 = n*8;
        freq_val = BOCcos_RA_L1(i, 3);
        freq_index = freq_val - 1558 + 1;
        resrec = recalc_threshold(BOCcos_RA_L1(i, 4), ntmin, ntmax, min_RA, max_RA);
        Qual_BoCcos_L1(m8, n8, freq_index, RA) = resrec;
    end
        for i = 1:size(BPSK_RA_L1, 1)
            n = BPSK_RA_L1(i, 1); n8 = n*8;
            freq_val = BPSK_RA_L1(i, 2);
            freq_index = freq_val - 1558 + 1;            
            resrec = recalc_threshold(BPSK_RA_L1(i, 3), ntmin, ntmax, min_RA, max_RA);
            Qual_BPSK_L1(n8, freq_index, RA) = resrec;
        end              

        
% %%%%%%%%%%%%%%%%%%%%%%%%%%% Accur_old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path_to_Accur = [pwd '/results/Accur_L1'];
% load([path_to_Accur '/BOCsin_Accur_L1.mat'], 'BOCsin_Accur_L1');
% load([path_to_Accur '/BOCcos_Accur_L1.mat'], 'BOCcos_Accur_L1');
% load([path_to_Accur '/BPSK_Accur_L1.mat'], 'BPSK_Accur_L1');
% 
% max_Accur = max([max(BOCsin_Accur_L1(:, 3));
%                  max(BOCcos_Accur_L1(:, 3));
%                  max(BPSK_Accur_L1(:, 2))]);
% 
% min_Accur = min([min(BOCsin_Accur_L1(:, 3));
%                  min(BOCcos_Accur_L1(:, 3));
%                  min(BPSK_Accur_L1(:, 2))]);
%               
% for i = 1:size(BOCsin_Accur_L1, 1)
%     m = BOCsin_Accur_L1(i, 1); m8 = m*8;
%     n = BOCsin_Accur_L1(i, 2); n8 = n*8;
%     resrec = recalc_threshold(BOCsin_Accur_L1(i, 3), ntmin, ntmax, min_Accur, max_Accur);
%     for freq_index = 1:fmax
%         Qual_BoCsin_L1(m8, n8, freq_index, Accur) = resrec;
%     end
% end
%     for i = 1:size(BOCcos_Accur_L1, 1)
%         m = BOCcos_Accur_L1(i, 1); m8 = m*8;
%         n = BOCcos_Accur_L1(i, 2); n8 = n*8;
%         resrec = recalc_threshold(BOCcos_Accur_L1(i, 3), ntmin, ntmax, min_Accur, max_Accur);
%         for freq_index = 1:fmax
%             Qual_BoCcos_L1(m8, n8, freq_index, Accur) = resrec;
%         end
%     end
%         for i = 1:size(BPSK_Accur_L1, 1)
%             n = BPSK_Accur_L1(i, 1); n8 = n*8;
%             resrec = recalc_threshold(BPSK_Accur_L1(i, 2), ntmin, ntmax, min_Accur, max_Accur);
%             for freq_index = 1:fmax
%                 Qual_BPSK_L1(n8, freq_index, Accur) = resrec;
%             end
%         end              
             
% %%%%%%%%%%%%%%%%%%%%%%%%%%% Accur_new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_accur = [pwd '/results/beta'];

load([path_to_results_accur '/Beta_BoCsin.mat'], 'Beta_BoCsin');
load([path_to_results_accur '/Beta_BoCcos.mat'], 'Beta_BoCcos');
load([path_to_results_accur '/Beta_BPSK.mat'], 'Beta_BPSK');

max_Accur = - 99999;
min_Accur = + 99999;

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    if max_Accur < (1 / Beta_BoCsin(m8, n8))
        max_Accur = 1 / Beta_BoCsin(m8, n8);
    end
    if min_Accur > (1 / Beta_BoCsin(m8, n8))
        min_Accur = 1 / Beta_BoCsin(m8, n8);
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        if max_Accur < (1 / Beta_BoCcos(m8, n8))
            max_Accur = 1 / Beta_BoCcos(m8, n8);
        end
        if min_Accur > (1 / Beta_BoCcos(m8, n8))
            min_Accur = 1 / Beta_BoCcos(m8, n8);
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;    
            if max_Accur < (1 / Beta_BPSK(n8))
                max_Accur = 1 / Beta_BPSK(n8);
            end
            if min_Accur > (1 / Beta_BPSK(n8))
                min_Accur = 1 / Beta_BPSK(n8);
            end
        end
for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    resrec = recalc_threshold(1 / Beta_BoCsin(m8, n8), ntmin, ntmax, min_Accur, max_Accur);
    for freq_index = 1:fmax
        Qual_BoCsin_L1(m8, n8, freq_index, Accur) = resrec;
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        resrec = recalc_threshold(1 / Beta_BoCcos(m8, n8), ntmin, ntmax, min_Accur, max_Accur);
        for freq_index = 1:fmax
            Qual_BoCcos_L1(m8, n8, freq_index, Accur) = resrec;
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            resrec = recalc_threshold(1 / Beta_BPSK(n8), ntmin, ntmax, min_Accur, max_Accur);
            for freq_index = 1:fmax
                Qual_BPSK_L1(n8, freq_index, Accur) = resrec;
            end
        end   

        
        
% % %%%%%%%%%%%%%%%%%%%%%%%%%%% Pomex_old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path_to_Pomex = [pwd '/results/Pomex_L1'];
% load([path_to_Pomex '/BOCsin_Pomex_L1.mat'], 'BOCsin_Pomex_L1');
% load([path_to_Pomex '/BOCcos_Pomex_L1.mat'], 'BOCcos_Pomex_L1');
% load([path_to_Pomex '/BPSK_Pomex_L1.mat'], 'BPSK_Pomex_L1');
% 
% max_Pomex = max([max(BOCsin_Pomex_L1(:, 3));
%                  max(BOCcos_Pomex_L1(:, 3));
%                  max(BPSK_Pomex_L1(:, 2))]);
% 
% min_Pomex = min([min(BOCsin_Pomex_L1(:, 3));
%                  min(BOCcos_Pomex_L1(:, 3));
%                  min(BPSK_Pomex_L1(:, 2))]);
% min_Pomex = 0;
%               
% for i = 1:size(BOCsin_Pomex_L1, 1)
%     m = BOCsin_Pomex_L1(i, 1); m8 = m*8;
%     n = BOCsin_Pomex_L1(i, 2); n8 = n*8;
%     resrec = recalc_threshold(BOCsin_Pomex_L1(i, 3), ntmin, ntmax, min_Pomex, max_Pomex);
%     for freq_index = 1:fmax
%         Qual_BoCsin_L1(m8, n8, freq_index, Pomex) = resrec;
%     end
% end
%     for i = 1:size(BOCcos_Pomex_L1, 1)
%         m = BOCcos_Pomex_L1(i, 1); m8 = m*8;
%         n = BOCcos_Pomex_L1(i, 2); n8 = n*8;
%         resrec = recalc_threshold(BOCcos_Pomex_L1(i, 3), ntmin, ntmax, min_Pomex, max_Pomex);
%         for freq_index = 1:fmax
%             Qual_BoCcos_L1(m8, n8, freq_index, Pomex) = resrec;
%         end
%     end
%         for i = 1:size(BPSK_Pomex_L1, 1)
%             n = BPSK_Pomex_L1(i, 1); n8 = n*8;
%             resrec = recalc_threshold(BPSK_Pomex_L1(i, 2), ntmin, ntmax, min_Pomex, max_Pomex);
%             for freq_index = 1:fmax
%                 Qual_BPSK_L1(n8, freq_index, Pomex) = resrec;
%             end
%         end       
        

% %%%%%%%%%%%%%%%%%%%%%%%%%%% Pomex_new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_results_pomex = [pwd '/results/pomexoyst'];

load([path_to_results_pomex '/Pomexoyst_BoCsin.mat'], 'Pomexoyst_BoCsin');
load([path_to_results_pomex '/Pomexoyst_BoCcos.mat'], 'Pomexoyst_BoCcos');
load([path_to_results_pomex '/Pomexoyst_BPSK.mat'], 'Pomexoyst_BPSK');

Pomexoyst_BoCsin = 10.^(Pomexoyst_BoCsin/10);
Pomexoyst_BoCcos = 10.^(Pomexoyst_BoCcos/10);
Pomexoyst_BPSK = 10.^(Pomexoyst_BPSK/10);

max_Pomex = - 99999;
min_Pomex = + 99999;

for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    if max_Pomex < Pomexoyst_BoCsin(m8, n8)
        max_Pomex = Pomexoyst_BoCsin(m8, n8);
    end
    if min_Pomex > Pomexoyst_BoCsin(m8, n8)
        min_Pomex = Pomexoyst_BoCsin(m8, n8);
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        if max_Pomex < Pomexoyst_BoCcos(m8, n8)
            max_Pomex = Pomexoyst_BoCcos(m8, n8);
        end
        if min_Pomex > Pomexoyst_BoCcos(m8, n8)
            min_Pomex = Pomexoyst_BoCcos(m8, n8);
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            m = 0; m8 = m*8;    
            if max_Pomex < Pomexoyst_BPSK(n8)
                max_Pomex = Pomexoyst_BPSK(n8);
            end
            if min_Pomex > Pomexoyst_BPSK(n8)
                min_Pomex = Pomexoyst_BPSK(n8);
            end
        end
for i = 1:size(BoCsin_Freq_L1_num, 1)
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    resrec = recalc_threshold(Pomexoyst_BoCsin(m8, n8), ntmin, ntmax, min_Pomex, max_Pomex);
    for freq_index = 1:fmax
        Qual_BoCsin_L1(m8, n8, freq_index, Pomex) = resrec;
    end
end
    for i = 1:size(BoCcos_Freq_L1_num, 1)
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        resrec = recalc_threshold(Pomexoyst_BoCcos(m8, n8), ntmin, ntmax, min_Pomex, max_Pomex);
        for freq_index = 1:fmax
            Qual_BoCcos_L1(m8, n8, freq_index, Pomex) = resrec;
        end
    end
        for i = 1:size(BPSK_Freq_L1_num, 1)
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            resrec = recalc_threshold(Pomexoyst_BPSK(n8), ntmin, ntmax, min_Pomex, max_Pomex);
            for freq_index = 1:fmax
                Qual_BPSK_L1(n8, freq_index, Pomex) = resrec;
            end
        end   
        
        

if Signal_OC 
% %%%%%%%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OC
C_Compl_OC = 0.4;
C_Search_OC = 0.15;
C_Accur_OC = 0.15;
C_RA_OC = 0.13;
C_Pomex_OC = 0.08;
C_MP_OC = 0.08;
C_InterNam_OC = 0.08;
C_InterIm_OC = 0.08;
C_IntraS_OC = 0.08;        
        
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Multi-objective optimization for L1OC signal\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!Antijamming<br>w=%.2f \n', C_Pomex_OC);
    fprintf('!Astronomy<br>w=%.2f  \n', C_RA_OC);
    fprintf('!Accurancy<br>w=%.2f  \n', C_Accur_OC);
    fprintf('!MultiPath<br>w=%.2f  \n', C_MP_OC);
    fprintf('!Inter from GPS<br>w=%.2f \n', C_InterNam_OC);
    fprintf('!Inter to GPS<br>w=%.2f \n', C_InterIm_OC);
    fprintf('!Search<br>w=%.2f\n', C_Search_OC);
    fprintf('!Complexity<br>w=%.2f\n', C_Compl_OC);
    fprintf('!Intra<br>w=%.2f\n', C_IntraS_OC);
    fprintf('!Norm');    
    
for i = 1:size(BoCsin_Freq_L1_num, 1)
    norm_s = 0;
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;
    if (m == 3.125)&&(n == 1.25)
%         n
    end
    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    fprintf('\n|- align="center"\n');
    fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Pomex));
    norm_s = norm_s + C_Pomex_OC*Qual_BoCsin_L1(m8, n8, freq_index, Pomex)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, RA));
    norm_s = norm_s + C_RA_OC*Qual_BoCsin_L1(m8, n8, freq_index, RA)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Accur));
    norm_s = norm_s + C_Accur_OC*Qual_BoCsin_L1(m8, n8, freq_index, Accur)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, MP));
    norm_s = norm_s + C_MP_OC*Qual_BoCsin_L1(m8, n8, freq_index, MP)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, InterNam));
    norm_s = norm_s + C_InterNam_OC*Qual_BoCsin_L1(m8, n8, freq_index, InterNam)^2;

    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, InterIm));
    norm_s = norm_s + C_InterIm_OC*Qual_BoCsin_L1(m8, n8, freq_index, InterIm)^2;

    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Search));
    norm_s = norm_s + C_Search_OC*Qual_BoCsin_L1(m8, n8, freq_index, Search)^2;

    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Compl));
    norm_s = norm_s + C_Compl_OC*Qual_BoCsin_L1(m8, n8, freq_index, Compl)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, IntraS));
    norm_s = norm_s + C_IntraS_OC*Qual_BoCsin_L1(m8, n8, freq_index, IntraS)^2;
    
    fprintf('|| %.4f ', sqrt(norm_s));

end

    for i = 1:size(BoCcos_Freq_L1_num, 1)
        norm_s = 0;
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;    
        fprintf('\n|- align="center"\n');
        fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Pomex));
        norm_s = norm_s + C_Pomex_OC*Qual_BoCcos_L1(m8, n8, freq_index, Pomex)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, RA));
        norm_s = norm_s + C_RA_OC*Qual_BoCcos_L1(m8, n8, freq_index, RA)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Accur));
        norm_s = norm_s + C_Accur_OC*Qual_BoCcos_L1(m8, n8, freq_index, Accur)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, MP));
        norm_s = norm_s + C_MP_OC*Qual_BoCcos_L1(m8, n8, freq_index, MP)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, InterNam));
        norm_s = norm_s + C_InterNam_OC*Qual_BoCcos_L1(m8, n8, freq_index, InterNam)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, InterIm));
        norm_s = norm_s + C_InterIm_OC*Qual_BoCcos_L1(m8, n8, freq_index, InterIm)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Search));
        norm_s = norm_s + C_Search_OC*Qual_BoCcos_L1(m8, n8, freq_index, Search)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Compl));
        norm_s = norm_s + C_Compl_OC*Qual_BoCcos_L1(m8, n8, freq_index, Compl)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, IntraS));
        norm_s = norm_s + C_IntraS_OC*Qual_BoCcos_L1(m8, n8, freq_index, IntraS)^2;

        fprintf('|| %.4f ', sqrt(norm_s));

    end
    
        for i = 1:size(BPSK_Freq_L1_num, 1)
            norm_s = 0;
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;    
            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Pomex));
            norm_s = norm_s + C_Pomex_OC*Qual_BPSK_L1(n8, freq_index, Pomex)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, RA));
            norm_s = norm_s + C_RA_OC*Qual_BPSK_L1(n8, freq_index, RA)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Accur));
            norm_s = norm_s + C_Accur_OC*Qual_BPSK_L1(n8, freq_index, Accur)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, MP));
            norm_s = norm_s + C_MP_OC*Qual_BPSK_L1(n8, freq_index, MP)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, InterNam));
            norm_s = norm_s + C_InterNam_OC*Qual_BPSK_L1(n8, freq_index, InterNam)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, InterIm));
            norm_s = norm_s + C_InterIm_OC*Qual_BPSK_L1(n8, freq_index, InterIm)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Search));
            norm_s = norm_s + C_Search_OC*Qual_BPSK_L1(n8, freq_index, Search)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Compl));
            norm_s = norm_s + C_Compl_OC*Qual_BPSK_L1(n8, freq_index, Compl)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, IntraS));
            norm_s = norm_s + C_IntraS_OC*Qual_BPSK_L1(n8, freq_index, IntraS)^2;

            fprintf('|| %.4f ', sqrt(norm_s));

        end 
fprintf('\n|} \n');
end

if ~Signal_OC 
% %%%%%%%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SC
C_Compl_SC = 1;
C_Search_SC = 4;
C_Accur_SC = 9;
C_RA_SC = 4;
C_Pomex_SC = 4;
C_MP_SC = 4;
C_InterNam_SC = 4;
C_InterIm_SC = 1;
C_IntraS_SC = 16^2;        
        
    fprintf('{| class="wikitable sortable" border="1" \n');
    fprintf('|+ Multi-objective optimization for L1SC signal\n');

    fprintf('|- align="center"\n');
    fprintf('!Signal \n');
    fprintf('!Antijamming<br>w=%.2f \n', C_Pomex_SC);
    fprintf('!Astronomy<br>w=%.2f  \n', C_RA_SC);
    fprintf('!Accurancy<br>w=%.2f  \n', C_Accur_SC);
    fprintf('!MultiPath<br>w=%.2f  \n', C_MP_SC);
    fprintf('!Inter from GPS<br>w=%.2f \n', C_InterNam_SC);
    fprintf('!Inter to GPS<br>w=%.2f \n', C_InterIm_SC);
    fprintf('!Search<br>w=%.2f\n', C_Search_SC);
    fprintf('!Complexity<br>w=%.2f\n', C_Compl_SC);
    fprintf('!Intra<br>w=%.2f\n', C_IntraS_SC);
    fprintf('!Norm');    
    
for i = 1:size(BoCsin_Freq_L1_num, 1)
    norm_s = 0;
    n = BoCsin_Freq_L1_num(i, 2); n8 = n*8;
    m = BoCsin_Freq_L1_num(i, 1); m8 = m*8;

    freq = BoCsin_Freq_L1_num(i, 3);
    freq_index = freq - 1558 + 1;    
    fprintf('\n|- align="center"\n');
    fprintf('|BoCsin(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Pomex));
    norm_s = norm_s + C_Pomex_SC*Qual_BoCsin_L1(m8, n8, freq_index, Pomex)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, RA));
    norm_s = norm_s + C_RA_SC*Qual_BoCsin_L1(m8, n8, freq_index, RA)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Accur));
    norm_s = norm_s + C_Accur_SC*Qual_BoCsin_L1(m8, n8, freq_index, Accur)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, MP));
    norm_s = norm_s + C_MP_SC*Qual_BoCsin_L1(m8, n8, freq_index, MP)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, InterNam));
    norm_s = norm_s + C_InterNam_SC*Qual_BoCsin_L1(m8, n8, freq_index, InterNam)^2;

    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, InterIm));
    norm_s = norm_s + C_InterIm_SC*Qual_BoCsin_L1(m8, n8, freq_index, InterIm)^2;

    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Search));
    norm_s = norm_s + C_Search_SC*Qual_BoCsin_L1(m8, n8, freq_index, Search)^2;

    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, Compl));
    norm_s = norm_s + C_Compl_SC*Qual_BoCsin_L1(m8, n8, freq_index, Compl)^2;
    
    fprintf('|| %.3f ', Qual_BoCsin_L1(m8, n8, freq_index, IntraS));
    norm_s = norm_s + C_IntraS_SC*Qual_BoCsin_L1(m8, n8, freq_index, IntraS)^2;
    
    fprintf('|| %.4f ', sqrt(norm_s));

end

    for i = 1:size(BoCcos_Freq_L1_num, 1)
        norm_s = 0;
        n = BoCcos_Freq_L1_num(i, 2); n8 = n*8;
        m = BoCcos_Freq_L1_num(i, 1); m8 = m*8;
        freq = BoCcos_Freq_L1_num(i, 3);
        freq_index = freq - 1558 + 1;    
        fprintf('\n|- align="center"\n');
        fprintf('|BoCcos(%.3f, %.3f) at %.0f<math>f_b</math> ', m, n, freq);

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Pomex));
        norm_s = norm_s + C_Pomex_SC*Qual_BoCcos_L1(m8, n8, freq_index, Pomex)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, RA));
        norm_s = norm_s + C_RA_SC*Qual_BoCcos_L1(m8, n8, freq_index, RA)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Accur));
        norm_s = norm_s + C_Accur_SC*Qual_BoCcos_L1(m8, n8, freq_index, Accur)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, MP));
        norm_s = norm_s + C_MP_SC*Qual_BoCcos_L1(m8, n8, freq_index, MP)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, InterNam));
        norm_s = norm_s + C_InterNam_SC*Qual_BoCcos_L1(m8, n8, freq_index, InterNam)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, InterIm));
        norm_s = norm_s + C_InterIm_SC*Qual_BoCcos_L1(m8, n8, freq_index, InterIm)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Search));
        norm_s = norm_s + C_Search_SC*Qual_BoCcos_L1(m8, n8, freq_index, Search)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, Compl));
        norm_s = norm_s + C_Compl_SC*Qual_BoCcos_L1(m8, n8, freq_index, Compl)^2;

        fprintf('|| %.3f ', Qual_BoCcos_L1(m8, n8, freq_index, IntraS));
        norm_s = norm_s + C_IntraS_SC*Qual_BoCcos_L1(m8, n8, freq_index, IntraS)^2;

        fprintf('|| %.4f ', sqrt(norm_s));

    end
    
        for i = 1:size(BPSK_Freq_L1_num, 1)
            norm_s = 0;
            n = BPSK_Freq_L1_num(i, 1); n8 = n*8;
            freq = BPSK_Freq_L1_num(i, 2);
            freq_index = freq - 1558 + 1;    
            fprintf('\n|- align="center"\n');
            fprintf('|BPSK(%.3f) at %.0f<math>f_b</math> ', n, freq);

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Pomex));
            norm_s = norm_s + C_Pomex_SC*Qual_BPSK_L1(n8, freq_index, Pomex)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, RA));
            norm_s = norm_s + C_RA_SC*Qual_BPSK_L1(n8, freq_index, RA)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Accur));
            norm_s = norm_s + C_Accur_SC*Qual_BPSK_L1(n8, freq_index, Accur)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, MP));
            norm_s = norm_s + C_MP_SC*Qual_BPSK_L1(n8, freq_index, MP)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, InterNam));
            norm_s = norm_s + C_InterNam_SC*Qual_BPSK_L1(n8, freq_index, InterNam)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, InterIm));
            norm_s = norm_s + C_InterIm_SC*Qual_BPSK_L1(n8, freq_index, InterIm)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Search));
            norm_s = norm_s + C_Search_SC*Qual_BPSK_L1(n8, freq_index, Search)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, Compl));
            norm_s = norm_s + C_Compl_SC*Qual_BPSK_L1(n8, freq_index, Compl)^2;

            fprintf('|| %.3f ', Qual_BPSK_L1(n8, freq_index, IntraS));
            norm_s = norm_s + C_IntraS_SC*Qual_BPSK_L1(n8, freq_index, IntraS)^2;

            fprintf('|| %.4f ', sqrt(norm_s));

        end 
fprintf('\n|} \n');
end

% %%%%%%%%%%%%%%%%%%%%%%%% Saving of results %%%%%%%%%%%%%%%%%%%%%%%
% save([path_to_prop '/Qual_BoCsin_L1.mat'], 'Qual_BoCsin_L1');
% save([path_to_prop '/Qual_BoCcos_L1.mat'], 'Qual_BoCcos_L1');
% save([path_to_prop '/Qual_BPSK_L1.mat'], 'Qual_BPSK_L1');
