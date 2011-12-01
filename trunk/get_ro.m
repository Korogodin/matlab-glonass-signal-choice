function ro_out = get_ro( m, n, Signal_Type, path_to_ro )
%/**
% Считываение из mat-файла АКФ
%
%@param Signal_Type задает BoCsin, BoCcos, BPSK - 1, 2, 3 соответственно
%@param m, n - Индексы m, n (BPSK <-> m = 0)
%@param path_to_ro - путь до mat-файлов
%*/

    try
        if Signal_Type == 1
            load([path_to_ro '/ro_BoCsin(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').mat'])
        elseif Signal_Type == 2
            load([path_to_ro '/ro_BoCcos(' sprintf('%.3f', m) ', ' sprintf('%.3f', n) ').mat'])
        elseif Signal_Type == 3
            load([path_to_ro '/ro_BoCsin(' sprintf('%.3f', 0) ', ' sprintf('%.3f', n) ').mat'])
        end
        ro_out = ro;
    catch exception
        ro_out = 0; % Если файла нет
    end
    
end

