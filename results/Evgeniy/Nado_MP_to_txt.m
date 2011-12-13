close all
clc
clear

load('results/Evgeniy/Nado_MP.mat', 'Nado_MP');    

fid = fopen('Nado_MP.txt', 'w');

fprintf(fid, 'First digit: 1 - BOCsin; 2 - BOCcos; 3 - BPSK \nSecond number - m index (0 for BPSK) \nThird - n\n\n');
for Signal_Type = 1:3
    for m8 = 1:80
        for n8 = 1:80
            if Nado_MP(Signal_Type, m8, n8)
                m = m8/8*(Signal_Type~=3);
                n = n8/8;
                fprintf(fid, '%.0f %.3f %.3f\n', Signal_Type, m, n);
            end
        end
    end
end

fclose(fid);