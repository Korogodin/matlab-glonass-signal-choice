clc
close all
clear

MP_BOCsin_results = [...
0.25 0.25   1.023
0.5 0.5   2.046
1 1   4.092
3.25 1.75  10.230
6.125 1.75  16.112
3.25 1.875  10.486
3.75 1.875  11.509
5.125 1.875  14.322
5.625 1.875  15.345
3.25 2  10.741
4 2  12.276
5 2  14.322
6 2  16.368
4.25 2.125  13.043
4.5 2.25  13.810
5.625 2.25  16.112
4.75 2.375  14.578
3.125 2.5  11.509
5 2.5  15.345
5.25 2.625  16.112];

MP_BOCcos_results = [...
0.25 0.25   1.023
3.25 1.75  10.230
3.25 1.875  10.486
5.125 1.875  14.322
3 2  10.230
3.25 2  10.741
3.375 2.25  11.509
3.125 2.5  11.509
2.75 2.75  11.253
2.875 2.875  11.764
3 3  12.276];

MP_BPSK_results = [...
0.5   1.023
0.75   1.534
1   2.046
1.25   2.557
1.5   3.069
2   4.092
2.5   5.115
3   6.138];


% %%%%%%%%%%%%%%%%%%%% MultiPath %%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_Evgeniy = pwd;
load([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
load([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
load([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');


for i = 1:size(MP_BOCsin_results, 1)
    m = MP_BOCsin_results(i, 1); m8 = m*8;
    n = MP_BOCsin_results(i, 2); n8 = n*8;
    if ~isnan(MP_BOCsin(m8, n8))
        fprintf('BOCsin(%.3f, %.3f) MP %f changed to %f \n', m, n, MP_BOCsin(m8, n8), MP_BOCsin_results(i, 3));
    end
    MP_BOCsin(m8, n8) = MP_BOCsin_results(i, 3);    
end
for i = 1:size(MP_BOCcos_results, 1)
    m = MP_BOCcos_results(i, 1); m8 = m*8;
    n = MP_BOCcos_results(i, 2); n8 = n*8;
    if ~isnan(MP_BOCcos(m8, n8))
        fprintf('BOCcos(%.3f, %.3f) MP %f changed to %f \n', m, n, MP_BOCcos(m8, n8), MP_BOCcos_results(i, 3));
    end
    MP_BOCcos(m8, n8) = MP_BOCcos_results(i, 3);    
end
for i = 1:size(MP_BPSK_results, 1)
    n = MP_BPSK_results(i, 1); n8 = n*8;
    if ~isnan(MP_BPSK(n8))
        fprintf('BPSK(%.3f) MP %f changed to %f \n', n, MP_BPSK(n8), MP_BPSK_results(i, 2));
    end
    MP_BPSK(n8) = MP_BPSK_results(i, 2);    
end

% save([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
% save([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
% save([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');


