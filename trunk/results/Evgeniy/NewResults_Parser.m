clc
clear
close all

path_to_Evgeniy = pwd;
load([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
load([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
load([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');

NewRes = [ ...
1 1 0.5  10.495
1 3.75 3.5   2.973
1 4 3.25   2.785
1 4 3.5   2.759
1 4.5 2.5   2.550
1 4.5 3.5   2.412
1 4.75 2.5   2.402
1 5 3   2.203
2 0.5 0.5  10.706
2 1 0.5   6.161
2 1 1   6.821
2 1.5 1   4.736
2 3.5 3.5   1.865
2 3.75 3.5   1.451
2 3.75 3.75   1.698
2 4 2.5   1.515
2 4 3.25   1.400
2 4 3.5   1.368
2 4 4   1.557
2 4.25 2.5   1.425
2 4.5 2.5   1.336
2 4.5 3.5   1.203
2 4.75 2.5   1.273
2 5 2.5   1.323
2 5 3   1.151
3 1.75 0  11.801
3 2.25 0  10.025
3 9.5 0   2.014
3 10 0   1.874
1 0.250 0.250 16.778
1 0.500 0.500 15.022
1 1.000 1.000 10.765
1 4.500 2.250 2.594
1 4.750 2.375 2.435
1 5.000 2.000 2.338
1 5.000 2.500 2.310
1 5.250 2.625 2.153
1 5.625 2.250 2.027
1 6.000 2.000 1.938
2 0.250 0.250 13.772
2 2.750 2.750 2.524
2 2.875 2.875 2.392
2 3.000 3.000 2.261
3 0.500 0   18.222
3 0.750 0   16.814
3 1.000 0   15.377
3 1.250 0   14.025
3 1.500 0   12.830
3 2.000 0   10.872
3 2.500 0   9.286
3 3.000 0   8.008];


BOCsin = 1; BOCcos = 2; BPSK = 3;
for j = 1:size(NewRes, 1)
    switch NewRes(j, 1)
        case BOCsin
            MP_BOCsin(NewRes(j, 2)*8, NewRes(j, 3)*8) = NewRes(j, 4);
        case BOCcos
            MP_BOCcos(NewRes(j, 2)*8, NewRes(j, 3)*8) = NewRes(j, 4);
        case BPSK
            MP_BPSK(NewRes(j, 2)*8) = NewRes(j, 4);
    end
end

save([path_to_Evgeniy '/MP_BOCsin.mat'], 'MP_BOCsin');
save([path_to_Evgeniy '/MP_BOCcos.mat'], 'MP_BOCcos');
save([path_to_Evgeniy '/MP_BPSK.mat'], 'MP_BPSK');



