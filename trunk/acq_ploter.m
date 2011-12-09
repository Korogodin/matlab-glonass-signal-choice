clear 
close all
clc

path_to_results = [pwd 'results/Evgeniy'];
load('results/Evgeniy/AcqBPSK.mat', 'AcqBPSK');
load('results/Evgeniy/AcqBOC.mat', 'AcqBOC');

n8max = 80;
n_arr = (1:n8max)/8;
hF = 0;
hF = figure(hF + 1);
plot(n_arr, AcqBOC)
xlabel('n')
ylabel('t, s');
grid on
title('Acquation time (1 correlator) for BOC(m, n) signals')


hF = 0;
hF = figure(hF + 1);
plot(n_arr, AcqBPSK, n_arr, AcqBOC)
xlabel('n')
ylabel('t, s');
grid on
title('Acquation time (1 correlator) for BOC(m, n) and BPSK(n) signals')

fprintf('<!DOCTYPE html> \n');
fprintf('<html><head><title>!DOCTYPE</title> \n');
fprintf('<meta charset="utf-8"> </head> <body> \n');   

fprintf('<table>\n');
fprintf('<tr><td><b>n</b></td><td>Time for BOC(m, n), hours</td><td>Time for BPSK(n), hours</td></tr>\n');
for n8 = 1:n8max
    
fprintf('<tr><td><b>%.3f</b></td><td>%.1f</td><td>%.2f</td></tr>\n', n8/8, AcqBOC(n8) / 3600, AcqBPSK(n8) / 3600);    
   
    
end
fprintf('</table>\n');
fprintf('</body> </html> \n');
