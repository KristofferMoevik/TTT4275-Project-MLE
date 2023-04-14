clear all;
close all;

%% Define parameters



SNR_decibel = -10:10:60;
SNR = db2mag(SNR_decibel);

N = 513;
k = [10,12,14,16,18,20];

var_e_omega = NaN(size(SNR,2), size(k,2));
var_e_phi = NaN(size(SNR,2), size(k,2));
CRLB_omega = NaN(size(SNR,2), size(k,2));
CRLB_phi = NaN(size(SNR,2), size(k,2));

for i = 1:size(SNR,2)
    for j = 1:size(k,2)
        [var_e_omega(i,j), var_e_phi(i,j), CRLB_omega(i,j), CRLB_phi(i,j)] = MLE_simulation(N,k(j),SNR(i));
    end
end

%% Plot
close all
h = figure;
A = axes;
surfl(k,SNR,var_e_phi); hold on;
surfl(k,SNR,CRLB_phi); 
xlabel('k')
ylabel('SNR')
zlabel('Phase variance')
set(A,'ZScale','log')
colormap pink    % change color map
%shading interp    % interpolate colors across lines and faces