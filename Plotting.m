clear;
close all;

%% Define parameters

N = 513;        % signal length
steps = 100;    % simulation length
SNR_dB = [-10 0 10 20 30 40 50 60];
k = [10 12 14 16 18 20];
n_SNR = length(SNR_dB);
n_k = length(k);

%% Plotting

% Frequency estimatation error variance 

load('sim_data_1.mat');

fig1 = figure;
semilogy(SNR_dB,omega_CRLB); hold on; grid on;
for i = 1:n_k
    plot(SNR_dB, omega_MLE_error_var(:,i)); hold on;
end
legend('CRLB','k = 10','k = 12','k = 14','k = 16','k = 18','k = 20');
title('Variance of frequency estimation error');
xlabel('SNR (dB)');
ylabel('Variance');

dims = [200 200 500 300]; % [x_pos, y_pos, x_brd, y_brd]
set(fig1, 'renderer', 'painters', 'position', dims, 'PaperPositionMode', 'auto');
%print('-dpng', '-r600', '/Users/johannesskaro/Documents/KYB 3.år/EDK/Project/omega_err_var_plot2');



% Phase estimation error variance

fig2 = figure;
semilogy(SNR_dB,phi_CRLB); hold on;  grid on;
for i = 1:n_k
    plot(SNR_dB, phi_MLE_error_var(:,i)); hold on;
end
legend('CRLB','k = 10','k = 12','k = 14','k = 16','k = 18','k = 20');
title('Variance of phase estimation error');
xlabel('SNR (dB)');
ylabel('Variance');

dims = [200 200 500 300]; % [x_pos, y_pos, x_brd, y_brd]
set(fig2, 'renderer', 'painters', 'position', dims, 'PaperPositionMode', 'auto');
%print('-dpng', '-r600', '/Users/johannesskaro/Documents/KYB 3.år/EDK/Project/phi_err_var_plot2');


% Frequency estimation opt error variance

fig3 = figure;
semilogy(SNR_dB,omega_CRLB); hold on;
semilogy(SNR_dB, omega_MLE_error_var(:,1)); hold on;
plot(SNR_dB, omega_MLE_error_var(:,6)); hold on;
plot(SNR_dB, omega_opt_error_var);
grid on;
xlabel('SNR (dB)');
ylabel('Variance');
title('Variance of frequency estimate errors with tuned estimate');
legend('CRLB','k = 10', 'k = 20', 'tuned');

dims = [200 200 500 300]; % [x_pos, y_pos, x_brd, y_brd]
set(fig3, 'renderer', 'painters', 'position', dims, 'PaperPositionMode', 'auto');
%print('-dpng', '-r600', '/Users/johannesskaro/Documents/KYB 3.år/EDK/Project/omega_opt_err_plot');


% Phase estimation opt error variance

fig4 = figure;
semilogy(SNR_dB, phi_MLE_error_var(:,1)); hold on;
plot(SNR_dB, phi_MLE_error_var(:,6)); hold on;
plot(SNR_dB, phi_opt_error_var);
grid on;
xlabel('SNR (dB)');
ylabel('Variance');
title('Variance of phase estimate errors with tuned estimate');
legend('k = 10', 'k = 20', 'tuned');

dims = [200 200 500 300]; % [x_pos, y_pos, x_brd, y_brd]
set(fig4, 'renderer', 'painters', 'position', dims, 'PaperPositionMode', 'auto');
%print('-dpng', '-r600', '/Users/johannesskaro/Documents/KYB 3.år/EDK/Project/phi_opt_err_plot');


% Frequency estimation opt error mean

fig5 = figure;
plot(SNR_dB, omega_MLE_error_mean(:,1)); hold on;
plot(SNR_dB, omega_MLE_error_mean(:,6)); hold on;
plot(SNR_dB, omega_opt_error_mean);
grid on;
xlabel('SNR (dB)');
ylabel('Hz');
title('Mean of frequency estimate errors with tuned estimate');
legend('k = 10', 'k = 20', 'tuned');

dims = [200 200 500 300]; % [x_pos, y_pos, x_brd, y_brd]
set(fig5, 'renderer', 'painters', 'position', dims, 'PaperPositionMode', 'auto');
%print('-dpng', '-r600', '/Users/johannesskaro/Documents/KYB 3.år/EDK/Project/omega_opt_err_mean_plot');


% Phase estimation opt error mean

fig6 = figure;
plot(SNR_dB, phi_MLE_error_mean(:,1)); hold on;
plot(SNR_dB, phi_MLE_error_mean(:,6)); hold on;
plot(SNR_dB, phi_opt_error_mean);
grid on;
xlabel('SNR (dB)');
ylabel('Rad');
title('Mean of phase estimate errors with tuned estimate');
legend('k = 10', 'k = 20', 'tuned');

dims = [200 200 500 300]; % [x_pos, y_pos, x_brd, y_brd]
set(fig6, 'renderer', 'painters', 'position', dims, 'PaperPositionMode', 'auto');
%print('-dpng', '-r600', '/Users/johannesskaro/Documents/KYB 3.år/EDK/Project/phi_opt_err_mean_plot');




