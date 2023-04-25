close all;
clear;

% Evaluate performance of MLE estimators and tune them

% Define parameters and variables
N = 513;        % signal length
steps = 500;    % simulation length
SNR_dB = [-10 0 10 20 30 40 50 60];
k = [10 12 14 16 18 20];
n_SNR = length(SNR_dB);
n_k = length(k);

% Allocate
omega_MLE_error_var = zeros(n_SNR,n_k);
phi_MLE_error_var = zeros(n_SNR,n_k);
omega_MLE_error_mean = zeros(n_SNR,n_k);
phi_MLE_error_mean = zeros(n_SNR,n_k);
omega_CRLB = zeros(n_SNR,1);
phi_CRLB = zeros(n_SNR,1);
omega_opt_error_var = zeros(n_SNR,1);
phi_opt_error_var = zeros(n_SNR,1);
omega_opt_error_mean = zeros(n_SNR,1);
phi_opt_error_mean = zeros(n_SNR,1);

% Simulation

for i = 1:n_SNR
    for j = 1:n_k

        % Find estimation error
        [omega_MLE_error_var(i,j), phi_MLE_error_var(i,j), omega_MLE_error_mean(i,j), phi_MLE_error_mean(i,j)] = MLE_error_stats(N,steps,SNR_dB(i),k(j));

    end
    
    % Find CRLB
    [omega_CRLB(i), phi_CRLB(i)] = calculate_CRLB(N,SNR_dB(i));
    
    % Find optimized estimate error variance
    [omega_opt_error_var(i), phi_opt_error_var(i), omega_opt_error_mean(i), phi_opt_error_mean(i)] = opt_error_stats(N,steps,SNR_dB(i));

    % Simulation progress
    fprintf('%d%%\n', round(100 * i/n_SNR));

end

% Save results
save('sim_data_2.mat', 'omega_MLE_error_var', 'phi_MLE_error_var', 'omega_CRLB', 'phi_CRLB', 'omega_opt_error_var', 'phi_opt_error_var', 'omega_opt_error_mean', 'phi_opt_error_mean', 'omega_MLE_error_mean', "phi_MLE_error_mean");


