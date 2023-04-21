clear;

%% Get estimators as start value

% Parameters

N = 513; 
k = 10;
SNR_dB = [-10 0 10 20 30 40 50 60];
n_SNR = length(SNR_dB);
F_s = 10^6;
f_0 = 10^5;
omega_0 = 2*pi*f_0;
T = 1/F_s;
A = 1;
phi = pi/8;
P = N*(N-1)/2;
n_0 = -P/N;
n = n_0:1:(n_0 + N-1);
M = 2^k;

% Allocate vectors

omega_0_vec = zeros(n_SNR,1);
phi_vec = zeros(n_SNR,1);
omega_MLE = zeros(n_SNR,1);
phi_MLE = zeros(n_SNR,1);
omega_opt = zeros(n_SNR,1);
phi_opt = zeros(n_SNR,1);


%% Optimzation

for k = 1:n_SNR
    
    % Signal x
    
    x = signal(N,SNR_dB(k));

    % FFT of x

    x_zeropad = [x zeros(1,M-size(x,2))]; % zero padding of x
    x_fft = fft(x_zeropad,M);
    [val, m_star] = max(x_fft,[],2,'linear');
    
    % Find MLE estimators
    
    omega_MLE(k) = (2 * pi * m_star) / (M * T);
    phi_MLE(k) = angle(exp(-1i * omega_MLE(k) * n_0 *  T) * val);
    MLE = [omega_MLE(k); phi_MLE(k)];
    
    % Estimated signal
    
    x_hat = A*exp(1i*(MLE(1)*n*T + MLE(2)));
    
    % Object function
    
    obj_func = @(MLE) norm(x - A*exp(1i*(MLE(1)*n*T + MLE(2))));
    obj_func_val = norm(x - A*exp(1i*(MLE(1)*n*T + MLE(2))));
    
    % Find minimum
    
    options = optimset('Display','iter');
    [MLE_opt, fval] = fminsearch(obj_func, MLE, options);
        
    omega_opt(k) = MLE_opt(1);
    phi_opt(k) = MLE_opt(2);

    omega_0_vec(k) = omega_0;
    phi_vec(k) = phi;

end 


%% Print

disp('Estimates before optimization:');
disp(['omega_MLE = ', num2str(MLE(1))]);
disp(['phi_MLE = ', num2str(MLE(2))]);
disp(['Obj func = ', num2str(obj_func_val)]);
disp(' ');
disp('Estimates after optimization:');
disp(['omega_opt = ', num2str(MLE_opt(1))]);
disp(['phi_opt = ', num2str(MLE_opt(2))]);
disp(['Obj func = ', num2str(fval)]);

%% Plotting

fig1 = figure;
plot(SNR_dB,omega_0_vec); hold on;
plot(SNR_dB,omega_MLE); hold on;
plot(SNR_dB,omega_opt);
grid on;
xlabel('SNR (dB)');
ylabel('Hz');
title('Tuning of frequency estimator');
legend('omega_0', 'omega_{MLE}', 'omega_{opt}');


fig2 = figure;
plot(SNR_dB,phi_vec); hold on;
plot(SNR_dB,phi_MLE); hold on;
plot(SNR_dB,phi_opt);
grid on;
xlabel('SNR (dB)');
ylabel('rad')
title('Tuning of phase estimator');
legend('phi', 'phi_{MLE}', 'phi_{opt}');






