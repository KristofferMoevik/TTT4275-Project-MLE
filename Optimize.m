clear;

%% Get estimators as start value

% Parameters

N = 513; 
k = 10;
SNR_dB = [-10 0 10 20 30 40 50 60];
F_s = 10^6;
T = 1/F_s;
A = 1;
phi = pi/8;
P = N*(N-1)/2;
n_0 = -P/N;
n = n_0:1:(n_0 + N-1);
M = 2^k;

% Signal x
    
x = signal(N,SNR_dB(1));

% FFT of x

x_zeropad = [x zeros(1,M-size(x,2))]; % zero padding of x
x_fft = fft(x_zeropad,M);
[val, m_star] = max(x_fft,[],2,'linear');
    
% Find MLE estimators
    
omega_MLE = (2 * pi * m_star) / (M * T);
phi_MLE = angle(exp(-1i * omega_MLE * n_0 *  T) * val);



%% Optimzation

MLE = [omega_MLE; phi_MLE];

% Estimated signal

x_hat = A*exp(1i*(MLE(1)*n*T + MLE(2)));

% Object function

obj_func = norm(abs(x-x_hat));

% Find minimum

options = optimset('Display','iter');

MLE_opt = fminsearch(@(MLE)obj_func, MLE, options);

fval = norm(abs(x - A*exp(1i*(MLE_opt(1)*n*T + MLE_opt(2)))));

%% Print

disp('Estimates before optimization:');
disp(['omega_MLE = ', num2str(MLE(1))]);
disp(['phi_MLE = ', num2str(MLE(2))]);
disp(['Obj func = ', num2str(obj_func)]);
disp(' ');
disp('Estimates after optimization:');
disp(['omega_opt = ', num2str(MLE_opt(1))]);
disp(['phi_opt = ', num2str(MLE_opt(2))]);
disp(['Obj func = ', num2str(fval)]);








