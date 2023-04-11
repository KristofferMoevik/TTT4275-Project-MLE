clear all;
close all;

%% Define parameters

F_s = 10^6;
T = 1/F_s;

f_0 = 10^5;
omega_0 = 2*pi*f_0;
phi = pi/8;
A = 1;

SNR_decibel = 1;
SNR = db2mag(SNR_decibel);
sigma_square = A^2/2*SNR;
sigma = sqrt(A^2./(2*SNR));


N = 513;
P = N*(N-1)/2;
Q = N*(N-1)*(2*N-1)/6;
n_0 = -P/N;
t = 0:T:T*(N-1);
k = [10,12,14,16,18,20];
M = 2^10;

%% Define signal

x = A* exp(1i*(omega_0*t + phi)) + normrnd(0,sigma,1,N) + 1i*normrnd(0,sigma,1,N);

%% Find CRLB
CRLB_freq = 12*sigma_square / A^2*T^2*N*(N^2-1);
CRLB_phase = 12*sigma_square*(n_0^2*N+2*n_0*P+Q) / (A^2*N^2*(N^2-1));

%% Fourier transform of x
% zero pad x
x = [x, zeros(1,M - size(x,2))];
Y = fft(x,M);
[val, m_star] = max(Y);

%% Find omega_hat
omega_hat = (2*pi*m_star) / (M*T);

%% Find phi_hat
phi_hat = angle(exp(-1i*omega_hat*n_0*T) * val);

%% Plot

plot(abs(Y));
