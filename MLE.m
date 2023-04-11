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

%% Define signal

x = A* exp(1i*(omega_0*t + phi)); %%+ sigma*randn(1,N) + 1i*sigma*randn(1,N);

plot(real(x)); hold on;
%plot(imag(x)); hold off;

w_r = normrnd(0,sigma,1,N);
w_i = 1i*normrnd(0,sigma,1,N);

x = A* exp(1i*(omega_0*t + phi)) + w_r + w_i;

%% Find CRLB
CRLB_freq = 12*sigma_square / A^2*T^2*N*(N^2-1);
CRLB_phase = 12*sigma_square*(n_0^2*N+2*n_0*P+Q) / (A^2*N^2*(N^2-1));

%% Fourier transform of x

M = 2^k(6);
x = [x zeros(1,M-size(x,2))];

Y = fft(x,M);
Y = Y .* [ones(1,0.5*M), zeros(1,0.5*M)];
[val, m_star] = max(Y);

%% Find estimates

omega_hat = (2*pi*m_star) / (M*T);
phi_hat = angle(exp(-1i*omega_hat*n_0*T) * val);


%% Estimation errors

e_omega = omega_0 - omega_hat;
e_phi = phi - phi_hat;

%% Plot

plot(abs(Y));

