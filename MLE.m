clear all;
close all;

%% Define parameters

F_s = 10^6;
T = 1/F_s;

f_0 = 10^5;
omega_0 = 2*pi*f_0;
phi = pi/8;
A = 1;

SNR_decibel = 0;
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

x = A* exp(1i*(omega_0*t + phi)) + normrnd(0,sigma,1,N) + 1i*normrnd(0,sigma,1,N);

%% Find CRLB
CRLB_freq = 12*sigma_square / A^2*T^2*N*(N^2-1);
CRLB_phase = 12*sigma_square*(n_0^2*N+2*n_0*P+Q) / (A^2*N^2*(N^2-1));


%% Plot

plot(abs(x));
