clear all;
close all;

%% Define parameters

F_s = 10^6;
T = 10e-6;

f_0 = 10^5;
omega_0 = 2*pi*f_0;
phi = pi/8;
A = 1;

SNR_decibel = 0;
SNR = db2mag(SNR_decibel);
sigma = sqrt(A^2./(2*SNR));

N = 513;
n_0 = -256;
t = 0:T:T*(N-1);
k = [10,12,14,16,18,20];

t1 = A* exp(1i*(omega_0*t + phi));
t2 = sigma*randn(N,1);
t3 = 1i*sigma*randn(N,1);

x = A* exp(1i*(omega_0*t + phi)); %%+ sigma*randn(1,N) + 1i*sigma*randn(1,N);

plot(real(x)); hold on;
%plot(imag(x)); hold off;