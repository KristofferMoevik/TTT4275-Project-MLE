clear; close all;

%% Define

F_s = 10^6;
T = 1/F_s;
f_0 = 10^5;
omega_0 = 2*pi*f_0;
phi = pi/8;
A = 1;

SNR_decibel = 20;
SNR = db2mag(SNR_decibel);

sigma_square = A^2/(2*SNR);
sigma = sqrt(A^2/(2*SNR));

N = 513;
P = N*(N-1)/2;
Q = N*(N-1)*(2*N-1)/6;
n_0 = -P/N;

n = n_0:1:(n_0 + N-1);
t = 0:T:T*(N-1);

k = [10,12,14,16,18,20];
M = 2^k(3);

sim_steps = 1000;



%% Simulations

% Pre allocate vectors

omega_MLE = zeros(sim_steps,1);
phi_MLE = zeros(sim_steps,1);
error_omega = zeros(sim_steps,1);
error_phi = zeros(sim_steps,1);


for i = 1:sim_steps

    % Create signal x

    w_r = normrnd(0,sigma,1,N);
    w_i = 1i*normrnd(0,sigma,1,N);
    x = A*exp(1i*(omega_0*n*T + phi)) + w_r + w_i;

    % Fourier transform of x

    x = [x zeros(1,M-size(x,2))]; % zero padding
    X = fft(x,M);
    [val, m_star] = max(abs(X),[],2,'linear');

    % Find estimates and errors

    omega_MLE(i) = (2 * pi * m_star) / (M * T);
    phi_MLE(i) = angle(exp(-1i * omega_MLE(i) * n_0 *  T) * val);

    error_omega(i) = omega_0 - omega_MLE(i);
    error_phi(i) = phi - phi_MLE(i);

    % Simulation progress

    if mod(i, sim_steps/100) == 0 && true
        fprintf('%i%%\n', 100 * i/sim_steps)
    end

end

var_error_omega = var(error_omega);
var_error_phi = var(error_phi);


%% CRLB

CRLB_omega = (12 * sigma_square) / (A^2 * T^2 * N * (N^2 - 1));
CRLB_phi = (12 * sigma_square * (n_0^2 * N + 2 * n_0 * P + Q)) / (A^2 * N^2 * (N^2 - 1));



%% Plot

fig1 = figure;
plot(abs(error_omega));
grid on;
xlabel('Steps');





