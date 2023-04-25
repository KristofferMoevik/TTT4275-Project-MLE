function [omega_opt, phi_opt] = optimize_MLE1(omega_MLE, phi_MLE, SNR_dB)

    N = 513;
    A = 1;
    P = N*(N-1)/2;
    n_0 = -P/N;
    n = n_0:1:(n_0 + N-1);
    F_s = 10^6;
    T = 1/F_s;
    
    MLE = [omega_MLE; phi_MLE];

    % Signal x
    x = signal(N, SNR_dB);

    % Object function
    obj_func = @(MLE) norm(x - A*exp(1i*(MLE(1)*n*T + MLE(2))));
       
    % Find minimum
    opt = fminunc(obj_func, MLE);

    omega_opt = opt(1);
    phi_opt = opt(2);
end
