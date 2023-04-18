function x = signal(N,SNR_dB)

    % Create signal x
    F_s = 10^6;
    T = 1/F_s;
    f_0 = 10^5;
    omega_0 = 2*pi*f_0;
    phi = pi/8;
    A = 1;
    P = N*(N-1)/2;
    n_0 = -P/N;
    n = n_0:1:(n_0 + N-1);
    SNR = db2mag(SNR_dB);
    sigma = sqrt(A^2/(2*SNR));
    
    w_r = normrnd(0,sigma,1,N);
    w_i = 1i*normrnd(0,sigma,1,N);
    x = A*exp(1i*(omega_0*n*T + phi)) + w_r + w_i;

end

