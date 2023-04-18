function [var_e_omega, var_e_phi, CRLB_omega, CRLB_phi] = MLE_simulation(N,k,SNR, sim_steps)


    F_s = 10^6;
    T = 1/F_s;
    
    f_0 = 10^5;
    omega_0 = 2*pi*f_0;
    phi = pi/8;
    A = 1;  

    P = N*(N-1)/2;
    Q = N*(N-1)*(2*N-1)/6;
    n_0 = -P/N;
    n = n_0:1:(n_0 + N-1);
    t = 0:T:T*(N-1);
    
    sigma_square = A^2/2*SNR;
    sigma = sqrt(A^2./(2*SNR));
    
    %% Find CRLB
    CRLB_omega = (12 * sigma_square) / (A^2 * T^2 * N * (N^2 - 1));
    CRLB_phi = (12 * sigma_square * (n_0^2 * N + 2 * n_0 * P + Q)) / (A^2 * N^2 * (N^2 - 1));

    
    omega_hat = [];
    phi_hat = [];
    omega_hat = [];
    e_omega = [];
    e_phi = [];
    for i=1:sim_steps
        %% Define signal
        w_r = normrnd(0,sigma,1,N);
        w_i = 1i*normrnd(0,sigma,1,N);
    
        x = A* exp(1i*(omega_0*n*T + phi)) + w_r + w_i;
    
        %% Fourier transform of x
    
        M = 2^k;
        x = [x zeros(1,M-size(x,2))];
    
        Y = fft(x,M);
        %Y = Y .* [ones(1,0.5*M), zeros(1,0.5*M)];
        [val, m_star] = max(Y);
    
        %% Find estimates
    
        omega_hat = [omega_hat; (2*pi*m_star) / (M*T)];
        phi_hat = [phi_hat; angle(exp(-1i*omega_hat(i)*n_0*T) * val)];
        omega_hat = [omega_hat; (2*pi*m_star) / (T*M)];

        %% Estimation errors
    
        e_omega = [e_omega; omega_0 - omega_hat(i)];
        e_phi = [e_phi; phi - phi_hat(i)];
    end
    var_e_omega = var(e_omega);
    var_e_phi = var(e_phi);

end

