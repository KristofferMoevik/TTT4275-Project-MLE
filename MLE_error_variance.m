function [omega_error_var, phi_error_var] = MLE_error_variance(N, steps, SNR_dB, k)

    %MLE_SIMULATION2 Simulate the performance of the FFT-based MLE
    
    % Define parameters
    F_s = 10^6;
    T = 1/F_s;
    f_0 = 10^5;
    omega_0 = 2*pi*f_0;
    phi = pi/8;
    P = N*(N-1)/2;
    n_0 = -P/N;
    M = 2^k;
    
    % Pre allocate vectors
    
    omega_MLE = zeros(steps,1);
    phi_MLE = zeros(steps,1);
    omega_error = zeros(steps,1);
    phi_error = zeros(steps,1);
    
    for i = 1:steps
    
        % Create signal x
    
        x = signal(N,SNR_dB);
    
        % Fourier transform of x
    
        x = [x zeros(1,M-size(x,2))]; % zero padding of x
        x_fft = fft(x,M);
        [val, m_star] = max(abs(x_fft),[],2,'linear');
    
        % Find estimates and errors
    
        omega_MLE(i) = (2 * pi * m_star) / (M * T);
        phi_MLE(i) = angle(exp(-1i * omega_MLE(i) * n_0 *  T) * val);
    
        omega_error(i) = omega_0 - omega_MLE(i);
        phi_error(i) = phi - phi_MLE(i);
    
    end
    
    omega_error_var = var(omega_error);
    phi_error_var = var(phi_error);

    if omega_error_var < 10^-20
        omega_error_var = 0;
    end
    if phi_error_var < 10^-20
        phi_error_var = 0;
    end
end




