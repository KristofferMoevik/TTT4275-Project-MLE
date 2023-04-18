clear;

%% Define

N = 513; 
k = 10;

% Create signal x
    
w_r = normrnd(0,sigma,1,N);
w_i = 1i*normrnd(0,sigma,1,N);
x = A*exp(1i*(omega_0*n*T + phi)) + w_r + w_i;
