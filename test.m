clear all;
%%Define parameters
F_s = 10e6;
T = 10e-6;
f_0 = 10e5;
w_0 = 2 * pi * f_0;
phi = pi / 8;
A = 1;
N = 513;
n = -256:1:N-257;
%n = 1:513;
size(n)
snr_db = 10;
s1 = sigma(snr_db);
s2 = A^2 / 2 ./ db2mag(snr_db);
mean = 0;
var = (A^2 / 2) ./ db2mag(snr_db);

%%Generate signals
figure(1)
x_martin = gen_signal_martin(w_0, n, A, T, phi, 0, sqrt(var));
plot(n, x_martin);
title('Martins signal')

figure(2)
x_ola = gen_signal_ola(w_0, n, A, T, phi, 0, sqrt(var), mean, N, s1);
plot(n, x_ola);
title('Our signal')


%%Function definitions
function x = gen_signal_martin(omega_0, n, A, T, phi, w_e, w_sigma)
v = normrnd(w_e, w_sigma, 1, length(n));
w = gen_noise_martin(length(n), w_e, w_sigma);
x = A * exp(1j * (omega_0 * n * T + phi + v));

end

function v = gen_noise_martin(N, w_e, w_sigma)
    v = normrnd(w_e, w_sigma, 1, N) + 1j * normrnd(w_e, w_sigma, 1, N);
end

function x = gen_signal_ola(w_0, n, A, T, phi, w_e, w_sigma, mean, N, s)
    v = s * randn(N, 1).' + mean;
    snr = 1 / (2 * s^2);
    snr_log = 10 * log10(snr);
    %snr_mag = db2mag(snr_log);
    log_x = 1i * ( w_0 * n * T + phi + v);
    x_p = w_0 * n * T + phi + v;
    x = A * exp( log_x);
end

function s = sigma(snr_db)
    [~, w] = size(snr_db);
    s = zeros(1,w);
    for i=1:w
        s(i) = 10^(-1/2*log10(2) - 1/20*snr_db(i));
    end
end