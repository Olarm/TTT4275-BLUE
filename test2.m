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
snr_db = 30;
s = sigma(snr_db);
s = A^2 / 2 ./ db2mag(snr_db);
mean = 1;
var = (A^2 / 2) ./ db2mag(snr_db);

figure(1)
subplot(3,1,1)
x_martin = gen_signal_martin(w_0, n, A, T, phi, 0, sqrt(var));
plot(n, x_martin);
title('Input signal affected by white noise')
grid on

subplot(3,1,2)
x_martin_angle = angle(x_martin);
plot(n, x_martin_angle);
title('Angle of noisy signal')
grid on

subplot(3,1,3)
x_martin_angle_unwrapped = unwrap(x_martin_angle);
plot(n, x_martin_angle_unwrapped);
title('Unwrapped angle of noisy signal')
grid on

figure(2)
subplot(3,1,1)
x_ola = gen_signal_ola(w_0, n, A, T, phi, 0, sqrt(var), mean, N, s);
plot(n, x_ola);
title('Input signal affected by white noise')
grid on

subplot(3,1,2)
x_ola_angle = angle(x_ola);
plot(n, x_ola_angle);
title('Angle of noisy signal')
grid on

subplot(3,1,3)
x_ola_angle_unwrapped = unwrap(x_ola_angle);
plot(n, x_ola_angle_unwrapped);
title('Unwrapped angle of noisy signal')
grid on




%%Function definitions
function x = gen_signal_martin(omega_0, n, A, T, phi, w_e, w_sigma)

w = gen_noise_martin(length(n), w_e, w_sigma);
x = A * exp(1j * (omega_0 * n * T + phi)) + w;

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