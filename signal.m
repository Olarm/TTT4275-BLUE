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

snr_db = -10;
s = sigma(snr_db);
mean = 0;

v = s * randn(N, 1).' + mean;
snr = 1 / (2 * s^2);

log_x = 1i * w_0 * n * T + phi + v;
x_p = w_0 * n * T + phi + v;
x = A * exp(log_x);
x_u = unwrap(x_p);

figure
subplot(2,1,1)
plot(n, x);
title(['SNR: ',  num2str(10*log10(snr)), ' db']);
legend('x')
subplot(2,1,2)
plot(n, x_u)
legend('unwrap(x)')
