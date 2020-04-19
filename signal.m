F_s = 10e6;
T = 10e-6;
f_0 = 10e5;
w_0 = 2 * pi * f_0;
phi = pi / 8;
A = 1;
N = 513;
n_0 = -256;
n = -256:1:N-257;

P = N * ( N - 1) / 2;
Q = N * ( N - 1) * (2*N - 1) / 6;

snr_db = 40;
s = sigma(snr_db);
%s = A^2 / 2 ./ db2mag(snr_db);
mean = 0;

v = s * randn(N, 1).' + mean;
snr = 1 / (2 * s^2);
snr_log = 10 * log10(snr);
snr_mag = db2mag(snr_log);

var_w_hat = 12*s^2 / (A^2 * T^2 * N ( N^2 - 1));
var_phi_hat = 12*s^2 * (n_0^2 * N + 2*n_0 * P + Q) / (A^2 * N^2 * ( N^2 - 1));

log_x = 1i * ( w_0 * n * T + phi + v );
x_p = w_0 * n * T + phi + v;
x = A * exp( log_x );
x_u = unwrap( x_p );


figure
subplot(2,1,1)
plot(n, x);
title(['SNR: ',  num2str( snr_log ), ' db']);
legend('x')
subplot(2,1,2)
plot(n, unwrap(angle(x)))
legend('unwrap(x)')
