F_s = 10e6;
T = 10e-6;
f_0 = 10e5;
w_0 = 2 * pi * f_0;
phi = pi / 8;
A = 1;
N = 513;
n_0 = -256;
n = -256:1:N-257;

H = [T * n', ones(length(n), 1)];
C = eye(N);
P = N * ( N - 1) / 2;
Q = N * ( N - 1) * (2*N - 1) / 6;

snr_db = 30;
s = sigma(snr_db);
%s = A^2 / 2 ./ db2mag(snr_db);
mean = 0;

[h, w] = size(snr_db);
var_w_hat = zeros(h, w);
var_phi_hat = zeros(h, w);
v = zeros(N, w);
snr = zeros(h, w);
snr_log = zeros(h, w);
snr_mag = zeros(h, w);
x = zeros(N, w);
log_x = zeros(N, w);
x_p = zeros(N, w);
x_u = zeros(N, w);
C_hat = zeros(N, w);
theta_hat = zeros(N, w);


for i=1:length(snr_db)
    var_w_hat = 12*s(i)^2 / (A^2 * T^2 * N * ( N^2 - 1));
    var_phi_hat = 12*s(i)^2 * (n_0^2 * N + 2*n_0 * P + Q) / (A^2 * N^2 * ( N^2 - 1));
    v(i, :) = s(i) * randn(N, 1).' + mean;
    snr(i) = 1 / (2 * s(i)^2);
    snr_log(i) = 10 * log10(snr(i));
    snr_mag(i) = db2mag(snr_log(i));
    
    log_x(i, :) = 1i * ( w_0 * n * T + phi + v(i, :));
    x_p(i, :) = w_0 * n * T + phi + v(i, :);
    x(i, :) = A * exp( log_x(i, :));
    x_u(i, :) = unwrap( x_p(i, :));

    C_hat(i, :) = s(i)^2 * (H' / (A * A' ) * H);
    theta_hat(i, :) = (H' / ( A * A' ) * H ) * H' / ( A * A' ) .* x(:, i);
end



figure
subplot(4,1,1)
plot(n, x);
title(['SNR: ',  num2str( snr_log ), ' db']);
legend('x')
subplot(4,1,2)
plot(n, unwrap(angle(x)))
legend('unwrap(x)')
subplot(4,1,3)
plot(n, var(theta_hat(1,:)))
title(['var: ',  num2str( var(theta_hat(1,:)))]);
legend('omega^')
subplot(4,1,4)
plot(n, var(theta_hat(2,:)))
title(['var: ',  num2str( var(theta_hat(2,:)))]);
legend('phi^')
