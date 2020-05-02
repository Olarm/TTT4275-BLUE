F_s = 10e6;
T = 10e-6;
f_0 = 10e5;
w_0 = 2 * pi * f_0;
phi = pi / 8;
A = 1;
N = 513;
n_0 = -256;
n = -256:1:N-257;
n = n';

H = [T * n, ones(length(n), 1)];
C = eye(N);
P = N * ( N - 1) / 2;
Q = N * ( N - 1) * (2*N - 1) / 6;

snr_db = [-10, 0, 10, 20, 30, 40];
s = sigma(snr_db);
var = (A^2 / 2) ./ db2mag(snr_db);
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
theta_hat = zeros(2, w);
martin = zeros(N, w);

for i=1:length(snr_db)
    v(:, i) = s(i) * randn(N, 1).' + mean;
    snr(i) = 1 / (2 * s(i)^2);
    snr_log(i) = 10 * log10(snr(i));
    snr_mag(i) = db2mag(snr_log(i));
    
    martin(:,i) = gen_signal_martin(w_0, n, A, T, phi, 0, sqrt(var(i)));
    
    log_x(:, i) = 1i * ( w_0 * n * T + phi + v(:, i));
    x(:, i) = A * exp( log_x(:, i));
    x_u(:, i) = unwrap(angle(x(:, i)));
    
    %C_hat(:, i) = s(i)^2 * (H' / (A * A' ) * H);
    theta_hat(:, i) = (H' / ( A * A' ) * H ) * H' / ( A * A' ) * martin(:, i);
end

%%% B

blue = BLUE_white(martin, H, C);
blue_var = inv(H' * inv(

CRLB_w_hat = 12.*s.^2 / (A^2 * T^2 * N * ( N^2 - 1));
CRLB_phi_hat = 12.*s.^2 * (n_0^2 * N + 2*n_0 * P + Q) / (A^2 * N^2 * ( N^2 - 1));

%%% C

phases = angle(martin(2:N,:)) - angle(martin(1:N-1,:));

%blue = BLUE_color(phases, H, C);


%%%% Oppg a

% figure
% plot(n, martin)
% legend('-10 db', '0 db', '10 db', '20 db', '30 db', '40 db')
% 
% figure;
% plot(n, unwrap(angle(martin)));
% grid on;
% title('Unwrapped angle of signals');
% ylabel('Angle [rad]');
% xlabel('SNR [dB]');
% legend('-10 db', '0 db', '10 db', '20 db', '30 db', '40 db')
% 
% 
% figure
% subplot(4,1,1)
% plot(n, x);
% title(['SNR: ',  num2str( snr_log ), ' db']);
% legend('x')
% subplot(4,1,2)
% plot(n, unwrap(angle(x)))
% legend('unwrap(x)')
% subplot(4,1,3)
% plot(n, var(theta_hat(1,:)))
% title(['var: ',  num2str( var(theta_hat(1,:)))]);
% legend('omega^')
% subplot(4,1,4)
% plot(n, var(theta_hat(2,:)))
% title(['var: ',  num2str( var(theta_hat(2,:)))]);
% legend('phi^')

%%% Oppg b
figure
subplot(2,2,1)
plot(snr_db, blue(2,:), snr_db, CRLB_phi_hat)
title('Phase estimate compared with CRLB')
grid on;
xlabel('SNR [dB]');
leg = legend({'BLUE var ($\hat{\phi}$)', 'CRLB var ($\phi$)'});
set(leg, 'interpreter', 'latex');

subplot(2,2,2)
plot(snr_db, blue(1,:), snr_db, CRLB_w_hat)
title('Frequency estimate compared with CRLB')
grid on;
xlabel('SNR [dB]');
leg = legend({'BLUE var ($\hat{\omega}$)', 'CRLB var ($\omega$)'});
set(leg, 'interpreter', 'latex');

subplot(2,2,3)
plot(snr_db, abs(blue(2,:) - CRLB_phi_hat))
title('Difference(BLUE phi, CRLB phi)')
grid on;
xlabel('SNR [dB]');

subplot(2,2,4)
plot(snr_db, abs(blue(1,:) - CRLB_w_hat))
title('Difference(BLUE omega, CRLB omega)')
grid on;
xlabel('SNR [dB]');


%%%% oppg c

% figure()
% plot(n(1:N-1), phases(:,1))

%%Function definitions
function x = gen_signal_martin(omega_0, n, A, T, phi, w_e, w_sigma)

w = gen_noise_martin(length(n), w_e, w_sigma);
x = A * exp(1j * (omega_0 * n * T + phi)) + w;

end

function v = gen_noise_martin(N, w_e, w_sigma)
    v = (normrnd(w_e, w_sigma, 1, N) + 1j * normrnd(w_e, w_sigma, 1, N))';
end

function est = BLUE_white(signal, H, C)
    M = ((H' * (C \ H)) \ H') / C;
    est = M * unwrap(angle(signal));
end

function est = BLUE_color(signal, H, C)
    M = ((H' * (C \ H)) \ H') / C;
    est = M(:, 1:512) * signal;
end
