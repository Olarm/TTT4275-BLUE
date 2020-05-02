
f1 = figure(1);
subplot(2,1,1)
plot(SNR, CRLB_omega, SNR, BLUE_omega)
tit = title({'Variance of $\hat{\omega}$ for CRLB and BLUE'});
set(tit, 'interpreter', 'latex');
grid on;
xlabel('SNR [dB]');
leg = legend({'CRLB var ($\hat{\omega}$)', 'BLUE var ($\hat{\omega}$)'});
set(leg, 'interpreter', 'latex');

subplot(2,1,2)
plot(SNR, BLUE_omega - CRLB_omega)
tit = title({'Difference in variance of $\hat{\omega}$ for CRLB and BLUE'});
set(tit, 'interpreter', 'latex');
grid on;
xlabel('SNR [dB]');

saveas(f1, 'figures/var_est_omega', 'epsc');

f2 = figure(2);

subplot(2,1,1)
plot(SNR, CRLB_phi, SNR, BLUE_phi)
tit = title({'Variance of $\hat{\phi}$ for CRLB and BLUE'});
set(tit, 'interpreter', 'latex');
grid on;
xlabel('SNR [dB]');
leg = legend({'CRLB var ($\hat{\phi}$)', 'BLUE var ($\hat{\phi}$)'});
set(leg, 'interpreter', 'latex');

subplot(2,1,2)
plot(SNR, BLUE_phi - CRLB_phi)
tit = title({'Difference in variance of $\hat{\phi}$ for CRLB and BLUE'});
set(tit, 'interpreter', 'latex');
grid on;
xlabel('SNR [dB]');

saveas(f2, 'figures/var_est_phi', 'epsc');

% f1 = figure(1);
% grid on;
% title('Variance of estimated omega for CRLB and BLUE');
% xlabel('SNR [dB]');
% ylabel('Variance of the estimate of omega from CRLB and BLUE');
% plot(SNR, CRLB_omega, SNR, BLUE_omega);
% legend('CRLB omega', 'BLUE omega');
% 
% saveas(f1, 'figures/var_est_omega', 'epsc');

% f2 = figure(2);
% grid on;
% title('Difference between the variance of the estimated omega');
% xlabel('SNR [dB]');
% ylabel('Difference of the variances of the estimate of omega from CRLB and BLUE');
% plot(SNR, BLUE_omega - CRLB_omega);
% legend('Omega: Difference between CRLB and BLUE');
% 
% saveas(f2, 'figures/diff_var_est_omega', 'epsc');

% f3 = figure(3);
% grid on;
% title('Plot of the variance of the estimated phi given by CRLB and BLUE');
% xlabel('SNR [dB]');
% ylabel('Variance of the estimate of phi from CRLB and BLUE');
% plot(SNR, CRLB_phi, SNR, BLUE_phi);
% legend('CRLB phi', 'BLUE phi');
% 
% saveas(f3, 'figures/var_est_phi', 'epsc');
% 
% f4 = figure(4);
% grid on;
% title('Difference between the variance of the estimated phi');
% xlabel('SNR [dB]');
% ylabel('Difference of the variances of the estimate of phi from CRLB and BLUE');
% plot(SNR, BLUE_phi - CRLB_phi);
% legend('Phi: Difference between CRLB and BLUE');
% 
% saveas(f4, 'figures/diff_var_est_phi', 'epsc');

figure(5);
grid on;
for i = 1:length(SNR)
   sig = A * exp(1j * (est(i, 1) * n' * T + mod(est(i, 2), pi)));
   y = unwrap(angle(sig));
   plot(n, y)
   if i == 1
       hold on;
   end
end
hold off;
Legend = cell(length(SNR), 1);
for i = 1:length(SNR)
    Legend{i} = strcat('Signal to Noise Ratio: ', num2str(SNR(i)));
end
legend(Legend);

f6 = figure(6);
title('Unwrapped angle of signals');
ylabel('Angle [rad]');
xlabel('SNR [dB]');
grid on;
for i = 1:length(SNR)
    sig = A * exp(1j * (est(i, 1) * n' * T + mod(est(i, 2), pi)));
    plot(n, unwrap(angle(x(i, :))));
    if i == 1
        hold on;
    end
end
hold off;

legend('-10 db SNR', '0 db SNR', '10 db SNR', '20 db SNR', '30 db SNR', '40 db SNR')
saveas(f6, 'figures/unwrapped_angles', 'epsc');

