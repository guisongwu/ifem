
figure(3);
clf(3);

mplot = 2;
nplot = 2;

sft_ext = @(x) sft(extend_mid(x));

d2Xidm2 = d2Xidm2_NW;

subplot(mplot, nplot, 1)
plot(sft(xbot), sft_ext(dXidm(:,1)), '-o', ...
     sft(xbot), sft_ext(dXidm_FD(:,1)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
legend('dXidm', 'dXidm FD');



subplot(mplot, nplot, 2)
icol = 1;
plot(sft(xbot), sft_ext(d2Xidm2(:,icol)), '-o', ...
     sft(xbot), sft_ext(d2Xidm2_FD(:,icol)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
legend('d2Xidm2', 'd2Xidm2 FD');


subplot(mplot, nplot, 3)
icol = 2;
plot(sft(xbot), sft_ext(d2Xidm2(:,icol)), '-o', ...
     sft(xbot), sft_ext(d2Xidm2_FD(:,icol)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
legend('d2Xidm2', 'd2Xidm2 FD');

subplot(mplot, nplot, 4)
icol = 10;
plot(sft(xbot), sft_ext(d2Xidm2(:,icol)), '-o', ...
     sft(xbot), sft_ext(d2Xidm2_FD(:,icol)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
legend('d2Xidm2 10', 'd2Xidm2 FD');


% subplot(mplot, nplot, 5)
% icol = 11;
% plot(sft(xbot), sft_ext(d2Xidm2(:,icol)), '-o', ...
%      sft(xbot), sft_ext(d2Xidm2_FD(:,icol)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
% legend('d2Xidm2 11', 'd2Xidm2 FD');



% subplot(mplot, nplot, 6)
% icol = 19;
% plot(sft(xbot), sft_ext(d2Xidm2(:,icol)), '-o', ...
%      sft(xbot), sft_ext(d2Xidm2_FD(:,icol)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
% legend('d2Xidm2 19', 'd2Xidm2 FD');


% subplot(mplot, nplot, 7)
% icol = 20;
% plot(sft(xbot), sft_ext(d2Xidm2(:,icol)), '-o', ...
%      sft(xbot), sft_ext(d2Xidm2_FD(:,icol)), '-x', 'LineWidth', 2, 'MarkerSize', 5);
% legend('d2Xidm2 20', 'd2Xidm2 FD');

