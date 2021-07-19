close all;
clear;

tvec = 1:10;

% Smax=8
MSE_genie_8 = load('MSE_vec_genie_8.mat','MSE_vec_avg').MSE_vec_avg;
MSE_kfcs_known_8 = load('MSE_vec_kfcs_known_8.mat','MSE_vec_avg').MSE_vec_avg;
MSE_kfcs_unknown_8 = load('MSE_vec_kfcs_unknown_8.mat','MSE_vec_avg').MSE_vec_avg;
MSE_regcs_8 = load('MSE_vec_regcs_8.mat','MSE_vec_avg').MSE_vec_avg;
figure;
plot(tvec,MSE_genie_8,'-.');
hold on;
plot(tvec,MSE_kfcs_known_8,'-o');
plot(tvec,MSE_kfcs_unknown_8, '-o');
% plot(tvec,MSE_regcs_8);

% legend('Genie-aided KF');
% legend('Genie-aided KF','KF-CS T1 known');
legend('Genie-aided KF','KF-CS T1 known','KF-CS T1 unknown');
% legend('Genie-aided KF','KF-CS T1 known','KF-CS T1 unknown','Regular CS');

title('S_{max}=8');
grid();
ylabel('MSE');
xlabel('Time');
ylim([0 0.4]);

% Smax=16
MSE_genie_16 = load('MSE_vec_genie_16.mat','MSE_vec_avg').MSE_vec_avg;
MSE_kfcs_known_16 = load('MSE_vec_kfcs_known_16.mat','MSE_vec_avg').MSE_vec_avg;
MSE_kfcs_unknown_16 = load('MSE_vec_kfcs_unknown_16.mat','MSE_vec_avg').MSE_vec_avg;
MSE_regcs_16 = load('MSE_vec_regcs_16.mat','MSE_vec_avg').MSE_vec_avg;
figure;
plot(tvec,MSE_genie_16,'-.');
hold on;
plot(tvec,MSE_kfcs_known_16,'-o');
plot(tvec,MSE_kfcs_unknown_16, '-o');
% plot(tvec,MSE_regcs_16);
grid();
ylabel('MSE');
xlabel('Time');
title('S_{max}=16');
ylim([0 5]);

% legend('Genie-aided KF');
% legend('Genie-aided KF','KF-CS T1 known');
legend('Genie-aided KF','KF-CS T1 known','KF-CS T1 unknown');
% legend('Genie-aided KF','KF-CS T1 known','KF-CS T1 unknown','Regular CS');


% Smax=25
MSE_genie_25 = load('MSE_vec_genie_25.mat','MSE_vec_avg').MSE_vec_avg;
MSE_kfcs_known_25 = load('MSE_vec_kfcs_known_25.mat','MSE_vec_avg').MSE_vec_avg;
MSE_kfcs_unknown_25 = load('MSE_vec_kfcs_unknown_25.mat','MSE_vec_avg').MSE_vec_avg;
MSE_regcs_25 = load('MSE_vec_regcs_25.mat','MSE_vec_avg').MSE_vec_avg;
figure;
plot(tvec,MSE_genie_25,'-.');
hold on;
plot(tvec,MSE_kfcs_known_25,'-o');
plot(tvec,MSE_kfcs_unknown_25, '-o');
% plot(tvec,MSE_regcs_25);
title('S_{max}=25');
grid();
ylabel('MSE');
xlabel('Time');
ylim([0 30]);
% legend('Genie-aided KF');
% legend('Genie-aided KF','KF-CS T1 known');
legend('Genie-aided KF','KF-CS T1 known','KF-CS T1 unknown');
% legend('Genie-aided KF','KF-CS T1 known','KF-CS T1 unknown','Regular CS');


% Full KF plot
MSE_fullkf_8 = load('MSE_vec_fullkf_8.mat','MSE_vec_avg').MSE_vec_avg;
MSE_fullkf_16 = load('MSE_vec_fullkf_16.mat','MSE_vec_avg').MSE_vec_avg;
MSE_fullkf_25 = load('MSE_vec_fullkf_25.mat','MSE_vec_avg').MSE_vec_avg;
figure;
plot(tvec,MSE_fullkf_8,'-d');
hold on;
plot(tvec,MSE_fullkf_16,'-o' );
plot(tvec,MSE_fullkf_25,'-s' );
grid();
ylabel('MSE');
xlabel('Time');
legend('S_{max}=8','S_{max}=16','S_{max}=25');
title('MSE of Full KF');
