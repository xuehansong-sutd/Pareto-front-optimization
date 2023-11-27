clear
clc
close all
% Load NN model into MATLAB
% load 40radbas_39radbas.mat
load 39radbas_29tansig.mat
% load 38radbas_37radbas.mat
% load 40tansig_39radbasn_expanded.mat
% load 38logsig_39logsig.mat
% Load Output and Input for test data
load Test_Output.mat
load Test_Parameters.mat
set(0,'DefaultAxesTitleFontWeight','normal');

FS = 20;

Predict_out = net(Parameter);
figure(1)

%% TC_Eff
subplot(1,2,1)
scatter(Predict_out(1,:),f(1,:).*100);
hold on
plot(f(1,:).*100,f(1,:).*100);
v=get(1,'currentaxes');
title('(a)')
xlabel('Predicted \eta_{T}')
ylabel('Actual \eta_{T}')
set(v,'fontsize',FS,'fontname','Times New Roman')
box on
axis square
axis tight
hold off
%% BC_Eff
subplot(1,2,2)
scatter(Predict_out(2,:),f(2,:).*100);
hold on
plot(f(2,:).*100,f(2,:).*100);
v=get(1,'currentaxes');
title('(b)')
xlabel('Predicted \eta_{B}')
ylabel('Actual \eta_{B}')
set(v,'fontsize',FS,'fontname','Times New Roman')
box on
axis square
axis tight
hold off
% %% T
% subplot(1,3,3)
% scatter(Predict_out(3,:),f(4,:));
% hold on
% plot(f(4,:),f(4,:));
% v=get(1,'currentaxes');
% title('(c)')
% xlabel('Predicted $\bar{T}$','Interpreter','Latex')
% ylabel('Actual $\bar{T}$','Interpreter','Latex')
% set(v,'fontsize',FS,'fontname','Times New Roman')
% box on
% axis square
% axis tight
% hold off