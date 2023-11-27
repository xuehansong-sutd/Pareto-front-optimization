clear
clc
close all
delete Layer_thickness.txt
delete EffTC_EffBC_T.txt
delete Layer_EffTEffB.txt

%% Load NN
% load 40radbas_39radbas.mat
load 39radbas_29tansig.mat

rng('default');
s = rng;

nvar = 5;

% Inputs = [l_F; l_H; l_P; l_E; l_B]; % loading input txt files as a row vector
lb = [100, 5, 420, 5,100].*1e-9; % l_F; l_H; l_P; l_E; l_B
ub = [350, 60, 700, 30,350].*1e-9; % ranges taken from standard grid search

% %% GA
% Set nondefault solver options
options = optimoptions('gamultiobj','CrossoverFraction',0.7,'PopulationSize',1e4,...
    'MaxTime',48*60*60,'UseParallel',true,'Display','iter','PlotFcn',...
    {'gaplotpareto', 'gaplotstopping'});

% Solve
[solution2,objectiveValue2] = gamultiobj(@(x)multiObjFcn(x,net),nvar,[],[],[],[],lb,ub,...
    [],[],options);

% Clear variables
clearvars options

%% Paretosearch
% options = optimoptions('paretosearch','ParetoSetSize',1000,...
%     'PlotFcn',{'psplotparetof' 'psplotparetox'},'Display','iter','MaxTime',48*60*60);
% % Solve
% [solution2,objectiveValue2] = paretosearch(@(x)multiObjFcn(x,net),nvar,[],[],[],[],lb,ub,...
%     [],options);
% 
% % Clear variables
% clearvars options

%%
function F = multiObjFcn(x,net)
tic
Output = net(x');
toc
Eff_TC = Output(1)./100;
Eff_BC = Output(2)./100;

Cost = 5.18e-3*(x(1)+x(5))*1e9 + 1.43e-3*x(2)*1e9 + 2.74e-3*x(3)*1e9 + 2.6e-1*x(4)*1e9;

F = [-Eff_TC; -Eff_BC; Cost];

%% Writing parametric search space to file

% Output_tot_str = '%e %e \n';
% fileID = fopen('EffTC_EffBC_T.txt','a');
% %     fprintf(fileID, 'Parameter space with current density\n\n');
% fprintf(fileID,Output_tot_str,[Eff_TC Eff_BC]);
% fclose(fileID);

% l_tot_str = '%e %e %e %e %e \n';
% Layer_thickness = x;
% fileID = fopen('Layer_thickness.txt','a');
% %     fprintf(fileID, 'Parameter space with current density\n\n');
% fprintf(fileID,l_tot_str,Layer_thickness);
% fclose(fileID);
% 
% Output_tot_str = '%e %e \n';
% fileID = fopen('EffTC_EffBC_T.txt','a');
% %     fprintf(fileID, 'Parameter space with current density\n\n');
% fprintf(fileID,Output_tot_str,[Eff_TC Eff_BC]);
% fclose(fileID);

tot_str = '%e %e %e %e %e %e %e %e \n';
Layer_thickness = x;
fileID = fopen('Layer_EffTEffB_Cost.txt','a');
%     fprintf(fileID, 'Parameter space with current density\n\n');
fprintf(fileID,tot_str,[Layer_thickness Eff_TC Eff_BC Cost]);
fclose(fileID);
end