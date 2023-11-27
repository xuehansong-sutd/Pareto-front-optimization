clear
clc
close all
set(0,'DefaultAxesTitleFontWeight','normal');

%% Load Simulated Data
load Parallel_sol_10k.mat
load Pareto_i_ph.mat
load Exp_i_ph.mat
load Exp_iR.mat
 
% Symbol size
MS = 8;
LW = 2;

% Parallel plot font size
p_FS = 20;
FS = 16;

% Color
green = [0,0.39,0];
blue = [0,0,0.55];
orange = [1,0.6,0.2];
grey = [0.75,0.75,0.75];
purple = [0.4,0,0.8];
red = [0.65,0.16,0.16];
black = [0,0,0];
pink = [1 0 1];

%% Load Experimental Data
Exp_CIS_stdal_Eff = 17.54;

Exp_TC_Eff_420 = 20.3;
Exp_BC_Eff_420 = 8.5;
Exp_Eff_420 = [Exp_TC_Eff_420 Exp_BC_Eff_420];
Exp_tand_Eff_420 = sum(Exp_Eff_420);
Exp_l_420 = [110 10 420 23 110];
Exp_T_420 = 59;

Exp_TC_Eff_700 = 19.7;
Exp_BC_Eff_700 = 7.7;
Exp_Eff_700 = [Exp_TC_Eff_700 Exp_BC_Eff_700];
Exp_tand_Eff_700 = sum(Exp_Eff_700);
Exp_l_700 = [110 10 700 23 110];
Exp_T_700 = 57;

Exp_TC_Eff_550 = 21.5;
Exp_BC_Eff_550 = 8.0;
Exp_tot_Eff_550 = Exp_TC_Eff_550+Exp_BC_Eff_550;
Exp_Eff_550 = [Exp_TC_Eff_550 Exp_BC_Eff_550];
Exp_l_550 = [110 10 550 23 110];
Exp_T_550 = 59;

Sim_TC_Eff_420 = 19.9;
Sim_BC_Eff_420 = 8.6;
Sim_Eff_420 = [Sim_TC_Eff_420 Sim_BC_Eff_420];
Sim_tand_Eff_420 = sum(Sim_Eff_420);
Sim_l_420 = [110 10 420 23 110];
Sim_T_420 = 63;
Sim_iT_420 = 204;
Sim_iphT_420 = i_ph_exp(1)*10;
Sim_iB_420 = 208.8;
Sim_iphB_420 = 208.8;
Sim_cost_420 = 8.39;

Sim_TC_Eff_550 = 21.2;
Sim_BC_Eff_550 = 8.2;
Sim_Eff_550 = [Sim_TC_Eff_550 Sim_BC_Eff_550];
Sim_tand_Eff_550 = sum(Sim_Eff_550);
Sim_l_550 = [110 10 550 23 110];
Sim_T_550 = 61;
Sim_iT_550 = 216;
Sim_iphT_550 = i_ph_exp(2);
Sim_iB_550 = 197.2;
Sim_iphB_550 = 197.3;
Sim_cost_550 = 8.74;

Sim_TC_Eff_700 = 20.0;
Sim_BC_Eff_700 = 7.9;
Sim_Eff_700 = [Sim_TC_Eff_700 Sim_BC_Eff_700];
Sim_tand_Eff_700 = sum(Sim_Eff_700);
Sim_l_700 = [110 10 700 23 110];
Sim_iT_700 = 225.6;
Sim_iphT_700 = i_ph_exp(3);
Sim_iB_700 = 189.8;
Sim_iphB_700 = 189.95;
Sim_cost_700 = 9.15;

% Number of Data points
n_pts = length(solution2);

% Tandem Device efficiency
Tandem_Eff = -(objectiveValue2(:,1)+objectiveValue2(:,2)).*100;
TC_Eff = -objectiveValue2(:,1).*100;
BC_Eff = -objectiveValue2(:,2).*100;
l_tot = objectiveValue2(:,3).*1e9;
Mat_cost = [0.005177497, 0.001434301, 0.002744, 0.264, 0.005177497];
Mat_cost_Pareto = solution2.*Mat_cost.*1e9;
% create array with solution and objective values and tandem_eff
solnObjTand = [solution2.*1e9, l_tot,sum(Mat_cost_Pareto,2),TC_Eff,BC_Eff,Tandem_Eff];

% sort by tandem efficiency
sortedSolnObjTand = sortrows(solnObjTand, 10);
Opt_Tand = sortedSolnObjTand(end,:);
scaled_Tand = rescale(sortedSolnObjTand(:,end));
hsv = rgb2hsv(parula);
cm_data=interp1(linspace(0,1,size(jet,1)),hsv,scaled_Tand);
cm_data=hsv2rgb(cm_data);

% %% Fourth subplot
% figure(1)
% subplot(1,3,[1 2 3])
% title('(a)')
% % tick labels for parallel plot
% Ticklabels = {'F','H','P','E','B','l tot','Mat C','Eff Perov','Eff CIGS', 'Eff tot'};
% 
% % create array with base and parallel plot data
% % base = [Exp_l_420,sum(Exp_l_420), Exp_Eff_420, Exp_Eff_420(1) + Exp_Eff_420(2)];
% % Exp_550 = [Exp_l_550,sum(Exp_l_550), Exp_Eff_550, Exp_Eff_550(1) + Exp_Eff_550(2)];
% % Exp_700 = [Exp_l_700,sum(Exp_l_700), Exp_Eff_700, Exp_Eff_700(1) + Exp_Eff_700(2)];
% base = [Sim_l_420,sum(Sim_l_420), Sim_cost_420, Sim_Eff_420, Sim_Eff_420(1) + Sim_Eff_420(2)];
% Exp_550 = [Sim_l_550,sum(Sim_l_550),Sim_cost_550 Sim_Eff_550, Sim_Eff_550(1) + Sim_Eff_550(2)];
% Exp_700 = [Sim_l_700,sum(Sim_l_700),Sim_cost_700, Sim_Eff_700, Sim_Eff_700(1) + Sim_Eff_700(2)];
% 
% parPlotData = [sortedSolnObjTand;Opt_Tand; base;Exp_550;Exp_700];
% % create parallel plot
% p = parallelplot(parPlotData(:,1:10));
% 
% % group each data individually
% p.GroupData = 1:n_pts+4;
% 
% % set colour and base data to black
% p.Color = [cm_data;pink;black;red;purple];
% LW = ones(1,n_pts);
% LW(1,end+1:end+4) = [3;3;3;3];
% p.LineWidth = LW;
% 
% Linestyle = cell(1,n_pts+4);
% for i = 1:n_pts+4
%     Linestyle(1,i) = {'-'};
% end
% for i = n_pts+2:n_pts+4
%     Linestyle(1,i) = {'-'};
% end
% p.LineStyle = Linestyle;
% 
% % hide legend and format axes
% p.LegendVisible = 'off';
% p.FontName = 'Times New Roman';
% p.CoordinateTickLabels = Ticklabels;
% p.FontSize = p_FS;
% 
% %% Fifth subplot
% figure(2)
% subplot(1,3,[1 2 3])
% solnObjTand = [solution2.*1e9,i_R(:,[4 5])./i_ph(1,:)',sum(i_R,2)./i_ph(1,:)',i_ph(1,:)',i_ph(2,:)', Tandem_Eff];
% sortedSolnObjTand = sortrows(solnObjTand, 11);
% sortedSolnObjTand = sortedSolnObjTand(:,1:10);
% Opt_Tand = sortedSolnObjTand(end,:);
% title('(e)')
% % tick labels for parallel plot
% Ticklabels = {'F','H','P','E','B','RII','RIII','Rtot','i_T','i_B'};
% 
% % create array with base and parallel plot data
% base = [Sim_l_420,i_R_Exp(1,[4 5])./i_ph_exp(1,1),sum(i_R_Exp(1,:))./i_ph_exp(1,1),Sim_iT_420,Sim_iB_420];
% Exp_550 = [Sim_l_550,i_R_Exp(2,[4 5])./i_ph_exp(1,2),sum(i_R_Exp(2,:))./i_ph_exp(1,2),Sim_iT_550,Sim_iB_550];
% Exp_700 = [Sim_l_700,i_R_Exp(1,[4 5])./i_ph_exp(1,3),sum(i_R_Exp(3,:))./i_ph_exp(1,3),Sim_iT_700,Sim_iB_700];
% 
% parPlotData = [sortedSolnObjTand;Opt_Tand; base;Exp_550;Exp_700];
% % create parallel plot
% p = parallelplot(parPlotData(:,1:10));
% 
% % group each data individually
% p.GroupData = 1:n_pts+4;
% 
% % set colour and base data to black
% p.Color = [cm_data;pink;black;red;purple];
% LW = ones(1,n_pts);
% LW(1,end+1:end+4) = [3;3;3;3];
% p.LineWidth = LW;
% 
% Linestyle = cell(1,n_pts+4);
% for i = 1:n_pts+4
%     Linestyle(1,i) = {'-'};
% end
% for i = n_pts+2:n_pts+4
%     Linestyle(1,i) = {'-'};
% end
% p.LineStyle = Linestyle;
% 
% % hide legend and format axes
% p.LegendVisible = 'off';
% p.FontName = 'Times New Roman';
% p.CoordinateTickLabels = Ticklabels;
% p.FontSize = p_FS;

%% Combined subplot
figure(3)
subplot(1,3,[1 2 3])
solnObjTand = [solution2.*1e9,l_tot,sum(Mat_cost_Pareto,2),i_R(:,[4 5])./i_ph(1,:)',sum(i_R,2)./i_ph(1,:)',i_ph(1,:)',i_ph(2,:)',TC_Eff,BC_Eff, Tandem_Eff];
sortedSolnObjTand = sortrows(solnObjTand, 15);
sortedSolnObjTand = sortedSolnObjTand(:,1:15);
Opt_Tand = sortedSolnObjTand(end,:);
title('(e)')
% tick labels for parallel plot
Ticklabels = {'F','H','P','E','B','ltot','MatC','II','III','tot','i_T','i_B','EffT','EffB','Efftot'};

% create array with base and parallel plot data
base = [Sim_l_420,sum(Sim_l_420), Sim_cost_420,i_R_Exp(1,[4 5])./i_ph_exp(1,1)',sum(i_R_Exp(1,:))./i_ph_exp(1,1)',i_ph_exp(1,1)',i_ph_exp(2,1)',Sim_Eff_420, Sim_Eff_420(1) + Sim_Eff_420(2)];
Exp_550 = [Sim_l_550,sum(Sim_l_550),Sim_cost_550,i_R_Exp(2,[4 5])./i_ph_exp(1,2)',sum(i_R_Exp(2,:))./i_ph_exp(1,2)',i_ph_exp(1,2)',i_ph_exp(2,2)',Sim_Eff_550, Sim_Eff_550(1) + Sim_Eff_550(2)];
Exp_700 = [Sim_l_700,sum(Sim_l_700),Sim_cost_700,i_R_Exp(3,[4 5])./i_ph_exp(1,3)',sum(i_R_Exp(3,:))./i_ph_exp(1,3)',i_ph_exp(1,3)',i_ph_exp(2,3)',Sim_Eff_700, Sim_Eff_700(1) + Sim_Eff_700(2)];

parPlotData = [sortedSolnObjTand;Opt_Tand; base;Exp_550;Exp_700];
% create parallel plot
p = parallelplot(parPlotData(:,1:15));

% group each data individually
p.GroupData = 1:n_pts+4;

% set colour and base data to black
p.Color = [cm_data;pink;black;red;purple];
LW = ones(1,n_pts);
LW(1,end+1:end+4) = [4;3;3;3];
p.LineWidth = LW;

Linestyle = cell(1,n_pts+4);
for i = 1:n_pts+4
    Linestyle(1,i) = {'-'};
end
for i = n_pts+2:n_pts+4
    Linestyle(1,i) = {'--'};
end
p.LineStyle = Linestyle;

% hide legend and format axes
p.LegendVisible = 'off';
p.FontName = 'Times New Roman';
p.CoordinateTickLabels = Ticklabels;
p.FontSize = p_FS;