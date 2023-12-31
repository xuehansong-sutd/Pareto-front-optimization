%% 
%   DESIGN OF SIMULATIONS (Random Sampling)
%
%
%
% Created: 5 May 2021

%%
clear
close all
clc
format short

%% Bounds for data set
L1 = [100*1e-9 350*1e-9]; % TC F/nm
L2 = [5*1e-9 60*1e-9]; % TC E/nm
L3 = [300*1e-9 1000*1e-9]; % TC P/nm
L4 = [5*1e-9 30*1e-9]; % TC P/nm
L6 = [100*1e-9 350*1e-9]; % TC B/nm  

%% 2nd data set 
L1_d = L1(2)-L1(1);
L2_d = L2(2)-L2(1); 
L3_d = L3(2)-L3(1); 
L4_d = L4(2)-L4(1); 
L6_d = L6(2)-L6(1); 

%% Create Random Sample
rng default
Num_Runs = 1000;
Num_Var = 5;
X1 = rand(Num_Runs,1);
X2 = rand(Num_Runs,1);
X3 = rand(Num_Runs,1);
X4 = rand(Num_Runs,1);
X6 = rand(Num_Runs,1);

L1_RAND = zeros(Num_Runs,1);
L2_RAND = zeros(Num_Runs,1);
L3_RAND = zeros(Num_Runs,1);
L4_RAND = zeros(Num_Runs,1);
L6_RAND = zeros(Num_Runs,1);


for i = 1:Num_Runs
    
L1_RAND(i,1) = L1(1)+X1(i,1)*L1_d;
L2_RAND(i,1) = L2(1)+X2(i,1)*L2_d;
L3_RAND(i,1) = L3(1)+X3(i,1)*L3_d;
L4_RAND(i,1) = L4(1)+X4(i,1)*L4_d;
L6_RAND(i,1) = L6(1)+X6(i,1)*L6_d;

end

Parameter = [L1_RAND L2_RAND L3_RAND L4_RAND L6_RAND];

%% Experimental point
L1E = 110*1e-9;
L2E = 10*1e-9;
L3E = 420*1e-9;
L4E = 23*1e-9;
L6E = 110*1e-9;

Parameter(end+1,:) = [L1E L2E L3E L4E L6E];

%% Specifying parametric space and Experimental point
Categorical = cell(Num_Runs+1,1);
Categorical(1:end-1,1) = {'Parametric Space'};
Categorical(end,1) = {'Experimental Point'};

%% Specifying parametric space and Experimental point
Categorical = cell(Num_Runs+1,1);
Categorical(1:end-1,1) = {'Parametric Space'};
Categorical(end,1) = {'Experimental Point'};

%% Figure 1: Combined parameter space (Parallel plot)

ax = gca;
% ax.YAxis.Exponent = 2;
p = parallelplot(Parameter);
% ax = gca;
% ax.CoordinateTickLabels = 'tex';
% Ticklabels = ["L1 / m","L3 / m","L6 / m","L7 / m","L8 / m",...
%                 "L9 / m","j_T / A m^-2","j_B / A m^-2","j_tot / A m^-2"];
Ticklabels = ["","","","",...
                "",""];

p.CoordinateTickLabels = Ticklabels;
p.GroupData = Categorical;
Green = [153/255,255/255,153/255];
Blue = [153/255,204/255,255/255];
Red = [255/255,153/255,153/255];
Black = [0,0,0];
Yellow = [255/255,255/255,102/255];

p.Color = [Red;Black]
p.LineWidth = [1,3];
p.LineStyle = {'--';':'};
p.FontName = 'Times New Roman';
p.FontSize = 20;

%% Figure 2
createfigure(Parameter,Categorical)

%% Export Parameter file
Parameter = Parameter';
save Parameters.mat Parameter

%%
function createfigure(data1, GroupData1)
%CREATEFIGURE(data1, GroupData1)
%  DATA1:  data
%  GROUPDATA1:  groupdata

%  Auto-generated by MATLAB on 29-Mar-2021 13:27:04

% Create figure
figure1 = figure('WindowState','maximized');

Green = [153/255,255/255,153/255];
Blue = [153/255,204/255,255/255];
Red = [255/255,153/255,153/255];
Black = [0,0,0];
Yellow = [255/255,255/255,102/255];

FontSize = 24

% Create parallelplot
parallelplot(figure1,data1,'GroupData',GroupData1,'Jitter',0.1,...
    'LineStyle',{'--',':'},...
    'LineWidth',[1 3],...
    'DataNormalization','range',...
    'CoordinateTickLabels',{'','','','','',''},...
    'FontName','Times New Roman',...
    'FontSize',FontSize,...
    'Color',[Red;Black]);

% Create textbox
annotation(figure1,'textbox',...
    [0.799374999999999 0.0471687886130276 0.0566666678935289 0.0601114661640423],...
    'String',{'{\itl}^{P}_{B} / nm'},...
    'LineStyle','none',...
    'FontSize',FontSize,...
    'FontName','Times New Roman');

% Create textbox
annotation(figure1,'textbox',...
    [0.672916666666666 0.0467706994410532 0.0575000012467305 0.0601114661640423],...
    'String',{'{\itl}^{H}_{B} / nm'},...
    'LineStyle','none',...
    'FontSize',FontSize,...
    'FontName','Times New Roman');

% Create textbox
annotation(figure1,'textbox',...
    [0.55125 0.0467706994410532 0.0566666678935289 0.0601114661640423],...
    'String',{'{\itl}^{F}_{B} / nm'},...
    'LineStyle','none',...
    'FontSize',FontSize,...
    'FontName','Times New Roman');

% Create textbox
annotation(figure1,'textbox',...
    [0.425625 0.0479649669569768 0.056666667893529 0.0601114661640423],...
    'String',{'{\itl}^{B}_{T} / nm'},...
    'LineStyle','none',...
    'FontSize',FontSize,...
    'FontName','Times New Roman');

% Create textbox
annotation(figure1,'textbox',...
    [0.297708333333333 0.0479649669569768 0.0562500012169282 0.0601114661640423],...
    'String',{'{\itl}^{P}_{T} / nm'},...
    'LineStyle','none',...
    'FontSize',FontSize,...
    'FontName','Times New Roman');

% Create textbox
annotation(figure1,'textbox',...
    [0.177083333333333 0.0475668777850023 0.0562500012169282 0.0601114661640423],...
    'String',{'{\itl}^{F}_{T} / nm'},...
    'LineStyle','none',...
    'FontSize',FontSize,...
    'FontName','Times New Roman');
end