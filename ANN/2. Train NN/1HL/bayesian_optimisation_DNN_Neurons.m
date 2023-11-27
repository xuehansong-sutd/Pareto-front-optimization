clear
clc
close all
%% Prepare Data
load Parameters.mat
load Output.mat
Inputs = Parameter';
Targets = Output;
%% Setting up variables
hiddenLayerSizeRange_1 = [1 70];
optimVars = [optimizableVariable('Layer1Size',hiddenLayerSizeRange_1,'Type','integer')
             optimizableVariable('TrainFcn1',{'logsig' 'tansig' ...
             'compet' 'elliotsig' ...
             'netinv' 'radbas' 'radbasn' 'softmax' 'tribas'},'Type','categorical')];
%% Bayesian Optimization
ObjFcn = makeObjFcn(Inputs, Targets);
BayesObject = bayesopt(ObjFcn,optimVars,...
    'MaxObj',2000,...
    'MaxTime',24*60*60,...
    'IsObjectiveDeterministic',false,...
    'UseParallel',true);
%% Evaluate Final Network
bestIdx = BayesObject.IndexOfMinimumTrace(end);
fileName = BayesObject.UserDataTrace{bestIdx};
load(fileName);
YPredicted = net(Inputs);
testError = perform(net,Targets,YPredicted);
testError;
valError;

%% Objective Function For Optimisation
function ObjFcn = makeObjFcn(XTrain,YTrain)
    ObjFcn = @valErrorFun;
      function [valError,cons,fileName] = valErrorFun(optVars)
          % Solve an Input-Output Fitting problem with a Neural Network
          % Choose a Training Function
          trainFcn = 'trainbr';  % Bayesian Regularization backpropagation.
          % Create a Fitting Network
          layer1_size = optVars.Layer1Size;
          maxEpochs = 4000;
          TrainFcn1 = char(optVars.TrainFcn1);
%           TrainFcn1 = char(optVars.TrainFcn1);
%           TrainFcn2 = char(optVars.TrainFcn2);
          hiddenLayerSizes = layer1_size;
          % Specifying activation function at each layer
          net = fitnet(hiddenLayerSizes,trainFcn);

          % Setup Division of Data for Training, Validation, Testing
          net.divideParam.trainRatio = 100/100;
          net.divideParam.valRatio = 0/100;
          net.divideParam.testRatio = 0/100;
          % Activation function for hidden layers
          net.layers{1}.transferFcn = TrainFcn1;  % Hidden layer 1
          net.layers{2}.transferFcn = 'purelin';  % Output layer
          % performance function
          net.performFcn = 'mse';
          net.performParam.normalization = 'none';
          % Train the Network
          idx = randperm(10001);
          Train_idx = idx(1:10001);
          Test_idx = idx(9000:10001);
          net.trainParam.showWindow = true;
          net.trainParam.showCommandLine = false;
          net.trainParam.epochs = maxEpochs;
          [net,~] = train(net,XTrain(:,Train_idx),YTrain(:,Train_idx));
          % Test the Network
          YPredicted = net(XTrain(:,Test_idx));
          valError = perform(net,YTrain(:,Test_idx),YPredicted);
          fileName = num2str(valError) + ".mat";
          save(fileName,'net','valError')
          cons = [];
      end
  end