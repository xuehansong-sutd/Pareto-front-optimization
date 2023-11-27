% Load NN model into MATLAB
load 66logsig.mat
% Load Output and Input for test data
load Test_Output.mat
load Test_Parameters.mat

Predict_out = net(X);
figure(1)
scatter(Y,Predict_out);
hold on
plot(Y,Y);