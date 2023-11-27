tic
Y([1 2],:) = f([1 2],:).*100;
Y([3 4 5 6 7],:) = log10(f([3 4 5 6 7],:));

X = Parameter;

Out_data = net(X);
perf = perform(net,Y,Out_data);
toc