%Import 3'-5'-AMP dataset
data = csvread("3-5-AMP.csv");

dataX = data(:, 1);
dataY = data(:, 2);
dataErr = data(:, 3);

%Initialize all known values
F6P = .1;       %mM
ATP = 2.3;      %mM
E1 = .12;      %uM
%E1 = .12 / 1000; %mM
KF6P = .11;     %mM
KATP = .42;     %mM
%kcat = .4;      %s-1
kcat = .4 * 60 * 60; %hr-1

r1 = kcat * E1 * (F6P / (KF6P + F6P)) * (ATP / (KATP + ATP)); %uM/hr

%Solve for W1 from AMP = 0 condition
W1 = dataY(1) / (r1 - dataY(1));
W1 = .0451058026219;

%Solve for W2fI value, assume fI <= 1 to obtain W2
W2fI = ((dataY/r1) + (dataY * W1 / r1) - W1) ./ (1 - (dataY / r1));
W2 = W2fI(6);

%Log normalize fI to generate Hill plot 
fI = W2fI ./ W2;
logTheta = log10(fI(2:5) ./ (1 - fI(2:5)));
logX = log10(dataX(2:5));

linX = log10((dataX(2):.01:dataX(5)));
lin = polyfit(logX, logTheta, 1);
linY = polyval(lin, linX);

close all
figure(1)
plot(logX, logTheta, 'o', linX, linY, '-');
legend('Experimental', 'Linear Fit');

%Initialize values from Hill plot
n = lin(1);
K = exp(-lin(2) / lin(1));

%Generate rate dataset from model using extracted constants
AMP = (0:.01:1).';
fI = (AMP ./ K).^n ./ (1 + (AMP ./ K).^n);
v = (W1 + W2*fI) ./ (1 + W1 + W2*fI);
modelY = r1 * v;

%Plot dataset and model with error bars
figure(2)
hold on
errorbar(dataX, dataY, dataErr);
plot(AMP, modelY);
legend('Experimental', 'Model');

