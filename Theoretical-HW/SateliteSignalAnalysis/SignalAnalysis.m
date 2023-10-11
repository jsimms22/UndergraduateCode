load datafile
X = VarName1; %time stored in ms
Y = VarName2; %signal strength stored in mW
%uncertainty of y at time x
uncertainty = diag(VarName3.^(-2));

%to establish to what iteration the for loop runs to
prompt = 'What is your largest desired n? ';
n = input(prompt);

%A matrix for when N = 0
A = ones(59,1);
%To calculate the A matrix to find all coefficient
for count = 1:n
    A(:,count+1) = X.^(-count);
end
%c = coefficients
c = ((A' * uncertainty * A)^(-1)) * (A' * uncertainty * Y);
c;
%Uncertainties of coefficients
cUncert = ((A' * uncertainty * A)^(-1)).^(1/2);

%Theoretical y's calculated
yFit = A * c;
chi = 0;

%To calculate the chi for "goodness" of fit
for count = 1:59
    chi = chi + ((yFit(count,1) - Y(count,1))^2 / 59);
end
chi;

yTest = 0;

%Testing prediction for .2 using linear least fit
for count = 1:n+1
    yTest =  yTest + c(count,1) * .2^(1 - count);
end
yTest;

%BELOW HERE IS FOR MAKING COMPARISON GRAPHS!

x = 0:.2:15;
x = x';
A = ones(76,1);
for count = 1:n
    A(:,count+1) = x.^(-count);
end
fx = A * c;
hold off
errorbar(X, Y, VarName3, 'k.')
hold on
plot(x, fx, 'rd-')