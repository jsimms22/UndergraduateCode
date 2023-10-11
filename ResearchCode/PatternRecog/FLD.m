prompt = 'Which set of data should be loaded (input 1 or 2)? ';
nn = input(prompt);
if nn == 1
    load hmwk1-dataset1.mat;
    disp('Data set 1 loaded.');
else
    load hmwk1-dataset2.mat;
    disp('Data set 2 loaded.');
end

%Trellis plot for features in all classes
n=1000; %This is the number of data points in each class
xx = [x1,x2]; %This concatenates the two classes into one array
yy = [ones(1,n),2*ones(1,n)]; %This assigns class labels to the
                              %data examples
figure; 
gplotmatrix(xx',[],yy','br');%This draws the trellis plot

m1 = mean(x1')'; 
m2 = mean(x2')';

% meanAll = mean(xx');
% meanAll = meanAll';
% 
% for ii = 1:4
%    x1(ii,:) = x1(ii,:) - meanAll(ii,1);
%    x2(ii,:) = x2(ii,:) - meanAll(ii,1);
% end

cov1 = cov(x1') .* n; %Covariance is same as scatter matrix but times 1/(n-1)
cov2 = cov(x2') .* n; %Where n is the number of points in dimension Di
S_W = cov1 + cov2;  %Within-Class Scatter Matrix
InvS = inv(S_W);  %Inverse of the within-class scatter

M = (m1 - m2); %diag matrix of the difference between means
w1 = InvS * M; %vector that the data will be projected onto

S_B = (m1 - m2) * (m1 - m2)'; %between scatter matrix

% Jw = InvS * S_B;

[w, Jw] = eigsort(InvS * S_B);
% for ii = 1:4
%     w(:,ii) = w(:,ii) .* Jw(ii,1);
% end

y1 = w(:,1)' * x1; %projecting x1 features onto axis w
y2 = w(:,1)' * x2; %Projecting x2 features onto axis w

%overlapping histogram of y1 and y2 features
figure;
histogram(y1); %gaussian-like distribution of y1 features
hold on
histogram(y2); %gaussian-like distribution of y2 features
hold off

