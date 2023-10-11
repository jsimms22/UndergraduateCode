load hmwk3.mat;

figure; % plotting data for a visual interpretation
plot3(x(1,:),x(2,:),x(3,:),'k.')
grid on
title('Original Data')
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
hold on;

%%  PCA Via Traditional Eigen Problem A * x = llambda * x

% mean-correcting data)
x = x';
M = mean(x)';
x = x';
% M = [mean(x(1,:));mean(x(2,:));mean(x(3,:))]; % calculating mean for each individual feature vector
xCorrected = ones(3,100); % pre-constructing matrix for mean-corrected data

for ii=1:3
   xCorrected(ii,:) = x(ii,:) - M(ii,1); % manipulating original data to be mean-corrected
end

figure; % plotting data for a visual interpretation
plot3(xCorrected(1,:),xCorrected(2,:),xCorrected(3,:),'k.')
grid on
title('Mean Corrected Data - Eigen Problem')
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
hold on;

% finding eigenvalues and eigenvectors
covariance = cov(xCorrected'); % covariance matrix of mean-corrected data
[eigenVect,eigenVal] = eig(covariance); % solving for eigenvectors and eigenvalues

for ii=1:3
   eigenVect(:,ii) = eigenVect(:,ii) .* eigenVal(ii,ii); % multiplying eigenvectors by their eigenvalues for graph analysis below
end

eigenVal = diag(eigenVal);              % extract diagonal elements of eigenvalue matrix
[eigenVal,index] = esort(eigenVal);     % put magnitude of eigenvalues in descending order
eigenVect(:,:) = eigenVect(:,index);  % sort eigenvectors likewise

% overlaying eigenvectors on plot of original data
plot3([0 eigenVect(1,1)],[0 eigenVect(2,1)],[0 eigenVect(3,1)],'r-');   % overlaying first eigenvector onto original data plot
plot3([0 eigenVect(1,2)],[0 eigenVect(2,2)],[0 eigenVect(3,2)],'g-');   % overlaying second eigenvector onto original data plot
plot3([0 eigenVect(1,3)],[0 eigenVect(2,3)],[0 eigenVect(3,3)],'b-');   % overlaying third eigenvector onto original data plot
hold off;

% calculating and graphing percent variance accounted for(%VAR); for both original data and PCs
totalVariance = 0; % preassigning a value;
percentVar = zeros(3,1); % pre-constructing a matrix for %VAR-original
for ii=1:3
    totalVariance = totalVariance + covariance(ii,ii); % summing all variances
end
for ii=1:3
    percentVar(ii,1) = covariance(ii,ii) / totalVariance; % calculating %VAR-original
end
percentVar = sort(percentVar,'descend'); % sorting values by descending value

totalVariancePC = sum(eigenVal);
percentVarPC = zeros(3,1); % pre-constructing a matrix for %VAR-PC
for ii=1:3
    percentVarPC(ii,1) = eigenVal(ii,1) / totalVariancePC; % calculating %VAR-PC
end
% percentVarPC = sort(percentVarPC,'descend'); % sorting values by descending value

figure;
% subplot(2,1,1);
% bar(percentVar');
% title('% VAR - original - Eigen Problem');

% subplot(2,1,2);
bar(percentVarPC');
title('% VAR - PC');

% reducing the dimensions of original data, constructing 2-d graph
featuredVects = zeros(3,2); % pre-constructing matrix vreating new basis
featuredVects(:,1) = eigenVect(:,1); % including eigenVect with largest variance for PC 1
featuredVects(:,2) = eigenVect(:,2); % including eigenVect with second largest variance for PC 2

finalData = featuredVects' * xCorrected; % projecting data matrix x onto new coordinate basis

figure;
plot(finalData(1,:),finalData(2,:),'k.'); % graph of data with new basis
title('Data with New Basis');
xlabel('x axis');
ylabel('y axis');

%% PCA via Single Value Decomposition C = U * Sigma * V'

% replotting reoriginal data to view PCA via SVD
figure; % plotting data for a visual interpretation
plot3(xCorrected(1,:),xCorrected(2,:),xCorrected(3,:),'k.')
grid on
title('Mean Corrected Data - SVD')
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
hold on;

[u,d,v] = svd(xCorrected); % using single value decomposition to find eigenvectors and eigenvalues

for ii=1:3
   u(:,ii) = u(:,ii) .* d(ii,ii); % multiplying eigenvectors by their eigenvalues for graph analysis below
end

d = diag(d); % extract diagonal elements of eigenvalue matrix
[d,index] = esort(d); % put magnitude of eigenvalues in descending order
u(:,:) = u(:,index); % sort eigenvectors likewise

% overlaying eigenvectors on plot of original data
plot3([0 u(1,1)],[0 u(2,1)],[0 u(3,1)],'r-');   % overlaying first eigenvector onto original data plot
plot3([0 u(1,2)],[0 u(2,2)],[0 u(3,2)],'g-');   % overlaying second eigenvector onto original data plot
plot3([0 u(1,3)],[0 u(2,3)],[0 u(3,3)],'b-');   % overlaying third eigenvector onto original data plot
hold off;

% reducing the dimensions of original data, constructing 2-d graph
featuredSVD_Vects = zeros(3,2); % pre-constructing matrix vreating new basis
featuredSVD_Vects(:,1) = u(:,1); % including eigenVect with largest variance for PC 1
featuredSVD_Vects(:,2) = u(:,2); % including eigenVect with second largest variance for PC 2

finalSVD_Data = featuredSVD_Vects' * xCorrected; % projecting data matrix x onto new coordinate basis

figure;
plot(finalSVD_Data(1,:),finalSVD_Data(2,:),'k.'); % graph of data with new basis
title('Data with New Basis - SVD');
xlabel('PC 1');
ylabel('PC 2');
