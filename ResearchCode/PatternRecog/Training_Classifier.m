function [ maskMatrixCat, w ] = Training_Classifier( X, mask_matrix, label )
%Training_Classifier description:
%   This function serves the purpose of using X which is a collection of
%   PET brain scans of patients with normal brain glucose metabolic
%   activity, and Alzheimer's Disease to train a classifier in order to
%   diagnose 

max = 128*128*91;
maskMatrixCat = ones(1,max) .* 2;
gg = 1;
for kk = 1:91
    for jj = 1:128
        for ii = 1:128
            maskMatrixCat(1,gg) = mask_matrix(ii,jj,kk);
%             % Testing for an incorrectly transformed cell
%             if maskMatrixCat(1,gg) == 2 
%                 disp('False conversion value at: (' + ii + ',' + jj + ',' + kk + ')' )
%             end
            gg = gg + 1;
        end
    end
end

maxDim = sum(maskMatrixCat);
xTransformed = ones(1,maxDim);
for ii = 1:42
    gg = 1;
    for jj = 1:max
        if maskMatrixCat(1,jj) == 1
            xTransformed(ii,gg) = X(ii,jj);
            gg = gg + 1;
        end
    end
end

xMean = mean(xTransformed);
xStdDev = std(xTransformed);
xTransformed = xTransformed';
for ii = 1:maxDim
    xTransformed(ii,:) = xTransformed(ii,:)' - xMean(1,ii);
    xTransformed(ii,:) = xTransformed(ii,:)' / xStdDev(1,ii);
end

[eigenVect,eigenVal] = svd(xTransformed,'econ'); % using single value decomposition to find eigenvectors and eigenvalues

for ii=1:3
   eigenVect(:,ii) = eigenVect(:,ii) .* eigenVal(ii,ii); % multiplying eigenvectors by their eigenvalues for graph analysis below
end

eigenVal = diag(eigenVal); % extract diagonal elements of eigenvalue matrix
[eigenVal,index] = esort(eigenVal); % put magnitude of eigenvalues in descending order
eigenVect(:,:) = eigenVect(:,index); % sort eigenvectors likewise

PC_Data = eigenVect' * xTransformed; % projecting data matrix x onto new coordinate basis

totalVariancePC = sum(eigenVal);
percentVarPC = zeros(42,1); % pre-constructing a matrix for %VAR-PC
for ii=1:42
    percentVarPC(ii,1) = eigenVal(ii,1) / totalVariancePC; % calculating %VAR-PC
end

figure;
subplot(2,1,1);
bar(percentVarPC');
title('% VAR - SVD');

figure;
gscatter(PC_Data(1,:),PC_Data(2,:),label);
title('PC 1 and PC 2');

cov1 = cov(PC_Data(1:22,:));
cov2 = cov(PC_Data(23:42,:));

S_W = cov1 + cov2;
InvS_W =  S_W ^ (-1);

PC_Data = PC_Data';
m1 = mean(PC_Data(1:22,:))';
m2 = mean(PC_Data(23:42,:))';
PC_Data = PC_Data';

S_B = (m1 - m2) * (m1 - m2)';

[w, Jw] = svd(InvS_W * S_B);

for ii = 1:4
    w(:,ii) = w(:,ii) .* Jw(ii,ii);
end

Jw = diag(Jw); % extract diagonal elements of eigenvalue matrix
[Jw,index] = esort(Jw); % put magnitude of eigenvalues in descending order
w(:,:) = w(:,index); % sort eigenvectors likewise
w = w(:,1);

yClass1 = w(:,1)' * PC_Data(:,1:22);  %projecting class 1 features onto w
yClass2 = w(:,1)' * PC_Data(:,23:42); %projecting class 2 features onto w

figure;
histogram(yClass1); %gaussian-like distribution of class 1
hold on
histogram(yClass2); %gaussian-like distribution of class 2
title('Blue - NL, Red - AL');
hold off
end

