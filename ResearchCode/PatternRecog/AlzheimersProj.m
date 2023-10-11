load('DATA.mat')
%%
%[ maskMatrixCat, w ] = Training_Classifier( X, mask_matrix, y );
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
%%
xMean = mean(xTransformed);
xStdDev = std(xTransformed);
% xTransformed = xTransformed;
for ii = 1:maxDim
    xTransformed(:,ii) = xTransformed(:,ii) - xMean(1,ii);
    xTransformed(:,ii) = xTransformed(:,ii) / xStdDev(1,ii);
end

[eigenVect,eigenVal] = svd(xTransformed','econ'); % using single value decomposition to find eigenvectors and eigenvalues

eigenVal = diag(eigenVal); % extract diagonal elements of eigenvalue matrix
[eigenVal,index] = esort(eigenVal); % put magnitude of eigenvalues in descending order
eigenVect(:,:) = eigenVect(:,index); % sort eigenvectors likewise

for ii=1:42
   eigenVect(:,ii) = eigenVect(:,ii) .* eigenVal(ii,1);
end
eigenVect = eigenVect(:,1:41);
PC_Data = eigenVect' * xTransformed'; % projecting data matrix x onto new coordinate basis

totalVariancePC = sum(eigenVal);
percentVarPC = zeros(42,1); % pre-constructing a matrix for %VAR-PC
for ii=1:42
    percentVarPC(ii,1) = eigenVal(ii,1) / totalVariancePC; % calculating %VAR-PC
end

figure;
% subplot(2,1,1);
bar(percentVarPC');
title('% VAR - SVD');

figure;
gscatter(PC_Data(1,:),PC_Data(2,:),y);
title('PC 1 and PC 2');

PC_Data = PC_Data(1:2,:);
PCsize = size(PC_Data);

cov1 = cov(PC_Data(:,1:22)');
cov2 = cov(PC_Data(:,23:42)');

S_W = cov1 + cov2;
InvS_W =  S_W ^ (-1);

PC_Data = PC_Data';
m1 = mean(PC_Data(1:22,:));
m2 = mean(PC_Data(23:42,:));
PC_Data = PC_Data';

S_B = (m1 - m2) * (m1 - m2)';

[w, Jw] = svd(InvS_W * S_B);

for ii = 1:PCsize(1,1)
    w(:,ii) = w(:,ii) .* Jw(ii,ii);
end

Jw = diag(Jw); % extract diagonal elements of eigenvalue matrix
[Jw,index] = esort(Jw); % put magnitude of eigenvalues in descending order
w(:,:) = w(:,index); % sort eigenvectors likewise

yClass1 = w' * PC_Data(:,1:22);  %projecting class 1 features onto w
yClass2 = w' * PC_Data(:,23:42); %projecting class 2 features onto w

figure;
histogram(yClass1(1,:)); %gaussian-like distribution of class 1
hold on
histogram(yClass2(1,:)); %gaussian-like distribution of class 2
title('Blue - NL, Red - AL');
hold off

DV = [0; 1];
% DV = DV .* 8e5;
% DV2 = [0; -1];
% DV2 = DV2 .* 8e5;

figure;
plot(yClass1(1,:),yClass1(2,:),'bo');
hold on
plot(yClass2(1,:),yClass2(2,:),'rx');
plot([0 DV(1,1)],[0 DV(2,1)],'k-');
% plot([0 DV2(1,1)],[0 DV2(2,1)],'k-');
xlabel('LDA 1');
ylabel('LDA 2');
title('Blue - NL, Red - AL');
hold off

%%
%[ label ] = PET_Classifier(X, maskMatrixCat, w);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% maxDim = sum(maskMatrixCat);
% xTransformed = ones(1,maxDim);
% for ii = 1:42
%     gg = 1;
%     for jj = 1:max
%         if maskMatrixCat(1,jj) == 1
%             xTransformed(ii,gg) = X(ii,jj);
%             gg = gg + 1;
%         end
%     end
% end
% 
% xMean = mean(xTransformed);
% xStdDev = std(xTransformed);
% for ii = 1:maxDim
%     xTransformed(:,ii) = xTransformed(:,ii) - xMean(1,ii);
%     xTransformed(:,ii) = xTransformed(:,ii) / xStdDev(1,ii);
% end
% 
% % [eigenVect,eigenVal] = svd(xTransformed','econ'); % using single value decomposition to find eigenvectors and eigenvalues
% % 
% % for ii=1:3
% %    eigenVect(:,ii) = eigenVect(:,ii) .* eigenVal(ii,ii); % multiplying eigenvectors by their eigenvalues for graph analysis below
% % end
% % 
% % eigenVal = diag(eigenVal); % extract diagonal elements of eigenvalue matrix
% % [eigenVal,index] = esort(eigenVal); % put magnitude of eigenvalues in descending order
% eigenVect(:,:) = eigenVect(:,index); % sort eigenvectors likewise
% 
% PC_Data = eigenVect' * xTransformed; % projecting data matrix x onto new coordinate basis

y1 = w' * PC_Data;

label = cell(42,1);

for ii = 1:42
    if y1(1,ii) < 0
        label{ii,1} = 'AL';
        string = ['Patient ', ii ,' has AL.'];
        disp(string)
    else 
        label{ii,1} = 'NL';
        string = ['Patient ', ii ,' does not have AL.'];
        disp(string)
    end
end

%%
TP = 0;
TN = 0;
FP = 0;
FN = 0;
for ii = 1:22
    if isequal(label{ii,1},'AL')
            FP = FP + 1;
    end
    if isequal(label{ii,1},'NL')
            TN = TN + 1;
    end
end
for ii = 23:42    
    if isequal(label{ii,1},'NL')    
            FN = FN + 1;
    end
    if isequal(label{ii,1},'AL')
            TP = TP + 1;
    end
end
conMat1 = [ TP FP; FN TN ];

misclaRate = (FP + FN) / (TP + TN + FP + FN);
truPosFract = TP / (TP + FN);
falPosFract = FP / (FP + TN);

%%

EigImg = w * DV;
EigImg = eigenVect(:,1:2) * EigImg;
mask_vec = mask_matrix(:)';
Img_vec = repmat(nan,1,length(mask_vec));
gg = 1;
for ii = 1:max
    if mask_vec(1,ii) == 1
        Img_vec(1,ii) = EigImg(gg,1);
        gg = gg + 1;
    end
end

Img = reshape(Img_vec,[128,128,91]);

figure;
imagesc(Img(:,:,67));
figure;
imagesc(Img(:,:,31));
figure;
imagesc(Img(:,:,48));
figure;
imagesc(Img(:,:,78));

NL1 = reshape(X(1,:),[128,128,91]);
AD2 = reshape(X(42,:),[128,128,91]);

% figure;
% imagesc(NL1(:,:,67));
figure;
imagesc(NL1(:,:,31));
% figure;
% imagesc(NL1(:,:,48));
% figure;
% imagesc(NL1(:,:,78));
% 
% figure;
% imagesc(AD2(:,:,67));
figure;
imagesc(AD2(:,:,31));
% figure;
% imagesc(AD2(:,:,48));
% figure;
% imagesc(AD2(:,:,78));






