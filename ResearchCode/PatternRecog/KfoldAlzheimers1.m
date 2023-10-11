load('DATA.mat')
% sizeX = size(X);
% gg = 1;
% jj = 1;
% set_1 = zeros(sizeX(1,1)*(2/3),sizeX(1,2)); label1 = cell(sizeX(1,1)*(2/3),1);
% set_2 = zeros(sizeX(1,1)*(1/3),sizeX(1,2)); label2= cell(sizeX(1,1)*(1/3),1);
% for ii = 3:3:sizeX(1,1)
%     set_1(gg,:) = X(ii-2,:); label1(gg,1) = y(ii-2,1);
%     set_1(gg+1,:) = X(ii-1,:); label1(gg+1,1) = y(ii-1,1);
%     set_2(jj,:) = X(ii,:); label2(jj,1) = y(ii,1);
%     gg = gg + 2;
%     jj = jj + 1;
% end
% size1 = size(set_1);
% size2 = size(set_2);
matri
for ii = 1:sizeX(1,1)/n:sizeX(1,1)
    for jj = 1:sizeX(1,1)
        zeros(sizeX(1,1)*(2/3),sizeX(1,2));
    end
        
end
    
   



%%
%[ maskMatrixCat, w ] = Training_Classifier( X, mask_matrix, y );
%Training_Classifier description:
%   This function serves the purpose of using X which is a collection of
%   PET brain scans of patients with normal brain glucose metabolic
%   activity, and Alzheimer's Disease to train a classifier in order to
%   diagnose 

maskMatrixCat = ones(1,sizeX(1,2)) .* 2;
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

maxDim = sum(maskMatrixCat(1,:));
set_1Transformed = ones(1,maxDim);
for ii = 1:size1(1,1)
    gg = 1;
    for jj = 1:sizeX(1,2)
        if maskMatrixCat(1,jj) == 1
            set_1Transformed(ii,gg) = set_1(ii,jj);
            gg = gg + 1;
        end
    end
end

xMean = mean(set_1Transformed);
xStdDev = std(set_1Transformed);
% xTransformed = xTransformed;
for ii = 1:maxDim
    set_1Transformed(:,ii) = set_1Transformed(:,ii) - xMean(1,ii);
    set_1Transformed(:,ii) = set_1Transformed(:,ii) / xStdDev(1,ii);
end

[eigenVect,eigenVal] = svd(set_1Transformed','econ'); % using single value decomposition to find eigenvectors and eigenvalues

eigenVal = diag(eigenVal); % extract diagonal elements of eigenvalue matrix
[eigenVal,index] = esort(eigenVal); % put magnitude of eigenvalues in descending order
eigenVect(:,:) = eigenVect(:,index); % sort eigenvectors likewise

for ii=1:size1(1,1)
   eigenVect(:,ii) = eigenVect(:,ii) .* eigenVal(ii,1);
end
eigenVect = eigenVect(:,1:27);
PC_Data = eigenVect' * set_1Transformed'; % projecting data matrix x onto new coordinate basis

PC_Data = PC_Data(1:2,:);
PCsize = size(PC_Data);

%%
cov1 = cov(PC_Data(:,1:15)');
cov2 = cov(PC_Data(:,16:28)');

S_W = cov1 + cov2;
InvS_W =  S_W ^ (-1);

PC_Data = PC_Data';
m1 = mean(PC_Data(1:15,:));
m2 = mean(PC_Data(16:28,:));
PC_Data = PC_Data';

S_B = (m1 - m2) * (m1 - m2)';

[w, Jw] = svd(InvS_W * S_B);

for ii = 1:PCsize(1,1)
    w(:,ii) = w(:,ii) .* Jw(ii,ii);
end

Jw = diag(Jw); % extract diagonal elements of eigenvalue matrix
[Jw,index] = esort(Jw); % put magnitude of eigenvalues in descending order
w(:,:) = w(:,index); % sort eigenvectors likewise

yClass1 = w' * PC_Data(:,1:15);  %projecting class 1 features onto w
yClass2 = w' * PC_Data(:,16:28); %projecting class 2 features onto w

figure;
histogram(yClass1(1,:)); %gaussian-like distribution of class 1
hold on
histogram(yClass2(1,:)); %gaussian-like distribution of class 2
title('Blue - NL, Red - AL');
hold off

figure;
plot(yClass1(1,:),yClass1(2,:),'ro');
hold on
plot(yClass2(1,:),yClass2(2,:),'bx');
xlabel('LDA 1');
ylabel('LDA 2');
hold off

%%
maxDim = sum(maskMatrixCat(1,:));
set_2Transformed = ones(1,maxDim);
for ii = 1:size2(1,1)
    gg = 1;
    for jj = 1:sizeX(1,2)
        if maskMatrixCat(1,jj) == 1
            set_2Transformed(ii,gg) = set_2(ii,jj);
            gg = gg + 1;
        end
    end
end

xMean = mean(set_2Transformed);
xStdDev = std(set_2Transformed);
% xTransformed = xTransformed;
for ii = 1:maxDim
    set_2Transformed(:,ii) = set_2Transformed(:,ii) - xMean(1,ii);
    set_2Transformed(:,ii) = set_2Transformed(:,ii) / xStdDev(1,ii);
end
% 
% [eigenVect,eigenVal] = svd(set_2Transformed','econ'); % using single value decomposition to find eigenvectors and eigenvalues
% 
% eigenVal = diag(eigenVal); % extract diagonal elements of eigenvalue matrix
% [eigenVal,index] = esort(eigenVal); % put magnitude of eigenvalues in descending order
% eigenVect(:,:) = eigenVect(:,index); % sort eigenvectors likewise
% 
% for ii=1:size2(1,1)
%    eigenVect(:,ii) = eigenVect(:,ii) .* eigenVal(ii,1);
% end
% eigenVect = eigenVect(:,1:12);
PC_Data = eigenVect' * set_2Transformed'; % projecting data matrix x onto new coordinate basis

PC_Data = PC_Data(1:2,:);

y1 = w' * PC_Data;

testlabel = cell(14,1);

for ii = 1:14
    if y1(1,ii) < 0
        testlabel{ii,1} = 'NL';
        string = ['Patient ', ii ,' does not have AL.'];
        disp(string)
    else 
        testlabel{ii,1} = 'AL';
        string = ['Patient ', ii ,' has AL.'];
        disp(string)
    end
end

TP = 0;
TN = 0;
FP = 0;
FN = 0;
for ii = 1:7
    if isequal(testlabel{ii,1},'AL')
            FP = FP + 1;
    end
    if isequal(testlabel{ii,1},'NL')
            TN = TN + 1;
    end
end
for ii = 8:14  
    if isequal(testlabel{ii,1},'NL')    
            FN = FN + 1;
    end
    if isequal(testlabel{ii,1},'AL')
            TP = TP + 1;
    end
end
conMat2 = [ TP FP; FN TN ];

%%
DV = [0;1];
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
imagesc(Img(:,:,57));