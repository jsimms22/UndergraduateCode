function [ label ] = PET_Classifier( data, mask, w )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

maxDim = sum(mask);
xTransformed = ones(1,maxDim);
for ii = 1:42
    gg = 1;
    for jj = 1:maxDim
        if mask(1,jj) == 1
            xTransformed(ii,gg) = data(ii,jj);
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

y1 = w' * PC_Data(:,1:42);

label = cell(42,1);

for ii = 1:42
    if y1(1,ii) < 0
        label{ii,1} = 'NL';
        string = ['Patient ', ii ,' does not have AL.'];
        disp(string)
    else 
        label{ii,1} = 'AL';
        string = ['Patient ', ii ,' has AL.'];
        disp(string)
    end
end

end

