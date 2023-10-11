function [phi,d] = eigsort(a)

%[phi,d] = eigsort(A)
%Returns eigenvector matrix and eigenvalue vector
%   in descending order of eigenvalues
%   A = input matrix
%   phi = eigenvector matrix
%   d = eigenvalue vector

[phi,d]=eig(a);         %solve eigenvector problem
d=diag(d);              %extract diagonal elements of eigenvalue matrix
[d,index]=esort(d);     %put magnitude of eigenvalues in descending order
phi(:,:)=phi(:,index);  %sort eigenvectors likewise
