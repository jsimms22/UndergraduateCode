load('DATA.mat')
sizeX = size(X);
n = 5;
rows = sizeX(1,1) / n;
rows = floor(rows);
A = cell(n+1,1);
gg = 1;
jj = 1;
for ii = rows:rows:sizeX(1,1)
    A{gg} = X(jj:ii,:);
    gg = gg + 1;
    jj = ii + 1;
end

if ne(rows,sizeX(1,1) / n)
    correction = rows * n + 1;
    A{n+1} = X(correction:sizeX(1,1),:);
end

y = kfoldfunctiontest(A);



