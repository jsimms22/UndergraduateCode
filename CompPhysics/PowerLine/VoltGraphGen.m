z = importdata('voltmap.dat');
imax = 51;
jmax = 51;
z = z .* (1/6e4);

z2 = importdata('voltmap2.dat');

dx = 1 /(imax - 1);
dy = 1 /(jmax - 1);

x = meshgrid(0:199/52:199);
y = meshgrid(0:199/52:199);
x = x';

figure;
surf(x,y,z);
figure;
surf(x,y,z2);
figure;
plot(x(:,1),z(:,27));
xlabel('Distance From Wire, dx');
ylabel('Percent Voltage');
%title('Orthongonal View of a Wire');

zed = zeros(53,53);
for ii = 1:53
    zed(:,ii) = z(:,27);
end
%y = meshgrid(0:dx:dx*52);
figure;
surf(x,y,zed);
xlabel('Distance From Wire, dx');
ylabel('Distance Along the Wire, dy');
zlabel('Percent Voltage');
%title('Electricity Along A Wire');
% legend('dx not equal to dy');


