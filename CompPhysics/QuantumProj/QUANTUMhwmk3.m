clear;
clc;
a = 1;
x = (0:.01:1);

psi = zeros(5,length(x));
psi(1,:) = (8/(3*pi*sqrt(a)))*sin((pi/a).*x);
psi(2,:) = (1/(sqrt(a)))*sin(((2*pi)/a).*x);
psi(3,:) = (8/(5*pi*sqrt(a)))*sin(((3*pi)/a).*x);
psi(4,:) = 0;
psi(5,:) = (-8/(21*pi*sqrt(a)))*sin(((5*pi)/a).*x);

psiTotal = zeros(1,length(x));
PSIxt = zeros(1,length(x));
for n = 1:5
    psiTotal(1,:) = psiTotal(1,:) +  psi(n,:);
end
PSIxt(1,:) = 2/sqrt(a) * sin(((2*pi)/(a)).*x);
figure;
hold on
plot(x,psiTotal,'-b');
plot(x,PSIxt,'-r');
legend('PsiSubnTotal','Psi(x,t=0)');
