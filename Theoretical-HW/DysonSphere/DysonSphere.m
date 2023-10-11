prompt = 'Insert number of steps : ';
numSteps = input(prompt);
numSteps = numSteps - 1;
timeStep = 1;

%1 AU = 1.496e+11 m
r = 1.496e11; %radius of Dyson sphere
rSun = 696.3e6; %radius of sun/star
%calculating velocity of sphere to model Earth-like gravity at equator of
    %sphere
f = 9.81 + .0059; % force desired at equator outset to account for sun's gravity
m = 1;
a = f / 1;
w = sqrt(a / r);

%let theta be from 0 to 2pi
thetaStep = 2 * pi / numSteps;
theta = meshgrid(0:thetaStep:2 * pi);
%let phi be from 0 to pi
phiStep = pi / numSteps;
phi = meshgrid(0:phiStep:pi);

%To create a dyson sphere
[i,j,k] = sphere(numSteps);
i = i .* r;
j = j .* r;
k = k .* r;
surf(i,j,k)
alpha 0.01
hold on
%to create the sun
[x,y,z] = sphere(numSteps);
x = x .* rSun;
y = y .* rSun;
z = z .* rSun;
surf(x,y,z)

%gravitational force of sun for a person
G = 6.67408e-11;
mSun = 1.989e30; %mass of our sun
mTest = 1; %mass of test subject on biosphere
fiSun = -1 * ((G * mSun * mTest) / r^3) .* i; %force of gravity due to sun on biosphere
fjSun = -1 * ((G * mSun * mTest) / r^3) .* j; %force of gravity due to sun on biosphere
fkSun = -1 * ((G * mSun * mTest) / r^3) .* k; %force of gravity due to sun on biosphere

aCent = w^2 .* sqrt(i.^2 + j.^2);
%z dotted with aCent
iCent = aCent .* cos(theta);
fiCent = iCent .* mTest;
%y dotted with aCent
jCent = aCent .* sin(theta);
fjCent = jCent .* mTest;
%z dottted with aCent is 0 because since the rotation is along
    %the z axis, there is no z acceleration
kCent = w^2 .* k .* 0;
fkCent = kCent .* mTest;

%Calculating net aforce m_test experiences on biosphere
fiNet = fiSun - fiCent;
fjNet = fjSun - fjCent;
fkNet = fkSun - fkCent;

%quiver of net force from "gravity" an object feels at position x along the
    %curve of the sphere
quiver3(i,j,k,fiNet,fjNet,fkNet)

%quiver3(x,y,z,fxSun,fySun,fzSun)
%quiver3(i,j,k,fiCent,fjCent,fkCent)
hold off







