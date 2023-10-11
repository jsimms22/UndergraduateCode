clear % clears Workspace (removes existing variables/values stored
clc % clears Command Window/terminal

% Everything below here is for the purpose of collecting data from a NI-DAQ
% device. A session will be created. Sequentially, channels are added to a
% list that MATLAB will collect data from, data is then collected, and after 
% the time alloted concludes for the session a plot will display data from
% vector matrix - 'data.' 'data' is a (total channels by (by scan rate) * (total time)
% sized matrix. Each column corresponds to the channels, the 1st column will 
% be for the 1st channel added and so on. Anything followed by a % is not
% code, but a comment providing explanation for the line or preceding line.
% A %% indicates a section break, click 'Run' to run all sections. 'Run and
% Advance' will run a single section and prepare the next section to run.

v = daq.getVendors(); % lists avaliable vendor ids, remove ";" to view results in the Command Window
d = daq.getDevices(); % lists avaliable devices, remove ";" to view results in the Command Window
%% Assigns Session Characteristics

% session characteristics
minutes = 2; % how many minutes the test will run
rate = 100;  % scans per second captured
             % insertSessionName.inputSingleScan into the Command Window 
             % will return a single scan as a 7 by 1 matrix array 
%% Assigns Channel Characteristics

% channel characteristics
range = [-5,5]; % range of channel measurement range [low,high]
                % Supported ranges: -/+20, -/+10, -/+ 5, -/+1 Volts.
                % Depends on the channel which range is supported.
                
coupling = 'DC'; % MATLAB'S default is '{AC}', can also use 'DC'.

terminalconfig = 'SingleEnded'; % modes: 'SingleEnded', 'NonReferencedSingleEnded',
                                % 'Differential', and 'PseudoDifferential'.
                                % default is 'Differential'.
                                % delete comment syntax to assign specific
                                % channels with new configuration below.
                                
%% Creates Sessions and Modifies Session Using Assigned Characteristics

s = daq.createSession('ni'); % s = session name
s.DurationInSeconds = minutes*60;
s.Rate = rate;

%% Adds Channels to the Session and Modifies Channels Using Assigned Characteristics

% to add an analog input channel use sessionName.addAnalogInputChannel
% ('device_Name','channel_ID','measurement_Type');
% Grounded channels are not allowed to be added.

% input channel 2
    ch2 = s.addAnalogInputChannel('TorqueDev','ai0','Voltage');
    ch2.Range = range;
    ch2.Name = 'Strain';
    ch2.Coupling = coupling;
    ch2.TerminalConfig = 'SingleEnded'; % can use 'differential'
% input channel 5
    ch5 = s.addAnalogInputChannel('TorqueDev','ai1','Voltage');
    ch5.Range = range;
    ch5.Name = 'Current';
    ch5.Coupling = coupling;
    ch5.TerminalConfig = 'Differential';
% input channel 8
    ch8 = s.addAnalogInputChannel('TorqueDev','ai2','Voltage');
    ch8.Range = range;
    ch8.Name = 'Horizontal';
    ch8.Coupling = coupling;
    ch8.TerminalConfig = 'SingleEnded'; % can use 'SingleEnded'
% input channel 11
    ch11 = s.addAnalogInputChannel('TorqueDev','ai3','Voltage');
    ch11.Range = range;
    ch11.Name = 'Vertical';
    ch11.Coupling = coupling;
    ch11.TerminalConfig = 'Differential';

% input channel 1 Thermocouple
    ch1_t = s.addAnalogInputChannel('ThermalDev','ai0','Thermocouple');
    ch1_t.Name = 'Thermal';
    ch1_t.ThermocoupleType = 'J';
    
%% Begins Collection of Data, Outputing of Data, and Blocks Command Window From Use Until Session Completes

% Function startForeground(sessionName) starts operation of the session
% while blocking MATLAB command line and allowing other command line input
data = startForeground(s);
time = (1/rate:1/rate:s.DurationInSeconds);

%-----------------------------------------%
% Converting all voltage signals to desired units

% Strain linear relationship to voltage
for ii = 1:length(data(:,1))
    data(ii,1) = data(ii,1)*166.67 - .0018;
end
% Vertical position linear relationship to voltage
for ii = 1:length(data(:,3));
    data(ii,3) = data(ii,3)*44.6559 - 3.0887;
end
% Horizontal position linear relationship to voltage
for ii = 1:length(data(:,4))
    data(ii,4) = data(ii,4)*30.7305 + .5165;
end
%-----------------------------------------%

%-----------------------------------------%
% Strain data - original and transformed variant

figure;
plot(time,data(:,1));
xlabel('Time (s)');
ylabel('torque (N * m)');

Y = fft(data(:,1)); % fast fourier transform of strain data
L = length(data(:,1)); % length of data vector
P2 = abs(Y/L); % scaled form of the fourier transform
P1 = P2(1:L/2+1); % single sided amplitude of the transform
P1(2:end-1) = 2*P1(2:end-1); % removing noninteger values
Fs = (time(1,2) - time(1,1))^(-1); % frequency of original data
f = Fs*(0:(L/2))/L; % frequency for transform
figure;
plot(f,P1);
xlabel('Frequency (Hz)');
ylabel('|P1|');
%-----------------------------------------%

%-----------------------------------------%
% Current data

figure;
plot(time,data(:,2));
xlabel('Time (s)');
ylabel('Current (Amps)');
%-----------------------------------------%

%-----------------------------------------%
% Horizontal data

figure;
plot(time,data(:,3));
xlabel('Time (s)');
ylabel('Horizontal Position (V)');
%-----------------------------------------%

%-----------------------------------------%
% Vertical data

figure;
plot(time,data(:,4));
xlabel('Time (s)');
ylabel('Vertical Position (V)');
%-----------------------------------------%

%-----------------------------------------%
% Thermocouple data

figure;
plot(time,data(:,5));
xlabel('Time (s)');
ylabel('Temperature (Deg C)');
%-----------------------------------------%
