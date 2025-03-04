clc
clear
close all

%% load motion data
[file path]=uigetfile('*.xlsx');
data=xlsread([path file]);
frame=data(:,1)-data(1,1);
FR=100; % frame rate, Hz
t=frame/FR; %time, sec
n=length(t);

%% 3D Markers
% pelvis markers, meters
LASIS=0.001*data(:,72:74);
RASIS=0.001*data(:,75:77);
LPSIS=0.001*data(:,78:80);
y_LPSIS=0.001*data(:,79);
RPSIS=0.001*data(:,81:83);
y_RPSIS=0.001*data(:,82);
OPSIS=0.5*(RPSIS+LPSIS); %PSIS center(sacrum)
y_OPSIS=0.5*(y_LPSIS+y_RPSIS); %PSIS y center(y sacrum)
%CPSIS=mean([RPSIS;LPSIS],1);

% foot markers
LHE=0.001*data(:,96:98);
x_LHE=0.001*data(:,96);
y_LHE=0.001*data(:,97);
z_LHE=0.001*data(:,98);

LTO=0.001*data(:,99:101);
x_LTO=0.001*data(:,99);
y_LTO=0.001*data(:,100);
z_LTO=0.001*data(:,101);

RHE=0.001*data(:,114:116);
x_RHE=0.001*data(:,114);
y_RHE=0.001*data(:,115);
z_RHE=0.001*data(:,116);

RTO=0.001*data(:,117:119);
x_RTO=0.001*data(:,117);
y_RTO=0.001*data(:,118);
z_RTO=0.001*data(:,119);

z_lfocent= 0.5*(z_LHE+z_LTO); % z left foot centre
z_rfocent= 0.5*(z_RHE+z_RTO); % z right foot centre



%% load force data
Fdata=xlsread([path file],2);
index=10:10:10*n;
Fdata=Fdata(index,:);

LGRF=Fdata(:,3:5); LFx=LGRF(:,1); LFy=LGRF(:,2); LFz=LGRF(:,3);
RGRF=Fdata(:,12:14); RFx=RGRF(:,1); RFy=RGRF(:,2); RFz=RGRF(:,3);
LM=Fdata(:,6:8); LMx=LM(:,1); LMy=LM(:,2); LMz=LM(:,3);
RM=Fdata(:,15:17); RMx=LM(:,1); RMy=LM(:,2); RMz=LM(:,3);

% for p= 1:length(LFz)
% LCOPx(p)=-0.01*(LMy(p)/LFz(p)); LCOPy(p)=0.01*(LMx(p)/LFz(p));
% RCOPx(p)=-0.01*(RMy(p)/RFz(p)); RCOPy(p)=0.01*(RMx(p)/RFz(p));
% end


one_COP=0.1*Fdata(:,9:11); one_COPx=one_COP(:,1); one_COPy=one_COP(:,2);
two_COP=0.1*Fdata(:,18:20); two_COPx=two_COP(:,1); two_COPy=two_COP(:,2);
three_COP=0.1*Fdata(:,27:29); three_COPx=three_COP(:,1); three_COPy=three_COP(:,2);

G=[0 0 -9.81];

%% Coordinate-Based Treadmill Algorithm_ EVENTS
% left heel-sacrum distance
Lheel=y_LHE-y_OPSIS;
% left toe-sacrum distance
Ltoe=-1*(y_LTO-y_OPSIS); % inverted

% right heel-sacrum distance
Rheel=y_RHE-y_OPSIS;
% right toe-sacrum distance
Rtoe=-1*(y_RTO-y_OPSIS); % inverted

%findpeaks/valleys left leg Events
[Lpks,flhs]=findpeaks(Lheel); %[peaks, Frames] left heel strike
figure; findpeaks(Lheel);
xlabel('frame');
ylabel('left heel strike');
Lhstimes=(flhs-1)/FR; % left heel strike times

[Lvlys,flto]=findpeaks(Ltoe); %[valleys, Frames] left toe off
figure; findpeaks(Ltoe);
xlabel('frame');
ylabel('left toe off');
Ltofftimes=(flto-1)/FR; % left toe off times

%findpeaks- right leg Events
[Rpks,frhs]=findpeaks(Rheel); %[peaks, Frames] right heel strike
figure; findpeaks(Rheel);
xlabel('frame');
ylabel('right heel strike');
Rhstimes=(frhs-1)/FR; % right heel strike times

[Rvlys,frto]=findpeaks(Rtoe); %[valleys, Frames] right toe off
figure; findpeaks(Rtoe);
xlabel('frame');
ylabel('right toe off');
Rtofftimes=(frto-1)/FR; % right toe off times
%% foot progression angle
twoD_LHE = 0.001 * data(:, 96:97); % Left heel x, y coordinates
twoD_LTO = 0.001 * data(:, 99:100); % Left toe x, y coordinates
twoD_RHE = 0.001 * data(:, 114:115); % Right heel x, y coordinates
twoD_RTO = 0.001 * data(:, 117:118); % Right toe x, y coordinates

progression_direction = [0, 1]; % Progression direction (along y-axis)

% Assuming frto and flto are defined and contain indices
Lmidsfr = frto + 1; % Indices for left mid-stance moments
Rmidsfr = flto + 1; % Indices for right mid-stance moments

% Initialize arrays to store the results
LFPAs = zeros(1, length(Lmidsfr));
RFPAs = zeros(1, length(Rmidsfr));

% Calculate FPA for left foot
for b = 1:length(Lmidsfr)
    Lmid = Lmidsfr(b);
    
    % Extract heel and toe positions at mid-stance
    Lheel = twoD_LHE(Lmid, :);
    Ltoe = twoD_LTO(Lmid, :);
    
    % Calculate the foot axis vector
    Lfoot_axis = Ltoe - Lheel;
    
    % Normalize vectors
    Lfoot_axis_normalized = Lfoot_axis / norm(Lfoot_axis);
    progression_direction_normalized = progression_direction / norm(progression_direction);
    
    % Calculate the dot product of the normalized vectors
    Ldot_product = dot(Lfoot_axis_normalized, progression_direction_normalized);
    
    % Calculate the angle in radians
    Langle_radians = acos(Ldot_product);
    
    % Convert the angle to degrees
    Langle_degrees = rad2deg(Langle_radians);
    
    % Determine the sign of the angle using the cross product
    % Extend 2D vectors to 3D by adding a zero z-component
    Lfoot_axis_3D = [Lfoot_axis_normalized, 0];
    progression_direction_3D = [progression_direction_normalized, 0];
    
    % Calculate the cross product
    cross_product = cross(Lfoot_axis_3D, progression_direction_3D);
    
    % Check the z-component of the cross product to determine the sign
    if cross_product(3) > 0
        Langle_degrees = -Langle_degrees;
    end
    
    % Store the calculated FPA
    LFPAs(b) = Langle_degrees;
end

% Calculate FPA for right foot
for v = 1:length(Rmidsfr)
    Rmid = Rmidsfr(v);
    
    % Extract heel and toe positions at mid-stance
    Rheel = twoD_RHE(Rmid, :);
    Rtoe = twoD_RTO(Rmid, :);
    
    % Calculate the foot axis vector
    Rfoot_axis = Rtoe - Rheel;
    
    % Normalize vectors
    Rfoot_axis_normalized = Rfoot_axis / norm(Rfoot_axis);
    progression_direction_normalized = progression_direction / norm(progression_direction);
    
    % Calculate the dot product of the normalized vectors
    Rdot_product = dot(Rfoot_axis_normalized, progression_direction_normalized);
    
    % Calculate the angle in radians
    Rangle_radians = acos(Rdot_product);
    
    % Convert the angle to degrees
    Rangle_degrees = rad2deg(Rangle_radians);
    
    % Determine the sign of the angle using the cross product
    % Extend 2D vectors to 3D by adding a zero z-component
    Rfoot_axis_3D = [Rfoot_axis_normalized, 0];
    progression_direction_3D = [progression_direction_normalized, 0];
    
    % Calculate the cross product
    cross_product = cross(Rfoot_axis_3D, progression_direction_3D);
    
    % Check the z-component of the cross product to determine the sign
    if cross_product(3) < 0
        Rangle_degrees = -Rangle_degrees;
    end
    
    % Store the calculated FPA
    RFPAs(v) = Rangle_degrees;
end