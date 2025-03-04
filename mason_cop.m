clc;
clear;
close all;

%% Initialize an empty matrix to hold COP path and variables data
LCOPx_rereg= []; RCOPx_rereg= []; LCOPy_rereg= []; RCOPy_rereg= [];
file_names = {};
Lcol_COP = []; Rcol_COP = [];
COPx_rereg= []; COPy_rereg= [];

LCOPx_vel_rereg = []; RCOPx_vel_rereg = []; LCOPy_vel_rereg = []; RCOPy_vel_rereg = [];
COPx_vel_rereg = []; COPy_vel_rereg = [];

PAPF_mx= []; Col_PAPF=[]; filename_PAPF=[];
d_copx_mx = []; d_copy_mx = [];
sstime_mx = []; dstime_mx = [];
v_copx_mx = []; v_copy_mx = [];
step_length_mx = []; step_width_mx = []; step_speed_mx=[]; step_time_mx=[];

%% read data
[file_list, path_n] = uigetfile('*.csv', 'Grab CSV Files', 'MultiSelect', 'on');

% Ensure file_list is a cell array (to handle single or multiple file selection)
if iscell(file_list) == 0
    file_list = {file_list};
end

% Loop through each selected CSV file
for u = 1:length(file_list)
    filename = file_list{u}; % Get the file name
    fullFilePath = fullfile(path_n, filename); % Combine path and file name

    % Read the CSV file as a cell array
    rawdata = readmatrix(fullFilePath);

% Remove the first 3 rows of rawdata
rawdata1 = rawdata(4:end, :);
% Extract the first column
col1 = rawdata1(:, 1);
% Find indices of the number 100 in column 1
indices100 = find(col1 == 100);

% Check the number of occurrences
if numel(indices100) == 1
    % Case 1: 100 appears one time
    % Force matrix: Remove rows from 1 row above the row where 100 is found to the end
    force = rawdata1(1:(indices100(1) - 2), :);
    
    % Motion matrix: Remove rows from 3 rows below the row where 100 is found to the start
    motion = rawdata1((indices100(1) + 4):end, :);
    
elseif numel(indices100) == 12
    % Case 2: 100 appears 12 times
    % Force matrix: Remove rows from 1 row above the 11th occurrence to the end
    force = rawdata1(1:(indices100(11) - 2), :);
    
    % Motion matrix: Remove rows from 3 rows below the 11th occurrence to the start
    motion = rawdata1((indices100(11) + 4):end, :);
else
    error('The number 100 was not found in the expected number of occurrences (1 or 12).');
end

%% load Markers
FR=100; % Motion Capture frame rate, Hz
frame=motion(:,1)-motion(1,1);
t=frame/FR; %time, sec
n=length(t);

% pelvis markers, meters
LASIS=0.001*motion(:,72:74);
RASIS=0.001*motion(:,75:77);
LPSIS=0.001*motion(:,78:80);
y_LPSIS=0.001*motion(:,79);
RPSIS=0.001*motion(:,81:83);
y_RPSIS=0.001*motion(:,82);
OPSIS=0.5*(RPSIS+LPSIS); %PSIS center(sacrum)
y_OPSIS=0.5*(y_LPSIS+y_RPSIS); %PSIS y center(y sacrum)
%CPSIS=mean([RPSIS;LPSIS],1);

% foot markers
LHE=0.001*motion(:,96:98);
x_LHE=0.001*motion(:,96);
y_LHE=0.001*motion(:,97);
z_LHE=0.001*motion(:,98);

LTO=0.001*motion(:,99:101);
x_LTO=0.001*motion(:,99);
y_LTO=0.001*motion(:,100);
z_LTO=0.001*motion(:,101);

RHE=0.001*motion(:,114:116);
x_RHE=0.001*motion(:,114);
y_RHE=0.001*motion(:,115);
z_RHE=0.001*motion(:,116);

RTO=0.001*motion(:,117:119);
x_RTO=0.001*motion(:,117);
y_RTO=0.001*motion(:,118);
z_RTO=0.001*motion(:,119);

z_lfocent= 0.5*(z_LHE+z_LTO); % z left foot centre
z_rfocent= 0.5*(z_RHE+z_RTO); % z right foot centre

%% Coordinate-Based Treadmill Algorithm_ EVENTS
if y_RHE(1,1)<0 || y_LHE(1,1)<0
  
% left heel-sacrum distance
Lheel=y_LHE-y_OPSIS;
% left toe-sacrum distance
Ltoe=-1*(y_LTO-y_OPSIS); % inverted

% right heel-sacrum distance
Rheel=y_RHE-y_OPSIS;
% right toe-sacrum distance
Rtoe=-1*(y_RTO-y_OPSIS); % inverted

else

% left heel-sacrum distance
Lheel=-(y_LHE-y_OPSIS);
% left toe-sacrum distance
Ltoe=(y_LTO-y_OPSIS); % inverted

% right heel-sacrum distance
Rheel=-(y_RHE-y_OPSIS);
% right toe-sacrum distance
Rtoe=(y_RTO-y_OPSIS); % inverted

end

%findpeaks/valleys left leg Events
[Lpks,flhs]=findpeaks(Lheel); %[peaks, Frames] left heel strike
%figure; findpeaks(Lheel);
% xlabel('frame');
% ylabel('left heel strike');
Lhstimes=(flhs-1)/FR; % left heel strike times

[Lvlys,flto]=findpeaks(Ltoe); %[valleys, Frames] left toe off
%figure; findpeaks(Ltoe);
% xlabel('frame');
% ylabel('left toe off');
Ltofftimes=(flto-1)/FR; % left toe off times

%findpeaks- right leg Events
[Rpks,frhs]=findpeaks(Rheel); %[peaks, Frames] right heel strike
%figure; findpeaks(Rheel);
% xlabel('frame');
% ylabel('right heel strike');
Rhstimes=(frhs-1)/FR; % right heel strike times

[Rvlys,frto]=findpeaks(Rtoe); %[valleys, Frames] right toe off
%figure; findpeaks(Rtoe);
% xlabel('frame');
% ylabel('right toe off');
Rtofftimes=(frto-1)/FR; % right toe off times

%% load force data
COP1x= force(:,9);
COP1y= force(:,10);
COP2x= force(:,18);
COP2y= force(:,19);
COP3x= force(:,27);
COP3y= force(:,28);
COP4x= force(:,36);
COP4y= force(:,37);

F1z= -force(:,5);
F2z= -force(:,14);
F3z= -force(:,23);
F4z= -force(:,32);

% Generate a Butterworth low-pass filter
Fs= 1000;
Fc= 10;
order= 4;
[b, a] = butter(order, Fc / (Fs / 2), 'low');
filtered_COP1x = filtfilt(b, a, COP1x); % Zero-phase filtering for no lag
filtered_COP1y = filtfilt(b, a, COP1y);
filtered_COP2x = filtfilt(b, a, COP2x);
filtered_COP2y = filtfilt(b, a, COP2y);
filtered_COP3x = filtfilt(b, a, COP3x);
filtered_COP3y = filtfilt(b, a, COP3y);
filtered_COP4x = filtfilt(b, a, COP4x);
filtered_COP4y = filtfilt(b, a, COP4y);

index=10:10:10*n;
Fdata=force(index,:);
F1_downsampled=Fdata(:,3:5); F1x_downsampled=F1_downsampled(:,1); F1y_downsampled=F1_downsampled(:,2); F1z_downsampled=F1_downsampled(:,3);
F2_downsampled=Fdata(:,12:14); F2x_downsampled=F2_downsampled(:,1); F2y_downsampled=F2_downsampled(:,2); F2z_downsampled=F2_downsampled(:,3);
F3_downsampled=Fdata(:,21:23); F3x_downsampled=F3_downsampled(:,1); F3y_downsampled=F3_downsampled(:,2); F3z_downsampled=F3_downsampled(:,3);
F4_downsampled=Fdata(:,30:32); F4x_downsampled=F4_downsampled(:,1); F4y_downsampled=F4_downsampled(:,2); F4z_downsampled=F4_downsampled(:,3);

% M1_downsampled=Fdata(:,6:8); M1x_downsampled=M1_downsampled(:,1); M1y_downsampled=M1_downsampled(:,2); M1z_downsampled=M1_downsampled(:,3);
% M2_downsampled=Fdata(:,15:17); M2x_downsampled=M2_downsampled(:,1); M2y_downsampled=M2_downsampled(:,2); M2z_downsampled=M2_downsampled(:,3);
% M3_downsampled=Fdata(:,24:26); M3x_downsampled=M3_downsampled(:,1); M3y_downsampled=M3_downsampled(:,2); M3z_downsampled=M3_downsampled(:,3);
% M4_downsampled=Fdata(:,33:35); M34_downsampled=M4_downsampled(:,1); M4y_downsampled=M4_downsampled(:,2); M4z_downsampled=M4_downsampled(:,3);
G=[0 0 -9.81];

%% % Step 1: Peak Anterior Propulsive Force (PAPF)
stance_F1 = find(F1z > 30); % Indices where V-GRF > threshold
stance_F2 = find(F2z > 30); % Indices where V-GRF > threshold
stance_F3 = find(F3z > 30); % Indices where V-GRF > threshold
stance_F4 = find(F4z > 30); % Indices where V-GRF > threshold

if y_RHE(1,1) < 0 || y_LHE(1,1) < 0
    F1y= -force(:,4);
    F2y= -force(:,13);
    F3y= -force(:,22);
    F4y= -force(:,31);
else
    F1y= force(:,4);
    F2y= force(:,13);
    F3y= force(:,22);
    F4y= force(:,31);
end

% Extract AP-GRF during the stance phase
AP_GRF1 = F1y(stance_F1);
AP_GRF2 = F2y(stance_F2);
AP_GRF3 = F3y(stance_F3);
AP_GRF4 = F4y(stance_F4);

propulsive_F1 = AP_GRF1 > 0; % Use only the positive values of AP-GRF (indicating propulsion)
propulsive_F2 = AP_GRF2 > 0;
propulsive_F3 = AP_GRF3 > 0;
propulsive_F4 = AP_GRF4 > 0;

%% Initialize re-registered matrices and PAPF
LCOPx_fullst1 = []; LCOPx_fullst2 = []; LCOPx_fullst3 = []; LCOPx_fullst4 = [];
LCOPy_fullst1 = []; LCOPy_fullst2 = []; LCOPy_fullst3 = []; LCOPy_fullst4 = [];
RCOPx_fullst1 = []; RCOPx_fullst2 = []; RCOPx_fullst3 = []; RCOPx_fullst4 = [];
RCOPy_fullst1 = []; RCOPy_fullst2 = []; RCOPy_fullst3 = []; RCOPy_fullst4 = [];
LPAPF1=[]; LPAPF2=[]; LPAPF3=[]; LPAPF4=[];
RPAPF1=[]; RPAPF2=[]; RPAPF3=[]; RPAPF4=[];
d_Lcopx1 = []; d_Lcopx2 = []; d_Lcopx3 = []; d_Lcopx4 = [];
d_Rcopx1 = []; d_Rcopx2 = []; d_Rcopx3 = []; d_Rcopx4 = [];
d_Lcopy1 = []; d_Lcopy2 = []; d_Lcopy3 = []; d_Lcopy4 = [];
d_Rcopy1 = []; d_Rcopy2 = []; d_Rcopy3 = []; d_Rcopy4 = [];
Lsstime1 = []; Lsstime2 = []; Lsstime3 = []; Lsstime4 = [];
Rsstime1 = []; Rsstime2 = []; Rsstime3 = []; Rsstime4 = [];
v_Lcopx1 = []; v_Lcopx2 = []; v_Lcopx3 = []; v_Lcopx4 = [];
v_Rcopx1 = []; v_Rcopx2 = []; v_Rcopx3 = []; v_Rcopx4 = [];
v_Lcopy1 = []; v_Lcopy2 = []; v_Lcopy3 = []; v_Lcopy4 = [];
v_Rcopy1 = []; v_Rcopy2 = []; v_Rcopy3 = []; v_Rcopy4 = [];
Lstep_length1 = []; Lstep_length2 = []; Lstep_length3 = []; Lstep_length4 = [];
Rstep_length1 = []; Rstep_length2 = []; Rstep_length3 = []; Rstep_length4 = [];
Lstep_time1 = []; Lstep_time2 = []; Lstep_time3 = []; Lstep_time4 = [];
Rstep_time1 = []; Rstep_time2 = []; Rstep_time3 = []; Rstep_time4 = [];
Lstep_speed1 = []; Lstep_speed2 = []; Lstep_speed3 = []; Lstep_speed4 = [];
Rstep_speed1 = []; Rstep_speed2 = []; Rstep_speed3 = []; Rstep_speed4 = [];
Lstep_width1 = []; Lstep_width2 = []; Lstep_width3 = []; Lstep_width4 = [];
Rstep_width1 = []; Rstep_width2 = []; Rstep_width3 = []; Rstep_width4 = [];
Ldstime1 = []; Ldstime2 = []; Ldstime3 = []; Ldstime4 = [];
Rdstime1 = []; Rdstime2 = []; Rdstime3 = []; Rdstime4 = [];

%% footcontact with forceplates
% alternative method: use round of stance_F instead of Fz_downsampled
% Force Plate 1
F1_footcontact = 10 * find(abs(F1z_downsampled) > 30);
if ~isempty(F1_footcontact)
    F1_heelstrike = F1_footcontact(1);
    F1_toeoff = F1_footcontact(end);
    frto1 = intersect(frto, F1_footcontact / 10);
    frhs1 = intersect(frhs, F1_footcontact / 10);
    flto1 = intersect(flto, F1_footcontact / 10);
    flhs1 = intersect(flhs, F1_footcontact / 10);
else
    F1_heelstrike = [];
    F1_toeoff = [];
    frto1 = [];
    frhs1 = [];
    flto1 = [];
    flhs1 = [];
end

% Force Plate 2
F2_footcontact = 10 * find(abs(F2z_downsampled) > 30);
if ~isempty(F2_footcontact)
    F2_heelstrike = F2_footcontact(1);
    F2_toeoff = F2_footcontact(end);
    frto2 = intersect(frto, F2_footcontact / 10);
    frhs2 = intersect(frhs, F2_footcontact / 10);
    flto2 = intersect(flto, F2_footcontact / 10);
    flhs2 = intersect(flhs, F2_footcontact / 10);
else
    F2_heelstrike = [];
    F2_toeoff = [];
    frto2 = [];
    frhs2 = [];
    flto2 = [];
    flhs2 = [];
end

% Force Plate 3
F3_footcontact = 10 * find(abs(F3z_downsampled) > 30);
if ~isempty(F3_footcontact)
    F3_heelstrike = F3_footcontact(1);
    F3_toeoff = F3_footcontact(end);
    frto3 = intersect(frto, F3_footcontact / 10);
    frhs3 = intersect(frhs, F3_footcontact / 10);
    flto3 = intersect(flto, F3_footcontact / 10);
    flhs3 = intersect(flhs, F3_footcontact / 10);
else
    F3_heelstrike = [];
    F3_toeoff = [];
    frto3 = [];
    frhs3 = [];
    flto3 = [];
    flhs3 = [];
end

% Force Plate 4
F4_footcontact = 10 * find(abs(F4z_downsampled) > 30);
if ~isempty(F4_footcontact)
    F4_heelstrike = F4_footcontact(1);
    F4_toeoff = F4_footcontact(end);
    frto4 = intersect(frto, F4_footcontact / 10);
    frhs4 = intersect(frhs, F4_footcontact / 10);
    flto4 = intersect(flto, F4_footcontact / 10);
    flhs4 = intersect(flhs, F4_footcontact / 10);
else
    F4_heelstrike = [];
    F4_toeoff = [];
    frto4 = [];
    frhs4 = [];
    flto4 = [];
    flhs4 = [];
end

%% COP path
if y_RHE(1,1) < 0 || y_LHE(1,1) < 0
% left
    if ~isempty(F1_heelstrike) && (y_LHE((F1_heelstrike)/10) > -0.025) && (y_LTO((F1_heelstrike)/10) < 0.565)
        LCOPx_fullst1 = filtered_COP1x(stance_F1); % COPx during stance time
        LCOPy_fullst1 = filtered_COP1y(stance_F1);
        LPAPF1 = max(AP_GRF1(propulsive_F1)); % Peak anterior propulsive force
        d_Lcopx1 = abs(filtered_COP1x(10*frto1) - filtered_COP1x(10*frhs1)); % Left ML-COP displacement
        d_Lcopy1 = abs(filtered_COP1y(10*frto1) - filtered_COP1y(10*frhs1)); % Left AP-COP displacement
        Lsstime1 = abs(frto1 - frhs1) / FR; % Left single stance time
        v_Lcopx1 = d_Lcopx1 / Lsstime1; % Left ML-COP velocity
        v_Lcopy1 = d_Lcopy1 / Lsstime1; % Left AP-COP velocity
        Lstep_length1 = abs(y_LHE(F1_heelstrike / 10) - y_RHE(frhs1));
        Lstep_time1= abs((F1_heelstrike / 10) - frhs1) / FR;
        Lstep_speed1= Lstep_length1/ Lstep_time1;
        Lstep_width1 = abs(x_LHE(F1_heelstrike / 10) - x_RHE(frhs1));
        Ldstime1 = abs((F1_heelstrike / 10) - frto1) / FR; % Left double stance time

    end

    if ~isempty(F2_heelstrike) && (y_LHE((F2_heelstrike)/10) > 0.580) && (y_LTO((F2_heelstrike)/10) < 1.170)
        LCOPx_fullst2 = filtered_COP2x(stance_F2); % COPx during stance time
        LCOPy_fullst2 = filtered_COP2y(stance_F2) - 605;
        LPAPF2 = max(AP_GRF2(propulsive_F2));
        d_Lcopx2 = abs(filtered_COP2x(10*frto2) - filtered_COP2x(10*frhs2));
        d_Lcopy2 = abs(filtered_COP2y(10*frto2) - filtered_COP2y(10*frhs2));
        Lsstime2 = abs(frto2 - frhs2) / FR;
        v_Lcopx2 = d_Lcopx2 / Lsstime2;
        v_Lcopy2 = d_Lcopy2 / Lsstime2;
        Lstep_length2 = abs(y_LHE(F2_heelstrike / 10) - y_RHE(frhs2));
        Lstep_time2 = abs((F2_heelstrike / 10) - frhs2) / FR;
        Lstep_speed2 = Lstep_length2 / Lstep_time2;
        Lstep_width2 = abs(x_LHE(F2_heelstrike / 10) - x_RHE(frhs2));
        Ldstime2 = abs((F2_heelstrike / 10) - frto2) / FR; % Left double stance time
    end

    if ~isempty(F3_heelstrike) && (y_LHE((F3_heelstrike)/10) > 1.18) && (y_LTO((F3_heelstrike)/10) < 1.775)
        LCOPx_fullst3 = filtered_COP3x(stance_F3); % COPx during stance time
        LCOPy_fullst3 = filtered_COP3y(stance_F3) - 1210; % COPx during stance time
        LPAPF3 = max(AP_GRF3(propulsive_F3));
        d_Lcopx3 = abs(filtered_COP3x(10*frto3) - filtered_COP3x(10*frhs3));
        d_Lcopy3 = abs(filtered_COP3y(10*frto3) - filtered_COP3y(10*frhs3));
        Lsstime3 = abs(frto3 - frhs3) / FR;
        v_Lcopx3 = d_Lcopx3 / Lsstime3;
        v_Lcopy3 = d_Lcopy3 / Lsstime3;
        Lstep_length3 = abs(y_LHE(F3_heelstrike / 10) - y_RHE(frhs3));
        Lstep_time3 = abs((F3_heelstrike / 10) - frhs3) / FR;
        Lstep_speed3 = Lstep_length3 / Lstep_time3;
        Lstep_width3 = abs(x_LHE(F3_heelstrike / 10) - x_RHE(frhs3));
        Ldstime3 = abs((F3_heelstrike / 10) - frto3) / FR; % Left double stance time
    end

    if ~isempty(F4_heelstrike) && (y_LHE((F4_heelstrike)/10) > 0.580) && (y_LTO((F4_heelstrike)/10) < 1.170)
        LCOPx_fullst4 = filtered_COP4x(stance_F4) - 405; % COPx during stance time
        LCOPy_fullst4 = filtered_COP4y(stance_F4) - 605; % COPx during stance time
        LPAPF4 = max(AP_GRF4(propulsive_F4));
        d_Lcopx4 = abs(filtered_COP4x(10*frto4) - filtered_COP4x(10*frhs4));
        d_Lcopy4 = abs(filtered_COP4y(10*frto4) - filtered_COP4y(10*frhs4));
        Lsstime4 = abs(frto4 - frhs4) / FR;
        v_Lcopx4 = d_Lcopx4 / Lsstime4;
        v_Lcopy4 = d_Lcopy4 / Lsstime4;
        Lstep_length4 = abs(y_LHE(F4_heelstrike / 10) - y_RHE(frhs4));
        Lstep_time4 = abs((F4_heelstrike / 10) - frhs4) / FR;
        Lstep_speed4 = Lstep_length4 / Lstep_time4;
        Lstep_width4 = abs(x_LHE(F4_heelstrike / 10) - x_RHE(frhs4));
        Ldstime4 = abs((F4_heelstrike / 10) - frto4) / FR; % Left double stance time
    end

    % right
    if ~isempty(F1_heelstrike) && (y_RHE((F1_heelstrike)/10) > -0.025) && (y_RTO((F1_heelstrike)/10) < 0.565)
        RCOPx_fullst1 = filtered_COP1x(stance_F1); % COPx during stance time
        RCOPy_fullst1 = filtered_COP1y(stance_F1);
        RPAPF1 = max(AP_GRF1(propulsive_F1)); % Peak anterior propulsive force
        d_Rcopx1 = abs(filtered_COP1x(10*flto1) - filtered_COP1x(10*flhs1)); % Right ML-COP displacement
        d_Rcopy1 = abs(filtered_COP1y(10*flto1) - filtered_COP1y(10*flhs1)); % Right AP-COP displacement
        Rsstime1 = abs(flto1 - flhs1) / FR; % Right single stance time
        v_Rcopx1 = d_Rcopx1 / Rsstime1; % Right ML-COP velocity
        v_Rcopy1 = d_Rcopy1 / Rsstime1; % Right AP-COP velocity
        Rstep_length1 = abs(y_RHE(F1_heelstrike / 10) - y_LHE(flhs1));
        Rstep_time1 = abs((F1_heelstrike / 10) - flhs1) / FR;
        Rstep_speed1 = Rstep_length1 / Rstep_time1;
        Rstep_width1 = abs(x_RHE(F1_heelstrike / 10) - x_LHE(flhs1));
        Rdstime1 = abs((F1_heelstrike / 10) - flto1) / FR; % Right double stance time
    end

    if ~isempty(F2_heelstrike) && (y_RHE((F2_heelstrike)/10) > 0.580) && (y_RTO((F2_heelstrike)/10) < 1.17)
        RCOPx_fullst2 = filtered_COP2x(stance_F2); % COPx during stance time
        RCOPy_fullst2 = filtered_COP2y(stance_F2) - 605;
        RPAPF2 = max(AP_GRF2(propulsive_F2));
        d_Rcopx2 = abs(filtered_COP2x(10*flto2) - filtered_COP2x(10*flhs2));
        d_Rcopy2 = abs(filtered_COP2y(10*flto2) - filtered_COP2y(10*flhs2));
        Rsstime2 = abs(flto2 - flhs2) / FR;
        v_Rcopx2 = d_Rcopx2 / Rsstime2;
        v_Rcopy2 = d_Rcopy2 / Rsstime2;
        Rstep_length2 = abs(y_RHE(F2_heelstrike / 10) - y_LHE(flhs2));
        Rstep_time2 = abs((F2_heelstrike / 10) - flhs2) / FR;
        Rstep_speed2 = Rstep_length2 / Rstep_time2;
        Rstep_width2 = abs(x_RHE(F2_heelstrike / 10) - x_LHE(flhs2));
        Rdstime2 = abs((F2_heelstrike / 10) - flto2) / FR; % Right double stance time
    end

    if ~isempty(F3_heelstrike) && (y_RHE((F3_heelstrike)/10) > 1.18) && (y_RTO((F3_heelstrike)/10) < 1.775)
        RCOPx_fullst3 = filtered_COP3x(stance_F3); % COPx during stance time
        RCOPy_fullst3 = filtered_COP3y(stance_F3) - 1210; % COPx during stance time
        RPAPF3 = max(AP_GRF3(propulsive_F3));
        d_Rcopx3 = abs(filtered_COP3x(10*flto3) - filtered_COP3x(10*flhs3));
        d_Rcopy3 = abs(filtered_COP3y(10*flto3)- filtered_COP3y(10*flhs3));
        Rsstime3 = abs(flto3 - flhs3) / FR;
        v_Rcopx3 = d_Rcopx3 / Rsstime3;
        v_Rcopy3 = d_Rcopy3 / Rsstime3;
        Rstep_length3 = abs(y_RHE(F3_heelstrike / 10) - y_LHE(flhs3));
        Rstep_time3 = abs((F3_heelstrike / 10) - flhs3) / FR;
        Rstep_speed3 = Rstep_length3 / Rstep_time3;
        Rstep_width3 = abs(x_RHE(F3_heelstrike / 10) - x_LHE(flhs3));
        Rdstime3 = abs((F3_heelstrike / 10) - flto3) / FR; % Right double stance time
    end

    if ~isempty(F4_heelstrike) && (y_RHE((F4_heelstrike)/10) > 0.580) && (y_RTO((F4_heelstrike)/10) < 1.17)
        RCOPx_fullst4 = filtered_COP4x(stance_F4) - 405; % COPx during stance time
        RCOPy_fullst4 = filtered_COP4y(stance_F4) - 605; % COPx during stance time
        RPAPF4 = max(AP_GRF4(propulsive_F4));
        d_Rcopx4 = abs(filtered_COP4x(10*flto4) - filtered_COP4x(10*flhs4));
        d_Rcopy4 = abs(filtered_COP4y(10*flto4) - filtered_COP4y(10*flhs4));
        Rsstime4 = abs(flto4 - flhs4) / FR;
        v_Rcopx4 = d_Rcopx4 / Rsstime4;
        v_Rcopy4 = d_Rcopy4 / Rsstime4;
        Rstep_length4 = abs(y_RHE(F4_heelstrike / 10) - y_LHE(flhs4));
        Rstep_time4 = abs((F4_heelstrike / 10) - flhs4) / FR;
        Rstep_speed4 = Rstep_length4 / Rstep_time4;
        Rstep_width4 = abs(x_RHE(F4_heelstrike / 10) - x_LHE(flhs4));
        Rdstime4 = abs((F4_heelstrike / 10) - flto4) / FR; % Right double stance time
    end

else

    if ~isempty(F3_heelstrike) && (y_LHE((F3_heelstrike)/10) < 1.815) && (y_LTO((F3_heelstrike)/10) > 1.24)
        LCOPx_fullst3 = 400 - filtered_COP3x(stance_F3);
        LCOPy_fullst3 = 1810 - filtered_COP3y(stance_F3);
        LPAPF3 = max(AP_GRF3(propulsive_F3)); % Peak anterior propulsive force
        d_Lcopx3 = abs(filtered_COP3x(10*frto3) - filtered_COP3x(10*frhs3)); % Left ML-COP displacement
        d_Lcopy3 = abs(filtered_COP3y(10*frto3) - filtered_COP3y(10*frhs3)); % Left AP-COP displacement
        Lsstime3 = abs(frto3 - frhs3) / FR; % Left single stance time
        v_Lcopx3 = d_Lcopx3 / Lsstime3; % Left ML-COP velocity
        v_Lcopy3 = d_Lcopy3 / Lsstime3; % Left AP-COP velocity
        Lstep_length3 = abs(y_LHE(F3_heelstrike / 10) - y_RHE(frhs3));
        Lstep_time3 = abs((F3_heelstrike / 10) - frhs3) / FR;
        Lstep_speed3 = Lstep_length3 / Lstep_time3;
        Lstep_width3 = abs(x_LHE(F3_heelstrike / 10) - x_RHE(frhs3));
        Ldstime3 = abs((F3_heelstrike / 10) - frto3) / FR; % Left double stance time
    end

    if ~isempty(F2_heelstrike) && (y_LHE((F2_heelstrike)/10) < 1.225) && (y_LTO((F2_heelstrike)/10) > 0.635)
        LCOPx_fullst2 = 400 - filtered_COP2x(stance_F2);
        LCOPy_fullst2 = 1205 - filtered_COP2y(stance_F2);
        LPAPF2 = max(AP_GRF2(propulsive_F2));
        d_Lcopx2 = abs(filtered_COP2x(10*frto2) - filtered_COP2x(10*frhs2));
        d_Lcopy2 = abs(filtered_COP2y(10*frto2) - filtered_COP2y(10*frhs2));
        Lsstime2 = abs(frto2 - frhs2) / FR;
        v_Lcopx2 = d_Lcopx2 / Lsstime2;
        v_Lcopy2 = d_Lcopy2 / Lsstime2;
        Lstep_length2 = abs(y_LHE(F2_heelstrike / 10) - y_RHE(frhs2));
        Lstep_time2 = abs((F2_heelstrike / 10) - frhs2) / FR;
        Lstep_speed2 = Lstep_length2 / Lstep_time2;
        Lstep_width2 = abs(x_LHE(F2_heelstrike / 10) - x_RHE(frhs2));
        Ldstime2 = abs((F2_heelstrike / 10) - frto2) / FR; % Left double stance time
    end

    if ~isempty(F1_heelstrike) && (y_LHE((F1_heelstrike)/10) < 0.625) && (y_LTO((F1_heelstrike)/10) > 0.035)
        LCOPx_fullst1 = 400 - filtered_COP1x(stance_F1);
        LCOPy_fullst1 = 600 - filtered_COP1y(stance_F1);
        LPAPF1 = max(AP_GRF1(propulsive_F1));
        d_Lcopx1 = abs(filtered_COP1x(10*frto1) - filtered_COP1x(10*frhs1));
        d_Lcopy1 = abs(filtered_COP1y(10*frto1) - filtered_COP1y(10*frhs1));
        Lsstime1 = abs(frto1 - frhs1) / FR;
        v_Lcopx1 = d_Lcopx1 / Lsstime1;
        v_Lcopy1 = d_Lcopy1 / Lsstime1;
        Lstep_length1 = abs(y_LHE(F1_heelstrike / 10) - y_RHE(frhs1));
        Lstep_time1 = abs((F1_heelstrike / 10) - frhs1) / FR;
        Lstep_speed1 = Lstep_length1 / Lstep_time1;
        Lstep_width1 = abs(x_LHE(F1_heelstrike / 10) - x_RHE(frhs1));
        Ldstime1 = abs((F1_heelstrike / 10) - frto1) / FR; % Left double stance time
    end

    if ~isempty(F4_heelstrike) && (y_LHE((F4_heelstrike)/10) < 1.225) && (y_LTO((F4_heelstrike)/10) > 0.635)
        LCOPx_fullst4 = 805 - filtered_COP4x(stance_F4); % COPx during stance time
        LCOPy_fullst4 = 1205 - filtered_COP4y(stance_F4); % COPx during stance time
        LPAPF4 = max(AP_GRF4(propulsive_F4));
        d_Lcopx4 = abs(filtered_COP4x(10*frto4) - filtered_COP4x(10*frhs4));
        d_Lcopy4 = abs(filtered_COP4y(10*frto4) - filtered_COP4y(10*frhs4));
        Lsstime4 = abs(frto4 - frhs4) / FR;
        v_Lcopx4 = d_Lcopx4 / Lsstime4;
        v_Lcopy4 = d_Lcopy4 / Lsstime4;
        Lstep_length4 = abs(y_LHE(F4_heelstrike / 10) - y_RHE(frhs4));
        Lstep_time4 = abs((F4_heelstrike / 10) - frhs4) / FR;
        Lstep_speed4 = Lstep_length4 / Lstep_time4;
        Lstep_width4 = abs(x_LHE(F4_heelstrike / 10) - x_RHE(frhs4));
        Ldstime4 = abs((F4_heelstrike / 10) - frto4) / FR; % Left double stance time
    end

    if ~isempty(F3_heelstrike) && (y_RHE((F3_heelstrike)/10) < 1.815) && (y_RTO((F3_heelstrike)/10) > 1.24)
        RCOPx_fullst3 = 400 - filtered_COP3x(stance_F3);
        RCOPy_fullst3 = 1810 - filtered_COP3y(stance_F3);
        RPAPF3 = max(AP_GRF3(propulsive_F3)); % Peak anterior propulsive force
        d_Rcopx3 = abs(filtered_COP3x(10*flto3) - filtered_COP3x(10*flhs3)); % Right ML-COP displacement
        d_Rcopy3 = abs(filtered_COP3y(10*flto3) - filtered_COP3y(10*flhs3)); % Right AP-COP displacement
        Rsstime3 = abs(flto3 - flhs3) / FR; % Right single stance time
        v_Rcopx3 = d_Rcopx3 / Rsstime3; % Right ML-COP velocity
        v_Rcopy3 = d_Rcopy3 / Rsstime3; % Right AP-COP velocity
        Rstep_length3 = abs(y_RHE(F3_heelstrike / 10) - y_LHE(flhs3));
        Rstep_time3 = abs((F3_heelstrike / 10) - flhs3) / FR;
        Rstep_speed3 = Rstep_length3 / Rstep_time3;
        Rstep_width3 = abs(x_RHE(F3_heelstrike / 10) - x_LHE(flhs3));
        Rdstime3 = abs((F3_heelstrike / 10) - flto3) / FR; % Right double stance time
    end

    if ~isempty(F2_heelstrike) && (y_RHE((F2_heelstrike)/10) < 1.225) && (y_RTO((F2_heelstrike)/10) > 0.635)
        RCOPx_fullst2 = 400 - filtered_COP2x(stance_F2);
        RCOPy_fullst2 = 1205 - filtered_COP2y(stance_F2);
        RPAPF2 = max(AP_GRF2(propulsive_F2));
        d_Rcopx2 = abs(filtered_COP2x(10*flto2) - filtered_COP2x(10*flhs2));
        d_Rcopy2 = abs(filtered_COP2y(10*flto2) - filtered_COP2y(10*flhs2));
        Rsstime2 = abs(flto2 - flhs2) / FR;
        v_Rcopx2 = d_Rcopx2 / Rsstime2;
        v_Rcopy2 = d_Rcopy2 / Rsstime2;
        Rstep_length2 = abs(y_RHE(F2_heelstrike / 10) - y_LHE(flhs2));
        Rstep_time2 = abs((F2_heelstrike / 10) - flhs2) / FR;
        Rstep_speed2 = Rstep_length2 / Rstep_time2;
        Rstep_width2 = abs(x_RHE(F2_heelstrike / 10) - x_LHE(flhs2));
        Rdstime2 = abs((F2_heelstrike / 10) - flto2) / FR; % Right double stance time
    end

    if ~isempty(F1_heelstrike) && (y_RHE((F1_heelstrike)/10) < 0.625) && (y_RTO((F1_heelstrike)/10) > 0.035)
        RCOPx_fullst1 = 400 - filtered_COP1x(stance_F1);
        RCOPy_fullst1 = 600 - filtered_COP1y(stance_F1);
        RPAPF1 = max(AP_GRF1(propulsive_F1));
        d_Rcopx1 = abs(filtered_COP1x(10*flto1) - filtered_COP1x(10*flhs1));
        d_Rcopy1 = abs(filtered_COP1y(10*flto1) - filtered_COP1y(10*flhs1));
        Rsstime1 = abs(flto1 - flhs1) / FR;
        v_Rcopx1 = d_Rcopx1 / Rsstime1;
        v_Rcopy1 = d_Rcopy1 / Rsstime1;
        Rstep_length1 = abs(y_RHE(F1_heelstrike / 10) - y_LHE(flhs1));
        Rstep_time1 = abs((F1_heelstrike / 10) - flhs1) / FR;
        Rstep_speed1 = Rstep_length1 / Rstep_time1;
        Rstep_width1 = abs(x_RHE(F1_heelstrike / 10) - x_LHE(flhs1));
        Rdstime1 = abs((F1_heelstrike / 10) - flto1) / FR; % Right double stance time
    end

    if ~isempty(F4_heelstrike) && (y_RHE((F4_heelstrike)/10) < 1.225) && (y_RTO((F4_heelstrike)/10) > 0.635)
        RCOPx_fullst4 = 805 - filtered_COP4x(stance_F4); % COPx during stance time
        RCOPy_fullst4 = 1205 - filtered_COP4y(stance_F4); % COPx during stance time
        RPAPF4 = max(AP_GRF4(propulsive_F4));
        d_Rcopx4 = abs(filtered_COP4x(10*flto4) - filtered_COP4x(10*flhs4));
        d_Rcopy4 = abs(filtered_COP4y(10*flto4) - filtered_COP4y(10*flhs4));
        Rsstime4 = abs(flto4 - flhs4) / FR;
        v_Rcopx4 = d_Rcopx4 / Rsstime4;
        v_Rcopy4 = d_Rcopy4 / Rsstime4;
        Rstep_length4 = abs(y_RHE(F4_heelstrike / 10) - y_LHE(flhs4));
        Rstep_time4 = abs((F4_heelstrike / 10) - flhs4) / FR;
        Rstep_speed4 = Rstep_length4 / Rstep_time4;
        Rstep_width4 = abs(x_RHE(F4_heelstrike / 10) - x_LHE(flhs4));
        Rdstime4 = abs((F4_heelstrike / 10) - flto4) / FR; % Right double stance time
    end
end    

%% reregister COP
nFrames = 101;
footlength = (1000 * sqrt(sum((LHE(1,:) - LTO(1,:)).^2)));
delta_t= 1/Fs;
LCOPx_vel_normalized1 = []; LCOPx_vel_normalized2 = []; LCOPx_vel_normalized3 = []; LCOPx_vel_normalized4 = [];
LCOPy_vel_normalized1 = []; LCOPy_vel_normalized2 = []; LCOPy_vel_normalized3 = []; LCOPy_vel_normalized4 = [];
RCOPx_vel_normalized1 = []; RCOPx_vel_normalized2 = []; RCOPx_vel_normalized3 = []; RCOPx_vel_normalized4 = [];
RCOPy_vel_normalized1 = []; RCOPy_vel_normalized2 = []; RCOPy_vel_normalized3 = []; RCOPy_vel_normalized4 = [];

% Initialize re-registered matrices
LCOPx_normalized1 = []; LCOPx_normalized2 = []; LCOPx_normalized3 = []; LCOPx_normalized4 = [];
LCOPy_normalized1 = []; LCOPy_normalized2 = []; LCOPy_normalized3 = []; LCOPy_normalized4 = [];
RCOPx_normalized1 = []; RCOPx_normalized2 = []; RCOPx_normalized3 = []; RCOPx_normalized4 = [];
RCOPy_normalized1 = []; RCOPy_normalized2 = []; RCOPy_normalized3 = []; RCOPy_normalized4 = [];

% 1
if exist('LCOPx_fullst1', 'var') && ~isempty(LCOPx_fullst1)
    % Step 1: Translate COP data
    LCOPx_translated1 = LCOPx_fullst1 - LCOPx_fullst1(1);
    LCOPy_translated1 = LCOPy_fullst1 - LCOPy_fullst1(1);
    
    % Combine translated COP data for PCA
    LCOP_translated1 = [LCOPx_translated1, LCOPy_translated1];
    
    % Step 2: Apply PCA to Align to Principal Axes
    [Lcoeff1, ~, ~] = pca(LCOP_translated1); % coeff1 contains PC1 and PC2 as columns
    LCOP_pca1 = LCOP_translated1 * Lcoeff1; % Align data to principal axes
    LCOP_nlf1 = (LCOP_pca1 / footlength) * 100;
    LCOPx_nlf1= -LCOP_nlf1(:, 2);
    LCOPy_nlf1= LCOP_nlf1(:, 1);

    % Output PCA-Aligned Data
    LCOPx_pca1 = -LCOP_pca1(:, 2); % New x-coordinates (aligned with PC2)
    LCOPy_pca1 = LCOP_pca1(:, 1); % New y-coordinates (aligned with PC1)
LCOPx_vel1 = zeros(size(LCOPx_pca1));
LCOPy_vel1 = zeros(size(LCOPy_pca1));
LCOPx_vel1(1) = (LCOPx_pca1(2) - LCOPx_pca1(1)) / delta_t;
LCOPy_vel1(1) = (LCOPy_pca1(2) - LCOPy_pca1(1)) / delta_t;

for i = 2:length(LCOPx_pca1)-1
    LCOPx_vel1(i) = (LCOPx_pca1(i+1) - LCOPx_pca1(i-1)) / (2 * delta_t);
    LCOPy_vel1(i) = (LCOPy_pca1(i+1) - LCOPy_pca1(i-1)) / (2 * delta_t);
end

LCOPx_vel1(end) = (LCOPx_pca1(end) - LCOPx_pca1(end-1)) / delta_t;
LCOPy_vel1(end) = (LCOPy_pca1(end) - LCOPy_pca1(end-1)) / delta_t;
    
    % Step 3: Normalize to 100% stance phase (101 frames)
    LframeIndices1 = linspace(1, length(LCOPx_pca1), nFrames);
    
    LCOPx_normalized1 = (interp1(1:length(LCOPx_nlf1), LCOPx_nlf1, LframeIndices1))';
    LCOPy_normalized1 = (interp1(1:length(LCOPy_nlf1), LCOPy_nlf1, LframeIndices1))';
    LCOPx_vel_normalized1 = (interp1(1:length(LCOPx_vel1), LCOPx_vel1, LframeIndices1))'; 
    LCOPy_vel_normalized1 = (interp1(1:length(LCOPy_vel1), LCOPy_vel1, LframeIndices1))';
end

if exist('RCOPx_fullst1', 'var') && ~isempty(RCOPx_fullst1)
    % Step 1: Translate COP data
    RCOPx_translated1 = RCOPx_fullst1 - RCOPx_fullst1(1);
    RCOPy_translated1 = RCOPy_fullst1 - RCOPy_fullst1(1);
    
    % Combine translated COP data for PCA
    RCOP_translated1 = [RCOPx_translated1, RCOPy_translated1];
    
    % Step 2: Apply PCA to Align to Principal Axes
    [Rcoeff1, ~, ~] = pca(RCOP_translated1); % coeff2 contains PC1 and PC2 as columns
    RCOP_pca1 = RCOP_translated1 * Rcoeff1; % Align data to principal axes
    RCOP_nlf1 = (RCOP_pca1 / footlength) * 100;
RCOPx_nlf1 = RCOP_nlf1(:, 2);
RCOPy_nlf1 = RCOP_nlf1(:, 1);

% Output PCA-Aligned Data
RCOPx_pca1 = RCOP_pca1(:, 2); % New x-coordinates (aligned with PC2)
RCOPy_pca1 = RCOP_pca1(:, 1); % New y-coordinates (aligned with PC1)
RCOPx_vel1 = zeros(size(RCOPx_pca1));
RCOPy_vel1 = zeros(size(RCOPy_pca1));
RCOPx_vel1(1) = (RCOPx_pca1(2) - RCOPx_pca1(1)) / delta_t;
RCOPy_vel1(1) = (RCOPy_pca1(2) - RCOPy_pca1(1)) / delta_t;

for i = 2:length(RCOPx_pca1)-1
    RCOPx_vel1(i) = (RCOPx_pca1(i+1) - RCOPx_pca1(i-1)) / (2 * delta_t);
    RCOPy_vel1(i) = (RCOPy_pca1(i+1) - RCOPy_pca1(i-1)) / (2 * delta_t);
end

RCOPx_vel1(end) = (RCOPx_pca1(end) - RCOPx_pca1(end-1)) / delta_t;
RCOPy_vel1(end) = (RCOPy_pca1(end) - RCOPy_pca1(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
RframeIndices1 = linspace(1, length(RCOPx_pca1), nFrames);

RCOPx_normalized1 = (interp1(1:length(RCOPx_nlf1), RCOPx_nlf1, RframeIndices1))';
RCOPy_normalized1 = (interp1(1:length(RCOPy_nlf1), RCOPy_nlf1, RframeIndices1))';
RCOPx_vel_normalized1 = (interp1(1:length(RCOPx_vel1), RCOPx_vel1, RframeIndices1))';
RCOPy_vel_normalized1 = (interp1(1:length(RCOPy_vel1), RCOPy_vel1, RframeIndices1))';
end

% 2
if exist('LCOPx_fullst2', 'var') && ~isempty(LCOPx_fullst2)
    % Step 1: Translate COP data
    LCOPx_translated2 = LCOPx_fullst2 - LCOPx_fullst2(1);
    LCOPy_translated2 = LCOPy_fullst2 - LCOPy_fullst2(1);

    % Combine translated COP data for PCA
    LCOP_translated2 = [LCOPx_translated2, LCOPy_translated2];

    % Step 2: Apply PCA to Align to Principal Axes
    [Lcoeff2, ~, ~] = pca(LCOP_translated2);
    LCOP_pca2 = LCOP_translated2 * Lcoeff2;
    LCOP_nlf2 = (LCOP_pca2 / footlength) * 100;

    LCOPx_nlf2= -LCOP_nlf2(:, 2);
LCOPy_nlf2= LCOP_nlf2(:, 1);

% Output PCA-Aligned Data
LCOPx_pca2 = -LCOP_pca2(:, 2); % New x-coordinates (aligned with PC2)
LCOPy_pca2 = LCOP_pca2(:, 1); % New y-coordinates (aligned with PC1)
LCOPx_vel2 = zeros(size(LCOPx_pca2));
LCOPy_vel2 = zeros(size(LCOPy_pca2));
LCOPx_vel2(1) = (LCOPx_pca2(2) - LCOPx_pca2(1)) / delta_t;
LCOPy_vel2(1) = (LCOPy_pca2(2) - LCOPy_pca2(1)) / delta_t;

for i = 2:length(LCOPx_pca2)-1
    LCOPx_vel2(i) = (LCOPx_pca2(i+1) - LCOPx_pca2(i-1)) / (2 * delta_t);
    LCOPy_vel2(i) = (LCOPy_pca2(i+1) - LCOPy_pca2(i-1)) / (2 * delta_t);
end

LCOPx_vel2(end) = (LCOPx_pca2(end) - LCOPx_pca2(end-1)) / delta_t;
LCOPy_vel2(end) = (LCOPy_pca2(end) - LCOPy_pca2(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
LframeIndices2 = linspace(1, length(LCOPx_pca2), nFrames);

LCOPx_normalized2 = (interp1(1:length(LCOPx_nlf2), LCOPx_nlf2, LframeIndices2))';
LCOPy_normalized2 = (interp1(1:length(LCOPy_nlf2), LCOPy_nlf2, LframeIndices2))';
LCOPx_vel_normalized2 = (interp1(1:length(LCOPx_vel2), LCOPx_vel2, LframeIndices2))'; 
LCOPy_vel_normalized2 = (interp1(1:length(LCOPy_vel2), LCOPy_vel2, LframeIndices2))';
end

% Process RCOP2
if exist('RCOPx_fullst2', 'var') && ~isempty(RCOPx_fullst2)
    % Step 1: Translate COP data
    RCOPx_translated2 = RCOPx_fullst2 - RCOPx_fullst2(1);
    RCOPy_translated2 = RCOPy_fullst2 - RCOPy_fullst2(1);

    % Combine translated COP data for PCA
    RCOP_translated2 = [RCOPx_translated2, RCOPy_translated2];

    % Step 2: Apply PCA to Align to Principal Axes
    [Rcoeff2, ~, ~] = pca(RCOP_translated2);
    RCOP_pca2 = RCOP_translated2 * Rcoeff2;
    RCOP_nlf2 = (RCOP_pca2 / footlength) * 100;

   RCOPx_nlf2= RCOP_nlf2(:, 2);
RCOPy_nlf2= RCOP_nlf2(:, 1);

% Output PCA-Aligned Data
RCOPx_pca2 = RCOP_pca2(:, 2); % New x-coordinates (aligned with PC2)
RCOPy_pca2 = RCOP_pca2(:, 1); % New y-coordinates (aligned with PC1)
RCOPx_vel2 = zeros(size(RCOPx_pca2));
RCOPy_vel2 = zeros(size(RCOPy_pca2));
RCOPx_vel2(1) = (RCOPx_pca2(2) - RCOPx_pca2(1)) / delta_t;
RCOPy_vel2(1) = (RCOPy_pca2(2) - RCOPy_pca2(1)) / delta_t;

for i = 2:length(RCOPx_pca2)-1
    RCOPx_vel2(i) = (RCOPx_pca2(i+1) - RCOPx_pca2(i-1)) / (2 * delta_t);
    RCOPy_vel2(i) = (RCOPy_pca2(i+1) - RCOPy_pca2(i-1)) / (2 * delta_t);
end

RCOPx_vel2(end) = (RCOPx_pca2(end) - RCOPx_pca2(end-1)) / delta_t;
RCOPy_vel2(end) = (RCOPy_pca2(end) - RCOPy_pca2(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
RframeIndices2 = linspace(1, length(RCOPx_pca2), nFrames);

RCOPx_normalized2 = (interp1(1:length(RCOPx_nlf2), RCOPx_nlf2, RframeIndices2))';
RCOPy_normalized2 = (interp1(1:length(RCOPy_nlf2), RCOPy_nlf2, RframeIndices2))';
RCOPx_vel_normalized2 = (interp1(1:length(RCOPx_vel2), RCOPx_vel2, RframeIndices2))'; 
RCOPy_vel_normalized2 = (interp1(1:length(RCOPy_vel2), RCOPy_vel2, RframeIndices2))';
end

    % 3
if exist('LCOPx_fullst3', 'var') && ~isempty(LCOPx_fullst3)
    % Step 1: Translate COP data
    LCOPx_translated3 = LCOPx_fullst3 - LCOPx_fullst3(1);
    LCOPy_translated3 = LCOPy_fullst3 - LCOPy_fullst3(1);

    % Combine translated COP data for PCA
    LCOP_translated3 = [LCOPx_translated3, LCOPy_translated3];

    % Step 2: Apply PCA to Align to Principal Axes
    [Lcoeff3, ~, ~] = pca(LCOP_translated3);
    LCOP_pca3 = LCOP_translated3 * Lcoeff3;
    LCOP_nlf3 = (LCOP_pca3 / footlength) * 100;

    LCOPx_nlf3= -LCOP_nlf3(:, 2);
LCOPy_nlf3= LCOP_nlf3(:, 1);

% Output PCA-Aligned Data
LCOPx_pca3 = -LCOP_pca3(:, 2); % New x-coordinates (aligned with PC2)
LCOPy_pca3 = LCOP_pca3(:, 1); % New y-coordinates (aligned with PC1)
LCOPx_vel3 = zeros(size(LCOPx_pca3));
LCOPy_vel3 = zeros(size(LCOPy_pca3));
LCOPx_vel3(1) = (LCOPx_pca3(2) - LCOPx_pca3(1)) / delta_t;
LCOPy_vel3(1) = (LCOPy_pca3(2) - LCOPy_pca3(1)) / delta_t;

for i = 2:length(LCOPx_pca3)-1
    LCOPx_vel3(i) = (LCOPx_pca3(i+1) - LCOPx_pca3(i-1)) / (2 * delta_t);
    LCOPy_vel3(i) = (LCOPy_pca3(i+1) - LCOPy_pca3(i-1)) / (2 * delta_t);
end

LCOPx_vel3(end) = (LCOPx_pca3(end) - LCOPx_pca3(end-1)) / delta_t;
LCOPy_vel3(end) = (LCOPy_pca3(end) - LCOPy_pca3(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
LframeIndices3 = linspace(1, length(LCOPx_pca3), nFrames);

LCOPx_normalized3 = (interp1(1:length(LCOPx_nlf3), LCOPx_nlf3, LframeIndices3))';
LCOPy_normalized3 = (interp1(1:length(LCOPy_nlf3), LCOPy_nlf3, LframeIndices3))';
LCOPx_vel_normalized3 = (interp1(1:length(LCOPx_vel3), LCOPx_vel3, LframeIndices3))'; 
LCOPy_vel_normalized3 = (interp1(1:length(LCOPy_vel3), LCOPy_vel3, LframeIndices3))';
end

% Process RCOP3
if exist('RCOPx_fullst3', 'var') && ~isempty(RCOPx_fullst3)
    % Step 1: Translate COP data
    RCOPx_translated3 = RCOPx_fullst3 - RCOPx_fullst3(1);
    RCOPy_translated3 = RCOPy_fullst3 - RCOPy_fullst3(1);

    % Combine translated COP data for PCA
    RCOP_translated3 = [RCOPx_translated3, RCOPy_translated3];

    % Step 2: Apply PCA to Align to Principal Axes
    [Rcoeff3, ~, ~] = pca(RCOP_translated3);
    RCOP_pca3 = RCOP_translated3 * Rcoeff3;
    RCOP_nlf3 = (RCOP_pca3 / footlength) * 100;

   RCOPx_nlf3= RCOP_nlf3(:, 2);
RCOPy_nlf3= RCOP_nlf3(:, 1);

% Output PCA-Aligned Data
RCOPx_pca3 = RCOP_pca3(:, 2); % New x-coordinates (aligned with PC2)
RCOPy_pca3 = RCOP_pca3(:, 1); % New y-coordinates (aligned with PC1)
RCOPx_vel3 = zeros(size(RCOPx_pca3));
RCOPy_vel3 = zeros(size(RCOPy_pca3));
RCOPx_vel3(1) = (RCOPx_pca3(2) - RCOPx_pca3(1)) / delta_t;
RCOPy_vel3(1) = (RCOPy_pca3(2) - RCOPy_pca3(1)) / delta_t;

for i = 2:length(RCOPx_pca3)-1
    RCOPx_vel3(i) = (RCOPx_pca3(i+1) - RCOPx_pca3(i-1)) / (2 * delta_t);
    RCOPy_vel3(i) = (RCOPy_pca3(i+1) - RCOPy_pca3(i-1)) / (2 * delta_t);
end

RCOPx_vel3(end) = (RCOPx_pca3(end) - RCOPx_pca3(end-1)) / delta_t;
RCOPy_vel3(end) = (RCOPy_pca3(end) - RCOPy_pca3(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
RframeIndices3 = linspace(1, length(RCOPx_pca3), nFrames);

RCOPx_normalized3 = (interp1(1:length(RCOPx_nlf3), RCOPx_nlf3, RframeIndices3))';
RCOPy_normalized3 = (interp1(1:length(RCOPy_nlf3), RCOPy_nlf3, RframeIndices3))';
RCOPx_vel_normalized3 = (interp1(1:length(RCOPx_vel3), RCOPx_vel3, RframeIndices3))'; 
RCOPy_vel_normalized3 = (interp1(1:length(RCOPy_vel3), RCOPy_vel3, RframeIndices3))';
end

% 4
if exist('LCOPx_fullst4', 'var') && ~isempty(LCOPx_fullst4)
    % Step 1: Translate COP data
    LCOPx_translated4 = LCOPx_fullst4 - LCOPx_fullst4(1);
    LCOPy_translated4 = LCOPy_fullst4 - LCOPy_fullst4(1);

    % Combine translated COP data for PCA
    LCOP_translated4 = [LCOPx_translated4, LCOPy_translated4];

    % Step 2: Apply PCA to Align to Principal Axes
    [Lcoeff4, ~, ~] = pca(LCOP_translated4);
    LCOP_pca4 = LCOP_translated4 * Lcoeff4;
    LCOP_nlf4 = (LCOP_pca4 / footlength) * 100;

    LCOPx_nlf4= -LCOP_nlf4(:, 2);
LCOPy_nlf4= LCOP_nlf4(:, 1);

% Output PCA-Aligned Data
LCOPx_pca4 = -LCOP_pca4(:, 2); % New x-coordinates (aligned with PC2)
LCOPy_pca4 = LCOP_pca4(:, 1); % New y-coordinates (aligned with PC1)
LCOPx_vel4 = zeros(size(LCOPx_pca4));
LCOPy_vel4 = zeros(size(LCOPy_pca4));
LCOPx_vel4(1) = (LCOPx_pca4(2) - LCOPx_pca4(1)) / delta_t;
LCOPy_vel4(1) = (LCOPy_pca4(2) - LCOPy_pca4(1)) / delta_t;

for i = 2:length(LCOPx_pca4)-1
    LCOPx_vel4(i) = (LCOPx_pca4(i+1) - LCOPx_pca4(i-1)) / (2 * delta_t);
    LCOPy_vel4(i) = (LCOPy_pca4(i+1) - LCOPy_pca4(i-1)) / (2 * delta_t);
end

LCOPx_vel4(end) = (LCOPx_pca4(end) - LCOPx_pca4(end-1)) / delta_t;
LCOPy_vel4(end) = (LCOPy_pca4(end) - LCOPy_pca4(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
LframeIndices4 = linspace(1, length(LCOPx_pca4), nFrames);

LCOPx_normalized4 = (interp1(1:length(LCOPx_nlf4), LCOPx_nlf4, LframeIndices4))';
LCOPy_normalized4 = (interp1(1:length(LCOPy_nlf4), LCOPy_nlf4, LframeIndices4))';
LCOPx_vel_normalized4 = (interp1(1:length(LCOPx_vel4), LCOPx_vel4, LframeIndices4))'; 
LCOPy_vel_normalized4 = (interp1(1:length(LCOPy_vel4), LCOPy_vel4, LframeIndices4))';
end

% Process RCOP4
if exist('RCOPx_fullst4', 'var') && ~isempty(RCOPx_fullst4)
    % Step 1: Translate COP data
    RCOPx_translated4 = RCOPx_fullst4 - RCOPx_fullst4(1);
    RCOPy_translated4 = RCOPy_fullst4 - RCOPy_fullst4(1);

    % Combine translated COP data for PCA
    RCOP_translated4 = [RCOPx_translated4, RCOPy_translated4];

    % Step 2: Apply PCA to Align to Principal Axes
    [Rcoeff4, ~, ~] = pca(RCOP_translated4);
    RCOP_pca4 = RCOP_translated4 * Rcoeff4;
    RCOP_nlf4 = (RCOP_pca4 / footlength) * 100;
RCOPx_nlf4= RCOP_nlf4(:, 2);
RCOPy_nlf4= RCOP_nlf4(:, 1);

% Output PCA-Aligned Data
RCOPx_pca4 = RCOP_pca4(:, 2); % New x-coordinates (aligned with PC2)
RCOPy_pca4 = RCOP_pca4(:, 1); % New y-coordinates (aligned with PC1)
RCOPx_vel4 = zeros(size(RCOPx_pca4));
RCOPy_vel4 = zeros(size(RCOPy_pca4));
RCOPx_vel4(1) = (RCOPx_pca4(2) - RCOPx_pca4(1)) / delta_t;
RCOPy_vel4(1) = (RCOPy_pca4(2) - RCOPy_pca4(1)) / delta_t;

for i = 2:length(RCOPx_pca4)-1
    RCOPx_vel4(i) = (RCOPx_pca4(i+1) - RCOPx_pca4(i-1)) / (2 * delta_t);
    RCOPy_vel4(i) = (RCOPy_pca4(i+1) - RCOPy_pca4(i-1)) / (2 * delta_t);
end

RCOPx_vel4(end) = (RCOPx_pca4(end) - RCOPx_pca4(end-1)) / delta_t;
RCOPy_vel4(end) = (RCOPy_pca4(end) - RCOPy_pca4(end-1)) / delta_t;

% Step 3: Normalize to 100% stance phase (101 frames)
RframeIndices4 = linspace(1, length(RCOPx_pca4), nFrames);

RCOPx_normalized4 = (interp1(1:length(RCOPx_nlf4), RCOPx_nlf4, RframeIndices4))';
RCOPy_normalized4 = (interp1(1:length(RCOPy_nlf4), RCOPy_nlf4, RframeIndices4))';
RCOPx_vel_normalized4 = (interp1(1:length(RCOPx_vel4), RCOPx_vel4, RframeIndices4))'; 
RCOPy_vel_normalized4 = (interp1(1:length(RCOPy_vel4), RCOPy_vel4, RframeIndices4))';
end

LCOPx_rereg= [LCOPx_rereg, LCOPx_normalized1, LCOPx_normalized2, LCOPx_normalized3, LCOPx_normalized4];
RCOPx_rereg= [RCOPx_rereg, RCOPx_normalized1, RCOPx_normalized2, RCOPx_normalized3, RCOPx_normalized4];
LCOPy_rereg= [LCOPy_rereg, LCOPy_normalized1, LCOPy_normalized2, LCOPy_normalized3, LCOPy_normalized4];
RCOPy_rereg= [RCOPy_rereg, RCOPy_normalized1, RCOPy_normalized2, RCOPy_normalized3, RCOPy_normalized4];

LCOPx_vel_rereg= [LCOPx_vel_rereg, LCOPx_vel_normalized1, LCOPx_vel_normalized2, LCOPx_vel_normalized3, LCOPx_vel_normalized4];
RCOPx_vel_rereg= [RCOPx_vel_rereg, RCOPx_vel_normalized1, RCOPx_vel_normalized2, RCOPx_vel_normalized3, RCOPx_vel_normalized4];
LCOPy_vel_rereg= [LCOPy_vel_rereg, LCOPy_vel_normalized1, LCOPy_vel_normalized2, LCOPy_vel_normalized3, LCOPy_vel_normalized4];
RCOPy_vel_rereg= [RCOPy_vel_rereg, RCOPy_vel_normalized1, RCOPy_vel_normalized2, RCOPy_vel_normalized3, RCOPy_vel_normalized4];

    LCOP_totalcol = size(LCOPx_normalized1, 2) + size(LCOPx_normalized2, 2) + size(LCOPx_normalized3, 2) + size(LCOPx_normalized4, 2);
    Lcol_COP = [Lcol_COP, LCOP_totalcol];
    RCOP_totalcol = size(RCOPx_normalized1, 2) + size(RCOPx_normalized2, 2) + size(RCOPx_normalized3, 2) + size(RCOPx_normalized4, 2);
    Rcol_COP = [Rcol_COP, RCOP_totalcol];

    PAPF_mx= [PAPF_mx, LPAPF1, LPAPF2, LPAPF3, LPAPF4,...
        RPAPF1, RPAPF2, RPAPF3, RPAPF4];

    PAPF_totalcol= size(LPAPF1, 2) + size(LPAPF2, 2) + size(LPAPF3, 2) + size(LPAPF4, 2) + ...
        size(RPAPF1, 2) + size(RPAPF2, 2) + size(RPAPF3, 2) + size(RPAPF4, 2);
    Col_PAPF= [Col_PAPF, PAPF_totalcol];

d_copx_mx = [d_copx_mx, d_Lcopx1, d_Lcopx2, d_Lcopx3, d_Lcopx4, d_Rcopx1, d_Rcopx2, d_Rcopx3, d_Rcopx4];
d_copy_mx = [d_copy_mx, d_Lcopy1, d_Lcopy2, d_Lcopy3, d_Lcopy4, d_Rcopy1, d_Rcopy2, d_Rcopy3, d_Rcopy4];
sstime_mx = [sstime_mx, Lsstime1, Lsstime2, Lsstime3, Lsstime4, Rsstime1, Rsstime2, Rsstime3, Rsstime4];
v_copx_mx = [v_copx_mx, v_Lcopx1, v_Lcopx2, v_Lcopx3, v_Lcopx4, v_Rcopx1, v_Rcopx2, v_Rcopx3, v_Rcopx4];
v_copy_mx = [v_copy_mx, v_Lcopy1, v_Lcopy2, v_Lcopy3, v_Lcopy4, v_Rcopy1, v_Rcopy2, v_Rcopy3, v_Rcopy4];

step_length_mx = [step_length_mx, Lstep_length1, Lstep_length2, Lstep_length3, Lstep_length4, ...
                  Rstep_length1, Rstep_length2, Rstep_length3, Rstep_length4];
step_time_mx = [step_time_mx, Lstep_time1, Lstep_time2, Lstep_time3, Lstep_time4, ...
                Rstep_time1, Rstep_time2, Rstep_time3, Rstep_time4];
step_speed_mx = [step_speed_mx, Lstep_speed1, Lstep_speed2, Lstep_speed3, Lstep_speed4, ...
                 Rstep_speed1, Rstep_speed2, Rstep_speed3, Rstep_speed4];
step_width_mx = [step_width_mx, Lstep_width1, Lstep_width2, Lstep_width3, Lstep_width4, ...
                 Rstep_width1, Rstep_width2, Rstep_width3, Rstep_width4];
dstime_mx = [dstime_mx, Ldstime1, Ldstime2, Ldstime3, Ldstime4, ...
             Rdstime1, Rdstime2, Rdstime3, Rdstime4];

   % Store file name in cell array
    [~, file_name, ~] = fileparts(filename);
    file_names{u} = file_name;
    
end

%% tables
% Add file names as column headers
Lheaders = [];
Rheaders = [];
for u = 1:length(file_names)
    Lheaders = [Lheaders, repmat({file_names{u}}, 1, Lcol_COP(u))];
    Rheaders = [Rheaders, repmat({file_names{u}}, 1, Rcol_COP(u))];
end

all_headers = [Lheaders, Rheaders];
COPx_rereg = [COPx_rereg, LCOPx_rereg, RCOPx_rereg];
COPx_rereg_wn = [all_headers; num2cell(COPx_rereg)]';
COPy_rereg= [COPy_rereg, LCOPy_rereg, RCOPy_rereg];
COPy_rereg_wn = [all_headers; num2cell(COPy_rereg)]';

COPx_vel_rereg = [COPx_vel_rereg, LCOPx_vel_rereg, RCOPx_vel_rereg];
COPx_vel_rereg_wn = [all_headers; num2cell(COPx_vel_rereg)]';
COPy_vel_rereg = [COPy_vel_rereg, LCOPy_vel_rereg, RCOPy_vel_rereg];
COPy_vel_rereg_wn = [all_headers; num2cell(COPy_vel_rereg)]';

filename_PAPF = repelem(file_names, Col_PAPF);
PAPF_mx = [filename_PAPF; num2cell(PAPF_mx)]';
d_copx_mx = [filename_PAPF; num2cell(d_copx_mx)]';
d_copy_mx = [filename_PAPF; num2cell(d_copy_mx)]';
sstime_mx = [filename_PAPF; num2cell(sstime_mx)]';
v_copx_mx = [filename_PAPF; num2cell(v_copx_mx)]';
v_copy_mx = [filename_PAPF; num2cell(v_copy_mx)]';
step_length_mx = [filename_PAPF; num2cell(step_length_mx)]';
step_time_mx = [filename_PAPF; num2cell(step_time_mx)]';
step_speed_mx = [filename_PAPF; num2cell(step_speed_mx)]';
step_width_mx = [filename_PAPF; num2cell(step_width_mx)]';
dstime_mx = [filename_PAPF; num2cell(dstime_mx)]';

%% Plot Results: COPx and COPy vs Frames
frames = 1:nFrames; % Frame indices

% COPx vs Frames
figure;
subplot(1,2,1);
hold on;
plot(frames, COPx_rereg, 'b', 'DisplayName', 'Force Plate 1'); % Force Plate 1
grid on;
title('COPx vs Frames');
xlabel('Frames');
ylabel('COPx (Lateral)');
legend show;

% COPy vs Frames
subplot(1,2,2);
hold on;
plot(frames, COPy_rereg, 'b', 'DisplayName', 'Force Plate 1'); % Force Plate 1
grid on;
title('COPy vs Frames');
xlabel('Frames');
ylabel('COPy (Anterior)');
legend show;


% velocity of COPx vs Frames
figure;
subplot(1,2,1);
hold on;
plot(frames, COPx_vel_rereg, 'b', 'DisplayName', 'Force Plate 1'); % Force Plate 1
grid on;
title('V_COPx vs Frames');
xlabel('Frames');
ylabel('V_COPx (Lateral)');
legend show;

% velocity of COPy vs Frames
subplot(1,2,2);
hold on;
plot(frames, COPx_vel_rereg, 'b', 'DisplayName', 'Force Plate 1'); % Force Plate 1
grid on;
title('V_COPy vs Frames');
xlabel('Frames');
ylabel('V_COPy (Anterior)');
legend show;

%% create mats and excels
headers = {'Participant', 'Task', 'd_copx', 'd_copy', 'sstime', 'v_copx', 'v_copy', 'CodeCondition', 'PAPF',...
    'step_length', 'step_width', 'step_time', 'step_speed', 'dstime'};% Extract the first column of d_copx as strings
first_column_words = string(d_copx_mx(:, 1));
% Initialize Participant column
Participant = regexp(first_column_words, '(?<=Over|DT)\d+', 'match', 'once'); % Extract numbers after Over or DT
Participant = str2double(Participant); % Convert to numeric (NaN for non-matching entries)
% Initialize CodeCondition
CodeCondition = zeros(size(first_column_words));
% Apply conditions
CodeCondition(contains(first_column_words, 'v', 'IgnoreCase', true)) = -0.5; % If the word contains 'v'
CodeCondition(contains(first_column_words, 'T', 'IgnoreCase', true)) = 0.5;  % If the word contains 't'
Participant= num2cell(Participant);
CodeCondition= num2cell(CodeCondition);
d_copx_mx(:,1)= cellstr(d_copx_mx(:,1));
% Combine headers and data, including Participant and CodeCondition
cop_variables_headers = [headers; Participant, d_copx_mx(:,1), d_copx_mx(:,2:end), d_copy_mx(:,2:end), sstime_mx(:,2:end), v_copx_mx(:,2:end), v_copy_mx(:,2:end),...
    CodeCondition, PAPF_mx(:,2:end), step_length_mx(:,2:end), step_width_mx(:,2:end),step_time_mx(:,2:end), step_speed_mx(:,2:end), dstime_mx(:,2:end)];

% Save matrices and their names
save('cop_variables.mat', 'cop_variables_headers');
save('copx_reregistered.mat', 'COPx_rereg_wn');
save('copy_reregistered.mat', 'COPy_rereg_wn');
save('copx_vel_reregistered.mat', 'COPx_vel_rereg_wn');
save('copy_vel_reregistered.mat', 'COPy_vel_rereg_wn');

writecell(COPx_rereg_wn, 'COP_reregistered.xlsx', 'Sheet', 'COPx');
writecell(COPy_rereg_wn, 'COP_reregistered.xlsx', 'Sheet', 'COPy');
writecell(COPx_vel_rereg_wn, 'COP_vel_reregistered.xlsx', 'Sheet', 'COPx');
writecell(COPy_vel_rereg_wn, 'COP_vel_reregistered.xlsx', 'Sheet', 'COPy');
writecell(cop_variables_headers, 'cop_variables.xlsx');
