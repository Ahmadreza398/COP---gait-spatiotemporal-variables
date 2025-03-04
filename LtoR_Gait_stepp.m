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

%% Markers
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

%% Spatiotemporal Gait Parameters
  %% step
if length(flhs) == length(frhs)
    % First condition
    if flhs(1) > frhs(1)

%     step_length = abs([y_RHE(frhs) - y_LHE(flhs); y_LHE(flhs(1:end-1)) - y_RHE(frhs(2:end))]);
    step_length = zeros(length(y_RHE(frhs) - y_LHE(flhs)) + length(y_LHE(flhs(1:end-1)) - y_RHE(frhs(2:end))), 1);
    step_length(1:2:end) =abs(y_RHE(frhs) - y_LHE(flhs));
    step_length(2:2:end) = abs(y_LHE(flhs(1:end-1)) - y_RHE(frhs(2:end)));

%     ti= abs([Rhstimes - Lhstimes; Lhstimes(1:end-1) - Rhstimes(2:end)]);
    ti = zeros(length(Rhstimes - Lhstimes)+ length(Lhstimes(1:end-1) - Rhstimes(2:end)), 1);
    ti(1:2:end)= abs(Rhstimes - Lhstimes);
    ti(2:2:end)= abs(Lhstimes(1:end-1) - Rhstimes(2:end));

    step_speed = step_length./ti;

%     step_width = abs([x_RHE(frhs) - x_LHE(flhs); x_LHE(flhs(1:end-1)) - x_RHE(frhs(2:end))]);
    step_width = zeros(length(x_RHE(frhs) - x_LHE(flhs))+length(x_LHE(flhs(1:end-1)) - x_RHE(frhs(2:end))), 1);
    step_width(1:2:end) = abs(x_RHE(frhs) - x_LHE(flhs));
    step_width(2:2:end) = abs(x_LHE(flhs(1:end-1)) - x_RHE(frhs(2:end)));
    else
%     step_length = abs([y_LHE(flhs) - y_RHE(frhs); y_RHE(frhs(1:end-1)) - y_LHE(flhs(2:end))]);
    step_length = zeros(length(y_LHE(flhs) - y_RHE(frhs)) + length(y_RHE(frhs(1:end-1)) - y_LHE(flhs(2:end))), 1);
    step_length(1:2:end) = abs(y_LHE(flhs) - y_RHE(frhs));
    step_length(2:2:end) = abs(y_RHE(frhs(1:end-1)) - y_LHE(flhs(2:end)));

%     ti= abs([Lhstimes - Rhstimes; Rhstimes(1:end-1) - Lhstimes(2:end)]);
    ti = zeros(length(Lhstimes - Rhstimes)+ length(Rhstimes(1:end-1) - Lhstimes(2:end)), 1);
    ti(1:2:end)= abs(Lhstimes - Rhstimes);
    ti(2:2:end)= abs(Rhstimes(1:end-1) - Lhstimes(2:end));

    step_speed = step_length./ti;

%     step_width = abs([x_LHE(flhs) - x_RHE(frhs); x_RHE(frhs(1:end-1)) - x_LHE(flhs(2:end))]);
    step_width = zeros(length(x_LHE(flhs) - x_RHE(frhs))+length(x_RHE(frhs(1:end-1)) - x_LHE(flhs(2:end))), 1);
    step_width(1:2:end) = abs(x_LHE(flhs) - x_RHE(frhs));
    step_width(2:2:end) = abs(x_RHE(frhs(1:end-1)) - x_LHE(flhs(2:end)));
    end
else
    % Second condition
    if flhs(1) > frhs(1)

%         step_length = abs([y_RHE(frhs(1:end-1)) - y_LHE(flhs); y_LHE(flhs) - y_RHE(frhs(2:end))]);
    step_length = zeros(length(y_RHE(frhs(1:end-1)) - y_LHE(flhs)) + length(y_LHE(flhs) - y_RHE(frhs(2:end))), 1);
    step_length(1:2:end) = abs(y_RHE(frhs(1:end-1)) - y_LHE(flhs));
    step_length(2:2:end) = abs(y_LHE(flhs) - y_RHE(frhs(2:end)));
%         ti= abs([Rhstimes(1:end-1) - Lhstimes; Lhstimes - Rhstimes(2:end)]);
    ti = zeros(length(Rhstimes(1:end-1) - Lhstimes)+ length(Lhstimes - Rhstimes(2:end)), 1);
    ti(1:2:end)= abs(Rhstimes(1:end-1) - Lhstimes);
    ti(2:2:end)= abs(Lhstimes - Rhstimes(2:end));

        step_speed = step_length./ti;

%         step_width = abs([x_RHE(frhs(1:end-1)) - x_LHE(flhs); x_LHE(flhs) - x_RHE(frhs(2:end))]);
    step_width = zeros(length(x_RHE(frhs(1:end-1)) - x_LHE(flhs))+length(x_LHE(flhs) - x_RHE(frhs(2:end))), 1);
    step_width(1:2:end) = abs(x_RHE(frhs(1:end-1)) - x_LHE(flhs));
    step_width(2:2:end) = abs(x_LHE(flhs) - x_RHE(frhs(2:end)));
    else
%         step_length = abs([y_LHE(flhs(1:end-1)) - y_RHE(frhs); y_RHE(frhs) - y_LHE(flhs(2:end))]);
    step_length = zeros(length(y_LHE(flhs(1:end-1)) - y_RHE(frhs)) + length(y_RHE(frhs) - y_LHE(flhs(2:end))), 1);
    step_length(1:2:end) = abs(y_LHE(flhs(1:end-1)) - y_RHE(frhs));
    step_length(2:2:end) = abs(y_RHE(frhs) - y_LHE(flhs(2:end)));

%         ti= abs([Lhstimes(1:end-1) - Rhstimes; Rhstimes - Lhstimes(2:end)]);
    ti = zeros(length(Lhstimes(1:end-1) - Rhstimes)+ length(Rhstimes - Lhstimes(2:end)), 1);
    ti(1:2:end)= abs(Lhstimes(1:end-1) - Rhstimes);
    ti(2:2:end)= abs(Rhstimes - Lhstimes(2:end));

        step_speed = step_length./ti;

%         step_width = abs([x_LHE(flhs(1:end-1)) - x_RHE(frhs); x_RHE(frhs) - x_LHE(flhs(2:end))]);
    step_width = zeros(length(x_LHE(flhs(1:end-1)) - x_RHE(frhs))+length(x_RHE(frhs) - x_LHE(flhs(2:end))), 1);
    step_width(1:2:end) = abs(x_LHE(flhs(1:end-1)) - x_RHE(frhs));
    step_width(2:2:end) = abs(x_RHE(frhs) - x_LHE(flhs(2:end)));

    end
end
