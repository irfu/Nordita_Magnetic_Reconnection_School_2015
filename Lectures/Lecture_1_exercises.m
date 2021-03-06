%% Nordita 2015, School Data Analysis
% 
% Lecture 1
%
%  Data representation
%  Data reading
%  Data plotting
%  Data access from databases
%  Data massage
%

%% Ex 1: Download and examine CDF files
mkdir Nordita2015
cd Nordita2015/
mkdir Event_20020330_1311
cd Event_20020330_1311

Tint = irf.tint('2002-03-30T13:11:30Z/2002-03-30T13:12:00Z'); % define event time interval
caa_download(Tint,'C?_CP_FGM_FULL');   % download FGM data from CSA
info = spdfcdfinfo('CAA/C1_CP_FGM_FULL/C1_CP_FGM_FULL__20020330_131130_20020330_131200_V140306.cdf'); % query metadata

%% Ex 2: Time
% Simple operations with time

utcT1 = '2002-03-04T09:30:00Z'; % UTC string
EpochTT1 = EpochTT(utcT1);      % New EpochTT object
EpochTT2 = EpochTT1 + 10;       % New EpochTT object offset by 10 sec
offset = EpochTT2 - EpochTT1;   % Offset between two times in sec
if EpochTT2>EpochTT1            % Compare times
  disp('larger!')
end 

EpochTT0 = EpochTT1 + (-5);     % New EpochTT with negative offset of 5 sec
if EpochTT0<EpochTT1            % Compare times
  disp('smaller!')
end 

EpochUnix1 = EpochUnix(utcT1);  % New EpochUnix
if EpochUnix1 == EpochTT1       % Compare times
  disp('equal!')
end
epochUnix = EpochUnix1.epochUnix; % double value of Unix epoch [sec]
ttns = EpochUnix1.ttns;         % int64 value of TT epoch [ns]
if ttns == EpochTT1.ttns        % Compare times
  disp('equal!')
end
disp(EpochUnix1.utc)            % convert to UTC string

TTarray = EpochTT0:1:EpochTT2;  % New time array, with 1 sec step
TintLim = ...                   % Time interval
  irf.tint('2002-03-04T09:30:00Z/2002-03-04T09:30:05Z');
[idxIn,TTarrayIn] = ...
  TTarray.tlim(TintLim);        % Limit TTarray by TintLim
[idxOut,TTarrayOut] = ...
  TTarray.tlim(TintLim,1);      % Limit TTarray by TintLim, XOR mode
flagBefore = TTarray < EpochTT1;% compare times

%% Ex 3.1: Artifical times series (TS)
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% showing exponentially growing wave and plot. It is good idead to get used 
% to using axis handles (variable 'h' in example). 

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T - T.start;                     % define relative time in s from start
x   = exp(0.001*t).*sin(2*pi*t/180);        % define function x(t)=exp(0.001(t-to))*sin(t-to)
TS1 = irf.ts_scalar(T,x);                   % define scalar TSeries object

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  

%% Ex 3.2: Multicomponent time series
% Generate data with two components and plot in the same figure.
% Add legend text in lower left corner
% As you notice irfu-matlab interprets some common names for variables, 
% i.e. B2 is assumed to be magnetic field measurement by Cluster 2

y = exp(0.001*t).*cos(2*pi*t/180);	% z(t)=exp(0.001(t-to))*cos(t)
F = irf.ts_vec_xy(T,[x y]);				  % B2 has two components, x & y
irf_plot(h,F)					              % plot in the same axis
irf_legend({'X','Y'},[0.02 0.02])	  % add legend text with the same colors as lines

%% Ex 4.1: Plot TSs, zoom in plots. 
% - Generate second data set that is a function of first.
%     Fnew = F*1.2 + 2	 
% - Plot both data sets in separate panels.
% - Zoom in to smaller 30min time interval .

Fnew = F*2 + 2;

h  = irf_plot(2,'newfigure');		% empty figure with 2 panels
irf_plot(h(1),F);
ylabel(h(1),'F');
irf_legend(h(1),{'X','Y'},[0.02 0.98],'fontsize',20)
irf_plot(h(2),Fnew);
ylabel(h(2),'F_{new} = F * 2 + 2');

Tint = irf.tint('2002-03-04T09:40:00Z/2002-03-04T10:20:00Z'); % 30 min interval
irf_zoom(h,'x',Tint)
irf_zoom(h,'y');

%% Ex 4.2: Plotting: Compare two TSs
% Compare component-wise two datasets.
% Add one more label row with hours from the beginning of time interval

h = irf_plot({F,Fnew},'comp');
ylabel(h(1),'F_X');
title(h(1),T.start.utc('yyyy-mm-dd'));
ylabel(h(2),'F_Y');
irf_legend(h(1),{'FNew','Fnew=F*2+2'},[0.02 0.98])

Hours = TSeries(T,(T-T.start)/3600); % hours from the beginning of the time interval
irf_timeaxis(h(2),Hours,{'hours'})
irf_timeaxis(h(end),'nodate');

%% Ex 5: dataobj and time series
D=dataobj('CAA/C1_CP_FGM_FULL/C1_CP_FGM_FULL__20020330_131130_20020330_131200_V140306.cdf');
display(D)                      % display varibles in dataobj
plot(D,'B_mag__C1_CP_FGM_FULL') % plot one of the varibles in dataobj

b1 = ...                        % extract one variable as simple array
  getmat(D,'B_vec_xyz_gse__C1_CP_FGM_FULL');
disp(size(b1))
irf_plot(b1)

B1 = ...                        % extract one variable as TS
  get_ts(D,'B_vec_xyz_gse__C1_CP_FGM_FULL');
display(B1)
irf_plot(B1)