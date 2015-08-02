%% Nordita 2015, School Data Analysis
% 
% Lecture 2 
% Single s/c boundary methods
%
% Minimum Variance Analysis
% De Hofmann-Teller frame
% 

%% Example 1 ideal boundary, B change in one component
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T - T.start;                % define relative time in s from start
Bx  = 5*tanh((t-1800)/60);      % define Bx jump 
By  = -3*tanh((t-1800)/60);      % define Bx jump 
Bz  = t*0; 
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);                   % define scalar TSeries object
TS1.units = 'nT';
TS1.userData.LABLAXIS = 'B';

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  
irf_minvar_gui(TS1);

%% Example XX ideal boundary, B change in two components
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T - T.start;                % define relative time in s from start
Bx  = 5*tanh((t-1800)/60);      % define Bx jump 
By  = -3*tanh((t-1800)/10);      % define Bx jump 
Bz  = t*0; 
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);                   % define scalar TSeries object
TS1.units = 'nT';
TS1.userData.LABLAXIS = 'B';

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  
irf_minvar_gui(TS1);

%% Example XX ideal boundary, B change in two components + noise
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T - T.start;                % define relative time in s from start
Bx  = 5*tanh((t-1800)/60);      % define Bx jump 
By  = -3*tanh((t-1800)/10);      % define Bx jump 
Bz  = t*0; 
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);                   % define scalar TSeries object
TS1 = TS1 + 1*rand(TS1.length,3);
TS1.units = 'nT';
TS1.userData.LABLAXIS = 'B';

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  
irf_minvar_gui(TS1);

%% Example XX, read time interval of data
Tint = irf.tint('2001-03-01T01:00:00Z/2001-04-01T11:00:00Z');
varName = 'sc_lat__C4_JP_AUX_PSE';
varMat = local.c_read(varName,Tint,'mat');
varTs = local.c_read(varName,Tint,'ts');
irf_plot(varTs);

%% Example XX, magnetopause crossing Rosenqvist2008



%% Example  2 Plot multicomponent data
% Generate data with two components and plot in the same figure.
% Add legend text in lower left corner
% As you notice irfu-matlab interprets some common names for variables, 
% i.e. B2 is assumed to be magnetic field measurement by Cluster 2

y = exp(0.001*t).*cos(2*pi*t/180);	% z(t)=exp(0.001(t-to))*cos(t)
F = TSeries(T,[x y]);							          % B2 has two components, x & y
irf_plot(h,F)					          % plot in the same axis
irf_legend({'X','Y'},[0.02 0.02])	% add legend text with the same colors as lines

%% Example 3 Work with data, zoom in plots. 
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

tint_zoom = irf_time([2008 03 01 10 20 00]) + [0 1800]; % 10 min interval
irf_zoom(h,'x',tint_zoom);
irf_zoom(h,'y');

%% Example  4 Compare two data  
% Compare component-wise two datasets.
% Add one more label row with hours from the beginning of time interval

h=irf_plot({F,Fnew},'comp');
ylabel(h(1),'F_X');
title(h(1),T.start.utc('yyyy-mm-dd'));
ylabel(h(2),'F_Y');
irf_legend(h(1),{'FNew','Fnew=F*2+2'},[0.02 0.98])

Hours=TSeries(T,(T-T.start)/3600); % hours from the beginning of the time interval
irf_timeaxis(h(2),Hours,{'hours'})
irf_timeaxis(h(end),'nodate');

%% %% Example XX time interval
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T - T.start;                % define relative time in s from start
Bx  = 5*tanh((t-1800)/60);      % define Bx jump 
By  = -3*tanh((t-1800)/10);      % define Bx jump 
Bz  = t*0; 
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);                   % define scalar TSeries object
TS1.units = 'nT';
TS1.userData.LABLAXIS = 'B';

tInt = irf.tint('2002-03-04T09:40:00Z/2002-03-04T09:50:00Z');

TS2 = TS1.tlim(tInt);

h   = irf_plot(1,'newfigure');			        % initialize figure
irf_plot(h,TS1);						                % plot times series  
irf_minvar_gui(TS1);

