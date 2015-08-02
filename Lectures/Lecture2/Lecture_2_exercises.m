%% Nordita 2015, School Data Analysis
% 
% Lecture 2 
% Single s/c boundary methods
%
% Minimum Variance Analysis
% De Hofmann-Teller frame
% 
% As real data we use Paschmann 2005 paper


%% Example 1 ideal boundary, B change in one component
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');% define time line as EpochTT object
t   = T - T.start;                    % define relative time in s from start
Bx  =  5*tanh((t-1800)/60);           % define Bx +-5nT jump and 1min width
By  = -3*tanh((t-1800)/60);           % define Bx -+3nT jump and 1min width
Bz  = t*0;                            % define Bz to be zero
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);   % define scalar TSeries object
TS1.units = 'nT';                     % units
TS1.userData.LABLAXIS = 'B';          % plot label axis

h = irf_plot(1,'newfigure');		  	  % initialize figure with one panel
irf_plot(h,TS1);						          % plot times series  
irf_minvar_gui(TS1);                  % run minimum variance analysis on time series

%% Example XX ideal boundary, B change in two components
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');% define time line as EpochTT object
t   = T - T.start;                    % define relative time in s from start
Bx  =  5*tanh((t-1800)/60);           % define Bx +-5nT jump and 1min width
By  = -3*tanh((t-1800)/10);           % define Bx -+3nT jump and 10s width
Bz  = t*0;                            % define Bz to be zero
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);   % define scalar TSeries object
TS1.units = 'nT';                     % units
TS1.userData.LABLAXIS = 'B';          % plot label axis

h = irf_plot(1,'newfigure');		  	  % initialize figure with one panel
irf_plot(h,TS1);						          % plot times series  
irf_minvar_gui(TS1);                  % run minimum variance analysis on time series

%% Example XX ideal boundary, B change in two components + noise
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');% define time line as EpochTT object
t   = T - T.start;                    % define relative time in s from start
Bx  =  5*tanh((t-1800)/60);           % define Bx jump +- 5nT and width 1 min
By  = -3*tanh((t-1800)/10);           % define Bx jump -+ 3nT and width 10s
Bz  = t*0;                            % define Bz as zero
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);   % define scalar TSeries object
TS1 = TS1 + 1*rand(TS1.length,3);     % add random noise of amplitude 1nT
TS1.units = 'nT';                     % specify units
TS1.userData.LABLAXIS = 'B';

h = irf_plot(1,'newfigure');		  	  % initialize figure with one panel
irf_plot(h,TS1);						          % plot times series  
irf_minvar_gui(TS1);                  % run minimum variance analysis on time series

%% Example XX, magnetopause crossing Rosenqvist2008
% http://www.cluster.rl.ac.uk/csdsweb-cgi/csdsweb_pick?P_TYPE=P1&YEAR=2001&MONTH=Jan&DAY=26&SUB_PLOT=S02

% tInt = irf.tint('2001-01-26T10:30:00Z/2001-01-26T11:00:00Z');
tInt = irf.tint('2001-07-05T04:30:00Z/2001-07-05T06:30:00Z');

%caa_download(tInt,'C1_CP_FGM_SPIN');
%caa_download(tInt,'C1_CP_EFW_L3_E3D_GSE');
caa_load C1_CP_FGM_SPIN
B1 = irf_get_data('B_vec_xyz_gse__C1_CP_FGM_SPIN','caa','ts');
irf_plot(1,'newfigure');
irf_plot(B1);
irf_minvar_gui(B1)

%% De Hofmann - Teller frame
