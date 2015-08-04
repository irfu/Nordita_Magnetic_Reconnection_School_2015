%% Nordita 2015, School Data Analysis
% 
% Lecture 2 
% Single s/c boundary methods
%
% Minimum Variance Analysis
% De Hofmann-Teller frame
% 
% As real data we use Paschmann 2005 paper

%% Change to working directory
cd /Users/andris/Dropbox/Projects/Nordita2015/Data/CAA_20010705_0430_20010705_0630


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
Bx  =  5*tanh((t-1800)/600);          % define Bx +-5nT jump and 10min width
By  = -3*tanh((t-1800)/60);           % define Bx -+3nT jump and 1min width
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
Bx  =  5*tanh((t-1800)/600);          % define Bx jump +- 5nT and width 10 min
By  = -3*tanh((t-1800)/60);           % define Bx jump -+ 3nT and width 1 min
Bz  = t*0;                            % define Bz as zero
TS1 = irf.ts_vec_xyz(T,[Bx By Bz]);   % define scalar TSeries object
TS1 = TS1 + 1*rand(TS1.length,3);     % add random noise of amplitude 1nT
TS1.units = 'nT';                     % specify units
TS1.userData.LABLAXIS = 'B';

h = irf_plot(1,'newfigure');		  	  % initialize figure with one panel
irf_plot(h,TS1);						          % plot times series  
irf_minvar_gui(TS1);                  % run minimum variance analysis on time series

%% Example XX, magnetopause crossing Paschmann 2005
% http://www.cluster.rl.ac.uk/csdsweb-cgi/csdsweb_pick?P_TYPE=P1&YEAR=2001&MONTH=Jan&DAY=26&SUB_PLOT=S02

% tInt = irf.tint('2001-01-26T10:30:00Z/2001-01-26T11:00:00Z');
tInt = irf.tint('2001-07-05T04:30:00Z/2001-07-05T06:30:00Z');

%caa_download(tInt,'C1_CP_FGM_SPIN');
%caa_download(tInt,'C1_CP_FGM_FULL');
%caa_download(tInt,'C1_CP_EFW_L3_E3D_GSE');
%caa_download(tInt,'C1_CP_CIS_HIA_ONBOARD_MOMENTS');

caa_load C1_CP_FGM_SPIN
B1 = irf_get_data('B_vec_xyz_gse__C1_CP_FGM_SPIN','caa','ts');
irf_plot(1,'newfigure');
irf_plot(B1);
irf_minvar_gui(B1)

%% De Hofmann - Teller frame
B1 = irf_get_data('B_vec_xyz_gse__C1_CP_FGM_SPIN','caa','ts');
V1 = c_caa_var_get('velocity_gse__C1_CP_CIS_HIA_ONBOARD_MOMENTS','caa','ts');
E1 = c_caa_var_get('E_Vec_xyz_GSE__C1_CP_EFW_L3_E3D_GSE','caa','ts');
E1vxb = irf_e_vxb(V1,B1);
V1exb = irf_e_vxb(E1,B1,-1);

VHT = irf_vht(E1,B1);
EHT = irf_e_vxb(VHT,B1);

