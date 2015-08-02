%% Nordita 2015, School Data Analysis
% 
% Lecture 3 
% Fluid equations, single s/c methods. Wal?n test. (Andris) 
%
% Walen test
% Hall fields
% 

%% Example 1 ideal rotational discontinuity
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

%% Example XX Magnetopause reconnection

%% Example Hall fields and boundary speed

