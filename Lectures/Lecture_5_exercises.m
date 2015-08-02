%% Nordita 2015, School Data Analysis
% 
% Lecture 5
% 5. Multi-s/c methods. Curlometer, gradients. (Andris)
%
% Gradient
% Curlometer
% Divergence
% Magnetic null
% 
% Huishan Fu poster during workshop comparing different methods  

%% Example 1 Generate multi-s/c time series of Harris current sheet
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

%% Make virtual satellite orbits through the current sheet for different s/c separation

%% Magnetic field gradient

%% Current

%% Divergence 

%% Generate 3D null and measure it

