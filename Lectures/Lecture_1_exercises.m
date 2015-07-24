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

%% Example 1  artifical times series
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

%% Example  time interval
Tint = irf.tint('2002-03-04T09:40:00Z/2002-03-04T10:20:00Z');
irf_zoom(h,'x',Tint)

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
title(h(1),T(1).utc('yyyy-mm-dd'));
ylabel(h(2),'F_Y');
irf_legend(h(1),{'FNew','Fnew=F*2+2'},[0.02 0.98])

Hours=TSeries(T,(T.t-T.t(1))/3600); % hours from the beginning of the time interval
irf_timeaxis(h(2),Hours,{'hours'})
irf_timeaxis(h(end),'nodate');
