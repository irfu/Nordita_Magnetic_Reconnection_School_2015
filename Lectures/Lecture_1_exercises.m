
% Nordita 2015, School Data Analysis
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
% Lets generate 5samples/s time series in interval 2008-03-01 10:00-11:00 UT,
%  - create exponentially growing wave 
%	    y(t)=exp(0.001(t-to))*sin(t-to) 
%  - plot it.

T = irf_time([2008 03 01 10 0 0],'vector>epochtt') + (0:.2:3600);
t = T.t-T.t(1);
y = exp(0.001*t).*sin(2*pi*t/180);	
Y = TSeries(T,y);

% plot
h = irf_plot(1,'newfigure'); % generate new figure with 1 panel
irf_plot(h,Y);						   % plot data

%% Example  2 Plot multicomponent data
% Generate data with two components and plot in the same figure.
% Add legend text in lower left corner
% As you notice irfu-matlab interprets some common names for variables, 
% i.e. B2 is assumed to be magnetic field measurement by Cluster 2

z = exp(0.001*t).*cos(2*pi*t/180);	% z(t)=exp(0.001(t-to))*cos(t)
F = TSeries(T,[y z]);							          % B2 has two components, y & z
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
