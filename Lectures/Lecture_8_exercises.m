%% Nordita 2015, School Data Analysis
% 
% Lecture 4 
% Using MMS QL data
% data files at Dropbox:
% https://www.dropbox.com/sh/5xwkqynebj628xn/AAB76xlJJ9RnlFkchTQ2foVVa?dl=0

% Load and Plot B data (DFG)
c_eval('Mms?_dfg_srvy_ql=dataobj(''mms?_dfg_srvy_ql*.cdf'');')
c_eval(['B? =mms.variable2ts(get_variable(Mms?_dfg_srvy_ql,''mms?_dfg_srvy_dmpa''));'...
  'B?.userData.LABLAXIS=''B''; B?.coordinateSystem=''dmpa'';'])

h = irf_pl_tx('B?');
irf_zoom(h,'x',irf.tint('2015-07-31T21:00:00Z/2015-07-31T23:00:00Z'))
irf_zoom(h,'y')
title(h(1),'MMS QL B DFG')
ylabel(h(1),'Bx [nT]'), ylabel(h(2),'By [nT]'), ylabel(h(3),'Bz [nT]')

%% Load and Plot E data (EDP=SDP+ADP)
c_eval('Mms?_edp_fast_ql=dataobj(''mms?_edp_fast_ql*.cdf'');')
c_eval(['E? =mms.variable2ts(get_variable(Mms?_edp_fast_ql,''mms?_edp_dce_xyz_dsl''));'...
  'E?.userData.LABLAXIS=''E''; E?.coordinateSystem=''dsl'';'])

h = irf_pl_tx('E?');
irf_zoom(h,'x',irf.tint('2015-07-31T02:40:00Z/2015-07-31T03:30:00Z'))
irf_zoom(h,'y')
title(h(1),'MMS QL E SDP+ADP')
ylabel(h(1),'Ex [mV/m]'), ylabel(h(2),'Ey [mV/m]'), ylabel(h(3),'Ez [mV/m]')

%% Plot E&B for MMS3, pilolarization event
h = irf_plot({E3.y,B3.z});
irf_zoom(h,'x',irf.tint('2015-07-31T21:00:00Z/2015-07-31T22:30:00Z'))
irf_zoom(h,'y')
title(h(1),'MMS3 QL E & B DSL/DMPA')
ylabel(h(1),'Ey [mV/m]'), ylabel(h(2),'Bz [nT]')