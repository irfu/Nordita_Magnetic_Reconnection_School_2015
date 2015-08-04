%% Nordita 2015, School Data Analysis
% 
% Lecture 4 
% Multi-s/c boundary methods : timing
%
% 
% As real data we use Paschmann 2005 paper

Tint = irf.tint('2002-03-30T13:11:30Z/2002-03-30T13:12:00Z');
cd Event_20020330_1311
c_eval('B?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''ts'');') % Load TSeries of B GSE
irf_pl_tx B?


%% 4 s/c timing
c_4_v_gui B?

%% Error estimate
c_eval('R?=c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''ts'');') % Load TSeries of B GSE

Tint1 = irf.tint('2002-03-30T13:11:43Z/2002-03-30T13:11:49Z'); % time interval, from data

[V, dV] = c_4_v_xcorr(Tint1,B1,B2,B3,B4,R1,R2,R3,R4);

%% 2 SC timing
% perform on 2 s/c pairs: 1-4, then for 2-3, comopare with 4-s/c data

irf_minvar_gui(B1)
n1 = ud.v3;

irf_minvar_gui(B4)
n4 = ud.v3;

fprintf('Angle between two normals %.2f deg\n',acosd(n1*n4'))

c_4_v_gui B? % align C1 and C4

%% Assignment use data from Paschmann 2005 paper