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

%% Boundary (Harris like) crossing, change in Bx component
% Lets generate 1sample each 5s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

T   = EpochTT('2002-03-04T09:30:00Z'):5 ...
     :EpochTT('2002-03-04T10:30:00Z');% define time line as EpochTT object
t   = T - T.start;                    % define relative time in s from start
t   = t - mean(t);                    % time zero in the middle of interval

Bx  =  5*tanh(t/60);                  % define Bx +-5nT jump and 1min width
By  = t*0;                            % define Bx to be zero
Bz  = t*0;                            % define Bz to be zero

B = irf.ts_vec_xyz(T,[Bx By Bz]);     % define TSeries object (vector)
B.units = 'nT';                       % units
B.userData.LABLAXIS = 'B';            % plot label axis

h = irf_plot(1,'newfigure');		  	  % initialize figure with one panel
irf_plot(h,B);						            % plot times series  
irf_legend(h,{'Bx','By','Bz'},[0,1]);

%% Define Harris current sheet

mu0 = 4*pi/1e7;

Ljy = 500e3;                          % 500km, half width of current sheet
B0  = 10;                             % Asymptotic magnetic field [nT]
vz = 1e3;                             % crossing current sheet at vz = 1km/s

Bx  = B0*tanh(t*vz/Ljy);              % define Bx jump +- 5nT and width 10 min
By  = t*0;                            % define Bx jump -+ 2nT and width 1 min
Bz  = t*0;                            % define Bz as zero

Jx  = t*0;                            % Narrow current sheet in X
Jy  = B0/Ljy*sech(t*vz/Ljy)/mu0;      % Harris current sheet in Y
Jz  = t*0;                            % zero current in Z directoin

B   = irf.ts_vec_xyz(T,[Bx By Bz]);   % define B as TSeries
J   = irf.ts_vec_xyz(T,[Jx Jy Jz]);   % define J as TSeries

h=irf_plot({B,J});
ylabel(h(1),'B [nT]')
ylabel(h(2),'J [n A/m^2]','interpreter','tex')
irf_legend(h(1),{'Bx','By','Bz'},[0.02,0.98]);
irf_legend(h(2),{'Jx','Jy','Jz'},[0.02,0.98]);

%% Double Harris current sheet, B change in two components

Ljy = 500e3;                          % 500km, half width of jy current sheet
Ljx = 50e3;                           % 50km,  half width of jx current sheet
B0x = 10;                             % Asymptotic Bx magnetic field [nT]
B0y = 3;                              % Asymptotic By magnetic field [nT]
vz = 1e3;                             % crossing current sheet at vz = 1km/s

Bx  = B0x*tanh(t*vz/Ljy);             % define Bx jump
By  = B0y*tanh(t*vz/Ljx);             % define By jump
Bz  = t*0;                            % define Bz as zero

Jx  = B0y/Ljx*sech(t*vz/Ljx).^2/mu0;   % Current sheet jx
Jy  = B0x/Ljy*sech(t*vz/Ljy).^2/mu0;   % Current sheet jy
Jz  = t*0;                            % zero current in Z directoin

B   = irf.ts_vec_xyz(T,[Bx By Bz]);   % define B as TSeries
J   = irf.ts_vec_xyz(T,[Jx Jy Jz]);   % define J as TSeries

h=irf_plot({B,J});
ylabel(h(1),'B [nT]')
ylabel(h(2),'J [n A/m^2]','interpreter','tex')
irf_legend(h(1),{'Bx','By','Bz'},[0.02,0.98]);
irf_legend(h(2),{'Jx','Jy','Jz'},[0.02,0.98]);

%% Minimum variance
irf_minvar_gui(B);                    % run minimum variance analysis on time series

Bgsm = irf_gse2gsm(B);
irf_minvar_gui(Bgsm);                 % run minimum variance analysis on time series
 
%%  Double Harris current sheet, B change in two components + noise
% Lets generate 5samples/s time series during 1h after 2002-03-04 09:30 UTC,
% including artificial boundary in the middle of the interval.

B0noise = 3;                          % noise amplitude
noisyTimeSeries = model.synthetic_time_series(...
	'fs',1,'f',0,'peakHalfWidth',0.3,...
	'timeInterval',1e4,'components',3);
Bnoisy = B + B0noise*noisyTimeSeries(1:B.length,2:4);% add random noise of amplitude 1nT
Bnoisy.units = 'nT';                  % specify units
Bnoisy.userData.LABLAXIS = 'B';


irf_minvar_gui(Bnoisy);              % run minimum variance analysis on time series

%% Example, multi s/c observations of current

% Define functions
B_ = @(x,y,z) [B0x*tanh(z/Ljy)     -B0y*tanh(z/Ljx)    0*x];
J_ = @(x,y,z) [B0y/Ljx*sech(z/Ljx).^2 B0x/Ljy*sech(z/Ljy).^2 0*x]/mu0;

L = 1250e3;               % s/c separation scale [m]
Rconf.dr1 = [0 0 0];    % C1 relative locations
Rconf.dr2 = [L 0 L/3];
Rconf.dr3 = [0 L L/2];
Rconf.dr4 = [0 0 L];

Rref = irf.ts_vec_xyz(T,t*[0 0 vz]); % satellite moves in Z with vz
Rref.units = 'm';
R.C1 = Rref + Rconf.dr1;
R.C2 = Rref + Rconf.dr2;
R.C3 = Rref + Rconf.dr3;
R.C4 = Rref + Rconf.dr4;
R.C  = Rref;
R.C.data = (R.C1.data+R.C2.data+R.C3.data+R.C4.data)/4;

clear B
B.C1 = R.C1;
B.C1.units = 'nT';B.C1.userData.LABLAXIS = 'B';
B.C1.data = B_(R.C1.data(:,1),R.C1.data(:,2),R.C1.data(:,3));
B.C2 = B.C1;
B.C3 = B.C1;
B.C4 = B.C1;
B.C2.data = B_(R.C2.data(:,1),R.C2.data(:,2),R.C2.data(:,3));
B.C3.data = B_(R.C3.data(:,1),R.C3.data(:,2),R.C3.data(:,3));
B.C4.data = B_(R.C4.data(:,1),R.C4.data(:,2),R.C4.data(:,3));

h=irf_pl_tx(B);
ylabel(h(1),'Bx [nT]');
ylabel(h(2),'By [nT]');
ylabel(h(3),'Bz [nT]');

%% Curlometer, current from 4 s/c measurements
curlB = c_4_grad(R,B,'curl');
jCurlometer = curlB * (mu0^(-1));
J = jCurlometer;
J.data = J_(R.C.data(:,1),R.C.data(:,2),R.C.data(:,3));
h=irf_plot({jCurlometer,J},'comp');
irf_zoom(h,'x',irf.tint(J.time.start,J.time.stop))  % all subplots the same time

%% Divergence B
divB = c_4_grad(R,B,'div');
jDiv = divB * (mu0^(-1));
relErr = jDiv;
relErr.data = jDiv.data ./ jCurlometer.abs.data;
irf_plot({jDiv,relErr});

%% De Hoffmann - Teller frame
vSpacecraft = [0 0 vz];
E.C1 = irf_e_vxb(vSpacecraft,B.C1);
h=irf_plot({E.C1,B.C1});
irf_legend(h(1),{'Ex','Ey','Ez'},[0.02,0.98]);
irf_legend(h(2),{'Bx','By','Bz'},[0.02,0.98]);
VHT = irf_vht(E.C1,B.C1);

%% De Hoffmann - Teller frame, 2 E field components

EE.C1 = E.C1;
EE.C1.data = E.C1.data + (rand(size(E.C1.data))-0.5);
BB.C1 = B.C1;
BB.C1.data = B.C1.data + (rand(size(B.C1.data))-0.5);
h=irf_plot({EE.C1,BB.C1});
VHT = irf_vht(EE.C1,BB.C1,2);


%% Example XX, magnetopause crossing Paschmann 2005
% http://www.cluster.rl.ac.uk/csdsweb-cgi/csdsweb_pick?P_TYPE=P1&YEAR=2001&MONTH=Jan&DAY=26&SUB_PLOT=S02

% tInt = irf.tint('2001-01-26T10:30:00Z/2001-01-26T11:00:00Z');
tInt = irf.tint('2001-07-05T04:30:00Z/2001-07-05T06:30:00Z');

if 0, 
	caa_download(tInt,'C1_CP_FGM_SPIN');
	caa_download(tInt,'C1_CP_FGM_5VPS');
	caa_download(tInt,'C1_CP_FGM_FULL');
	caa_download(tInt,'C1_CP_EFW_L2_E3D_GSE');
	caa_download(tInt,'C1_CP_EFW_L3_E3D_GSE');
	caa_download(tInt,'C1_CP_CIS_HIA_ONBOARD_MOMENTS');
	caa_download(tInt,'C1_CP_CIS_HIA_HS_1D_PEF');
	caa_download(tInt,'C1_CP_RAP_ESPCT6');
	caa_download(tInt,'C1_CP_PEA_PITCH_SPIN_DEFlux');
end 

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


