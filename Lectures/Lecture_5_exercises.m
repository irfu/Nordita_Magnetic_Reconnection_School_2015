%% Nordita 2015, School Data Analysis, Multi-s/c methods, Curlometer, gradients. 
% 
% Huishan Fu poster during the workshop comparing different methods  


%% Example, multi s/c observations of current

T   = EpochTT('2002-03-04T09:30:00Z'):5 ...
     :EpochTT('2002-03-04T10:30:00Z');% define time line as EpochTT object
t   = T - T.start;                    % define relative time in s from start
t   = t - mean(t);                    % time zero in the middle of interval

Ljy = 500e3;                          % 500km, half width of jy current sheet
Ljx = 50e3;                           % 50km,  half width of jx current sheet
B0x = 10;                             % asymptotic Bx magnetic field [nT]
B0y = 3;                              % asymptotic By magnetic field [nT]
vz  = 1e3;                            % crossing current sheet at vz = 1km/s

mu0 = 4*pi/1e7;
% Define functions
B_  = @(x,y,z) [B0x*tanh(z/Ljy)     -B0y*tanh(z/Ljx)    0*x];
J_  = @(x,y,z) [B0y/Ljx*sech(z/Ljx).^2 B0x/Ljy*sech(z/Ljy).^2 0*x]/mu0;

L = 150e3;                % s/c separation scale [m]
Rconf.dr1 = [0 0 0];      % C1 relative locations
Rconf.dr2 = [L 0 L/3];    % C2 -=-
Rconf.dr3 = [0 L L/2];
Rconf.dr4 = [0 0 L];

Rref = irf.ts_vec_xyz(T,t*[0 0 vz]); % satellite moves in Z with vz
Rref.units = 'm';
R.C1 = Rref + Rconf.dr1;             % C1 position
R.C2 = Rref + Rconf.dr2;             % C2 position 
R.C3 = Rref + Rconf.dr3;
R.C4 = Rref + Rconf.dr4;
R.C  = Rref;                         % s/c tetrahedron center 
R.C.data = (R.C1.data+R.C2.data+R.C3.data+R.C4.data)/4;

clear B
B.C1 = R.C1;
B.C1.units = 'nT';
B.C1.userData.LABLAXIS = 'B';
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
curlB       = c_4_grad(R,B,'curl'); 
jCurlometer = curlB * (mu0^(-1));   % Current in units nA/m^2
J           = jCurlometer;
J.data      = J_(R.C.data(:,1),R.C.data(:,2),R.C.data(:,3)); % theoretical current

h = irf_plot({jCurlometer,J},'comp');
irf_zoom(h,'x',irf.tint(J.time.start,J.time.stop))  % all subplots the same time
ylabel(h(1),'J_X')
title(h(1),['Spacecraft separation ' num2str(L/1e3) ' km.']);
irf_legend(h(1),{'J_{curlometer}','J_{theory}'},[0.02 0.95]);
irf_legend(h(1),['Ljx=' num2str(Ljx/1e3) 'km.'] ,[0.98 0.95]);
ylabel(h(2),'J_Y')
irf_legend(h(2),['Ljy=' num2str(Ljy/1e3) 'km.'] ,[0.98 0.95]);
ylabel(h(3),'J_Z')

%% Divergence B
divB = c_4_grad(R,B,'div'); 
relErr = divB;              % divB/|curlB|
relErr.data = divB.data ./ curlB.abs.data;
h=irf_plot({divB,relErr});
title(h(1),'Multi-s/c analysis');
ylabel(h(1),'\nabla\cdot B','interpreter','tex')
ylabel(h(2),'\nabla\cdot B/\nabla\times B','interpreter','tex')
