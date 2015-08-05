% Calculates theta_k over short intervals using MVA and plots data

ic = 2; %Spacecraft number 1--4

Bisr2 = c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');
Bibm=c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','ts');
SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');

Bibm.data(:,1) = Bibm.data(:,1)-nanmean(Bibm.data(:,1));
Bibm.data(:,2) = Bibm.data(:,2)-nanmean(Bibm.data(:,2));
Bibm.data(:,3) = Bibm.data(:,3)-nanmean(Bibm.data(:,3));
Bibm.data(isnan(Bibm.data)) = 0;

tlimit = irf.tint(Bibm.time.start.utc,Bibm.time.stop.utc);
Bisr2 = Bisr2.tlim(tlimit);
SCpos = SCpos.tlim(tlimit);
Bisr2 = Bisr2.resample(Bibm.time);
SCpos = SCpos.resample(Bibm.time);

%rotate entire burst mode interval into field-aligned coordinates
Bfac = Bibm;
Bmag = sqrt(Bisr2.data(:,1).^2+Bisr2.data(:,2).^2+Bisr2.data(:,3).^2);
Rpar = Bisr2.data./[Bmag Bmag Bmag];
Rperpy = irf_cross(Rpar,SCpos.data);
Rmag = sqrt(Rperpy(:,1).^2+Rperpy(:,2).^2+Rperpy(:,3).^2);
Rperpy = Rperpy./[Rmag Rmag Rmag];
Rperpx = irf_cross(Rperpy, Rpar);
Rmag = sqrt(Rperpx(:,1).^2+Rperpx(:,2).^2+Rperpx(:,3).^2);
Rperpx = Rperpx./[Rmag Rmag Rmag];

Bfac.data(:,3)=  Rpar(:,1).*Bibm.data(:,1)+  Rpar(:,2).*Bibm.data(:,2)+  Rpar(:,3).*Bibm.data(:,3);
Bfac.data(:,1)=Rperpx(:,1).*Bibm.data(:,1)+Rperpx(:,2).*Bibm.data(:,2)+Rperpx(:,3).*Bibm.data(:,3);
Bfac.data(:,2)=Rperpy(:,1).*Bibm.data(:,1)+Rperpy(:,2).*Bibm.data(:,2)+Rperpy(:,3).*Bibm.data(:,3);

%define intervals 
time = irf_time(tlimit(2),'epochtt>epoch')-irf_time(tlimit(1),'epochtt>epoch');
dt = time/length(Bibm.data(:,1));

%change this to approximate whistler frequency
whistfreq = 100;

% set each interval length to 10 wave periods, can change this different
% values to see how MVA analysis is affected
inttime = 10.0/whistfreq;
numtvals = floor(inttime/dt);
numintervals = floor(time/inttime);

%define arrays of theta and l_int/l_min
thetak = zeros(numintervals,1);
lintlmin = zeros(numintervals,1);

smallttseries = Bibm.time([floor(numtvals/2)+(0:numintervals-1)*numtvals]);

for ii=1:numintervals;
    Bibms = Bibm.data([((ii-1)*numtvals+1):1:ii*numtvals],:);
    Bisr2s = Bisr2.data([((ii-1)*numtvals+1):1:ii*numtvals],:);
    Bibms(:,1) = Bibms(:,1)-nanmean(Bibms(:,1));
    Bibms(:,2) = Bibms(:,2)-nanmean(Bibms(:,2));
    Bibms(:,3) = Bibms(:,3)-nanmean(Bibms(:,3));
    Bmag = sqrt(Bisr2s(:,1).^2+Bisr2s(:,2).^2+Bisr2s(:,3).^2);
    Bmag = nanmean(Bmag);
    Bvec = nanmean(Bisr2s)/Bmag;
    
    %min variance analysis of the short interval
    MVArow1 = [nanmean(Bibms(:,1).^2)-nanmean(Bibms(:,1))^2 nanmean(Bibms(:,1).*Bibms(:,2))-nanmean(Bibms(:,1))*nanmean(Bibms(:,2)) ...
        nanmean(Bibms(:,1).*Bibms(:,3))-nanmean(Bibms(:,1))*nanmean(Bibms(:,3))];
    MVArow2 = [nanmean(Bibms(:,2).*Bibms(:,1))-nanmean(Bibms(:,1))*nanmean(Bibms(:,1)) nanmean(Bibms(:,2).^2)-nanmean(Bibms(:,2))^2 ...
        nanmean(Bibms(:,2).*Bibms(:,3))-nanmean(Bibms(:,2))*nanmean(Bibms(:,3))];
    MVArow3 = [nanmean(Bibms(:,3).*Bibms(:,1))-nanmean(Bibms(:,3))*nanmean(Bibms(:,1)) ... 
        nanmean(Bibms(:,3).*Bibms(:,2))-nanmean(Bibms(:,3))*nanmean(Bibms(:,2)) nanmean(Bibms(:,3).^2)-nanmean(Bibms(:,3))^2];    
    MVAmat = [MVArow1; MVArow2; MVArow3];

    [MVAvecs,MVAeigs] = eig(MVAmat);
    vmin = MVAvecs(:,1);
    
    thetaktemp = acosd(abs(Bvec*vmin));
    lintlmintemp = MVAeigs(2,2)/MVAeigs(1,1);
    thetak(ii) = thetaktemp;
    lintlmin(ii) = lintlmintemp;
end
    
thetak = TSeries(smallttseries,thetak);
lintlmin = TSeries(smallttseries,lintlmin);

h=irf_plot(4,'newfigure'); 

h(1)=irf_panel('Bisr2');
irf_plot(h(1),Bibm);
ylabel(h(1),'B_{ISR2} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Bfac');
irf_plot(h(2),Bfac);
ylabel(h(2),'B_{FAC} (nT)','Interpreter','tex');
irf_legend(h(2),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('theta');
irf_plot(h(3),thetak);
ylabel(h(3),'\theta_{k} (deg)','Interpreter','tex');
irf_zoom(h(3),'y',[0 90]);
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('lintlmin');
irf_plot(h(4),lintlmin);
ylabel(h(4),'\lambda_{int}/\lambda_{min}','Interpreter','tex');
set(h(4),'yscale','log');
irf_zoom(h(4),'y',[1 200]);
set(h(4),'ytick',[1e0 1e1 1e2 1e3]);
irf_legend(h(4),'(d)',[0.99 0.98],'color','k');


irf_plot_axis_align(1,h(1:4));
irf_zoom(h(1:4),'x',tlimit);

xwidth = 0.86;
ywidth = 0.22;
set(h(1),'position',[0.12 0.99-ywidth xwidth ywidth]);
set(h(2),'position',[0.12 0.99-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.12 0.99-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.12 0.99-4*ywidth xwidth ywidth]);

%irf_timeaxis(h(1:4),'nodate');
