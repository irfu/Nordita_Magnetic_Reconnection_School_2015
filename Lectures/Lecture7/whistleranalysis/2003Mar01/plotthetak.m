% Calculates theta_k over short intervals using MVA and plots data

ic = 4; %Spacecraft number 1--4

Bisr2 = c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');
Bibm=c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','mat');
SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');
Bibm(:,2) = Bibm(:,2)-nanmean(Bibm(:,2));
Bibm(:,3) = Bibm(:,3)-nanmean(Bibm(:,3));
Bibm(:,4) = Bibm(:,4)-nanmean(Bibm(:,4));
Bibm(isnan(Bibm)) = 0;

[~,idx]=irf_tlim(Bisr2,[Bibm(1,1) Bibm(end,1)]);
Bisr2 = Bisr2(idx,:);
Bisr2 = irf_resamp(Bisr2,Bibm);

Bfac = irf_convert_fac(Bibm,Bisr2,SCpos);

%define intervals 
time = Bibm(end,1)-Bibm(1,1);
dt = time/length(Bibm(:,1));
whistfreq = 100;

% set each interval length to 10 wave periods
inttime = 10.0/whistfreq;
numtvals = floor(inttime/dt);
numintervals = floor(time/inttime);

%define arrays of theta and l_int/l_min
thetak = zeros(numintervals,2);
lintlmin = zeros(numintervals,2);


for ii=1:numintervals;
    Bibms = Bibm([((ii-1)*numtvals+1):1:ii*numtvals],:);
    Bisr2s = Bisr2([((ii-1)*numtvals+1):1:ii*numtvals],:);
    Bibms(:,2) = Bibms(:,2)-nanmean(Bibms(:,2));
    Bibms(:,3) = Bibms(:,3)-nanmean(Bibms(:,3));
    Bibms(:,4) = Bibms(:,4)-nanmean(Bibms(:,4));
    Bmag = sqrt(Bisr2s(:,2).^2+Bisr2s(:,3).^2+Bisr2s(:,4).^2);
    Bmag = nanmean(Bmag);
    Bvec = nanmean(Bisr2s(:,[2:4]))/Bmag;
    
    %min variance analysis of the short interval
    MVArow1 = [mean(Bibms(:,2).^2)-mean(Bibms(:,2))^2 mean(Bibms(:,2).*Bibms(:,3))-mean(Bibms(:,2))*mean(Bibms(:,3)) mean(Bibms(:,2).*Bibms(:,4))-mean(Bibms(:,2))*mean(Bibms(:,4))];
    MVArow2 = [mean(Bibms(:,3).*Bibms(:,2))-mean(Bibms(:,3))*mean(Bibms(:,2)) mean(Bibms(:,3).^2)-mean(Bibms(:,3))^2 mean(Bibms(:,3).*Bibms(:,4))-mean(Bibms(:,3))*mean(Bibms(:,4))];
    MVArow3 = [mean(Bibms(:,4).*Bibms(:,2))-mean(Bibms(:,4))*mean(Bibms(:,2)) mean(Bibms(:,4).*Bibms(:,3))-mean(Bibms(:,4))*mean(Bibms(:,3)) mean(Bibms(:,4).^2)-mean(Bibms(:,4))^2];    
    MVAmat = [MVArow1; MVArow2; MVArow3];

    [MVAvecs,MVAeigs] = eig(MVAmat);
    vmin = MVAvecs(:,1);
    
    thetaktemp = acosd(abs(Bvec*vmin));
    lintlmintemp = MVAeigs(2,2)/MVAeigs(1,1);
    ttemp = mean(Bibms(:,1));
    thetak(ii,1) = ttemp;
    thetak(ii,2) = thetaktemp;
    lintlmin(ii,1) = ttemp;
    lintlmin(ii,2) = lintlmintemp;
end
    
h=irf_plot(4,'newfigure'); 

h(1)=irf_panel('Bisr2');
irf_plot(h(1),Bibm);
ylabel(h(1),'B_{ISR2} (nT)');
irf_legend(h(1),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Bfac');
irf_plot(h(2),Bfac);
ylabel(h(2),'B_{FAC} (nT)');
irf_legend(h(2),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('theta');
irf_plot(h(3),thetak);
ylabel(h(3),'\theta_{k} (deg)');
irf_zoom(h(3),'y',[0 90]);
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

tint1 = [Bibm(1,1) Bibm(end,1)];

h(4)=irf_panel('lintlmin');
irf_plot(h(4),lintlmin);
hold(h(4),'on');
irf_plot([tint1(1) 10; tint1(2) 10],'r');
hold(h(4),'off');
ylabel(h(4),'\lambda_{int}/\lambda_{min}');
set(h(4),'yscale','log');
irf_zoom(h(4),'y',[1 200]);
set(h(4),'ytick',[1e0 1e1 1e2 1e3]);
irf_legend(h(4),'(d)',[0.99 0.98],'color','k');

irf_plot_axis_align(1,h(1:4));
irf_zoom(h(1:4),'x',tint1);

xwidth = 0.86;
ywidth = 0.18;
set(h(1),'position',[0.12 0.99-ywidth xwidth ywidth]);
set(h(2),'position',[0.12 0.99-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.12 0.99-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.12 0.99-4*ywidth xwidth ywidth]);

irf_timeaxis(h(1:3),'nodate');
