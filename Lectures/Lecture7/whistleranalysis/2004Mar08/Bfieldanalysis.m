% Script to plot search coil magnetic fields from L2 internal burst mode.
% Fields are transformed into field-aligned and minimum variance
% coordinates. Estimates wave-normal angle from minimum variance analysis.
% Plots time series of B and hodograms.

ic = 4; %Spacecraft number 1--4

%Read in Cluster data: Search coil B, FGM B, and SC position.
Bisr2 = c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');
SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');
Bibm=c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','ts');

fprintf('Interval burst mode interval \n')
fprintf(strcat(Bibm.time.start.utc,'--',Bibm.time.stop.utc,'\n\n'))
tlimit = irf.tint(Bibm.time.start.utc,Bibm.time.stop.utc);
Bisr2 = Bisr2.tlim(tlimit);
SCpos = SCpos.tlim(tlimit);

%Resample data to the same frequency as the search coil
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

% select time interval of interest
tint1 = EpochUnix([iso2epoch('2004-03-08T15:55:33.20Z') iso2epoch('2004-03-08T15:55:33.80Z')]);

%extract shorter time interval
Bibms=Bibm.tlim(tint1);
Bfacs=Bfac.tlim(tint1);


%remove background offset fields
Bibms.data(:,1) = Bibms.data(:,1)-nanmean(Bibms.data(:,1));
Bibms.data(:,2) = Bibms.data(:,2)-nanmean(Bibms.data(:,2));
Bibms.data(:,3) = Bibms.data(:,3)-nanmean(Bibms.data(:,3));
Bfacs.data(:,1) = Bfacs.data(:,1)-nanmean(Bfacs.data(:,1));
Bfacs.data(:,2) = Bfacs.data(:,2)-nanmean(Bfacs.data(:,2));
Bfacs.data(:,3) = Bfacs.data(:,3)-nanmean(Bfacs.data(:,3));

% Calculate average direction of B over the short time interval tint1
Bisr2s=Bisr2.tlim(tint1);
Bmag = sqrt(Bisr2s.data(:,1).^2+Bisr2s.data(:,2).^2+Bisr2s.data(:,3).^2);
Bmag = nanmean(Bmag);
Bvec = nanmean(Bisr2s.data(:,[1:3]))/Bmag;


% Simple minimum variance analysis
% Construct MVA matrix
MVArow1 = [nanmean(Bibms.data(:,1).^2)-nanmean(Bibms.data(:,1))^2 nanmean(Bibms.data(:,1).*Bibms.data(:,2))-nanmean(Bibms.data(:,1))*nanmean(Bibms.data(:,2)) ...
    nanmean(Bibms.data(:,1).*Bibms.data(:,3))-nanmean(Bibms.data(:,1))*nanmean(Bibms.data(:,3))];
MVArow2 = [nanmean(Bibms.data(:,2).*Bibms.data(:,1))-nanmean(Bibms.data(:,1))*nanmean(Bibms.data(:,1)) nanmean(Bibms.data(:,2).^2)-nanmean(Bibms.data(:,2))^2 ...
    nanmean(Bibms.data(:,2).*Bibms.data(:,3))-nanmean(Bibms.data(:,2))*nanmean(Bibms.data(:,3))];
MVArow3 = [nanmean(Bibms.data(:,3).*Bibms.data(:,1))-nanmean(Bibms.data(:,3))*nanmean(Bibms.data(:,1)) ... 
    nanmean(Bibms.data(:,3).*Bibms.data(:,2))-nanmean(Bibms.data(:,3))*nanmean(Bibms.data(:,2)) nanmean(Bibms.data(:,3).^2)-nanmean(Bibms.data(:,3))^2];    
MVAmat = [MVArow1; MVArow2; MVArow3];

[MVAvecs,MVAeigs] = eig(MVAmat);
vmin = MVAvecs(:,1);
vint = MVAvecs(:,2);
vmax = MVAvecs(:,3);
MVAeigs = [MVAeigs(1,1) MVAeigs(2,2) MVAeigs(3,3)]; 

fprintf(strcat('MVA lambda_max/lambda_int = ',num2str(MVAeigs(3)/MVAeigs(2)),'\n'))
fprintf(strcat('MVA lambda_int/lambda_min = ',num2str(MVAeigs(2)/MVAeigs(1)),'\n'))

%Calculate theta_k 
thetak = acosd(abs(Bvec*vmin));
fprintf(strcat('Wave normal angle theta_k = ',num2str(thetak),'\n'))


% Coordinate transformation from ISR2 coordinates to Min Var coordinate
% system
Bibmmva = Bibms;
Bibmmva.data(:,1) = Bibms.data(:,1)*vmax(1)+Bibms.data(:,2)*vmax(2)+Bibms.data(:,3)*vmax(3); 
Bibmmva.data(:,2) = Bibms.data(:,1)*vint(1)+Bibms.data(:,2)*vint(2)+Bibms.data(:,3)*vint(3);
Bibmmva.data(:,3) = Bibms.data(:,1)*vmin(1)+Bibms.data(:,2)*vmin(2)+Bibms.data(:,3)*vmin(3);



%some plots of magnetic field waveforms

%Hodograms in Minimum variance coordinates
if 1, 
fn=figure;
set(fn,'Position',[10 400 600 200])
    h(1)=axes('position',[0.07 0.2 0.25 0.75]); % [x y dx dy]
    h(2)=axes('position',[0.4 0.2 0.25 0.75]); % [x y dx dy]
    h(3)=axes('position',[0.73 0.2 0.25 0.75]); % [x y dx dy]
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 

plot(h(1),Bibmmva.data(:,3),Bibmmva.data(:,1));
ylabel(h(1),'B_{max} (nT)');
xlabel(h(1),'B_{min} (nT)');
irf_zoom(h(1),'x',[-0.3 0.3])
irf_zoom(h(1),'y',[-0.3 0.3])

plot(h(2),Bibmmva.data(:,3),Bibmmva.data(:,2));
ylabel(h(2),'B_{int} (nT)');
xlabel(h(2),'B_{min} (nT)');
irf_zoom(h(2),'x',[-0.3 0.3])
irf_zoom(h(2),'y',[-0.3 0.3])

plot(h(3),Bibmmva.data(:,2),Bibmmva.data(:,1));
ylabel(h(3),'B_{max} (nT)');
xlabel(h(3),'B_{int} (nT)');
irf_zoom(h(3),'x',[-0.3 0.3])
irf_zoom(h(3),'y',[-0.3 0.3])

if 0, %set to 1 to save figure
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters minvarhodograms.eps;
end
 
end

% B waveform in ISR2 coordinates and field-aligned coordiantes. 
% Panel (c)

if 1,
h=irf_plot(3,'newfigure'); 

h(1)=irf_panel('Bisr2');
irf_plot(h(1),Bibms);
ylabel(h(1),'B_{ISR2} (nT)');
irf_legend(h(1),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Bfac');
irf_plot(h(2),Bfacs);
ylabel(h(2),'B_{FAC} (nT)');
irf_legend(h(2),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('Bfac2');
irf_plot(h(3),Bfacs);
ylabel(h(3),'B_{FAC} (nT)');
irf_legend(h(3),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

xwidth = 0.86;
ywidth = 0.18;

tint2 = EpochUnix([iso2epoch('2004-03-08T15:55:33.45Z') iso2epoch('2004-03-08T15:55:33.51Z')]);

irf_plot_axis_align(1,h(1:3))
irf_zoom(h(1:2),'x',tint1);
irf_zoom(h(3),'x',tint2);

set(h(1),'position',[0.12 0.99-ywidth xwidth ywidth]);
set(h(2),'position',[0.12 0.99-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.12 0.99-3*ywidth-0.05 xwidth ywidth]);

irf_pl_mark(h(1:2),irf_time(tint2,'epochtt>epoch')',[255 255 0]/255)
irf_plot_zoomin_lines_between_panels(h(2),h(3));

irf_timeaxis(h(1:2),'nodate');

if 0,
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters whistlerfac.eps;
end

end

% Hodogram in field-aligned coordinates
if 1,
fn=figure;
set(fn,'Position',[10 400 250 250])
    h(1)=axes('position',[0.2 0.2 0.75 0.75]); % [x y dx dy]
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 

plot(h(1),Bfacs.data(:,2),Bfacs.data(:,1),Bfacs.data(1,2),Bfacs.data(1,1),'ko',Bfacs.data(end,2),Bfacs.data(end,1),'kx');
ylabel(h(1),'B_x (nT)');
xlabel(h(1),'B_y (nT)');

if 0,
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters whistlerhodogram.eps;
end

end
