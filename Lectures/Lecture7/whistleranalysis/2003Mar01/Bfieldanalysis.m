ic = 2; %Spacecraft number 1--4

Bisr2 = c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');
SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');
Bibm=c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','mat');

bmint = irf_time([Bibm(1,1), Bibm(end,1)],'epoch>iso');
fprintf('Interval burst mode interval \n')
fprintf(strcat(bmint(1,:),'--',bmint(2,:),'\n'))

Bisr2=irf_tlim(Bisr2,[Bibm(1,1) Bibm(end,1)]);
Bisr2 = irf_resamp(Bisr2,Bibm);

%set undefined points to zero
Bibm(isnan(Bibm)) = 0;

% select time interval of interest
tint1 = [irf_time([2003 03 01 04 02 20.5],'vector6>epoch') irf_time([2003 03 01 04 02 21.5],'vector6>epoch')]; % short time interval within internal burst interval
%tint1 = [Bibm(1,1) Bibm(end,1)]; %entire internal burst mode interval



%extract shorter time interval
Bibms=irf_tlim(Bibm,tint1);

%remove background field
Bibms(:,2) = Bibms(:,2)-nanmean(Bibms(:,2));
Bibms(:,3) = Bibms(:,3)-nanmean(Bibms(:,3));
Bibms(:,4) = Bibms(:,4)-nanmean(Bibms(:,4));

Bisr2s=irf_tlim(Bisr2,tint1);

Bmag = sqrt(Bisr2s(:,2).^2+Bisr2s(:,3).^2+Bisr2s(:,4).^2);
Bmag = nanmean(Bmag);

% Average vector of B over time interval tint
Bvec = nanmean(Bisr2s(:,[2:4]))/Bmag;

% Simple minimum variance analysis
% Construct MVA matrix
MVArow1 = [mean(Bibms(:,2).^2)-mean(Bibms(:,2))^2 mean(Bibms(:,2).*Bibms(:,3))-mean(Bibms(:,2))*mean(Bibms(:,3)) mean(Bibms(:,2).*Bibms(:,4))-mean(Bibms(:,2))*mean(Bibms(:,4))];
MVArow2 = [mean(Bibms(:,3).*Bibms(:,2))-mean(Bibms(:,3))*mean(Bibms(:,2)) mean(Bibms(:,3).^2)-mean(Bibms(:,3))^2 mean(Bibms(:,3).*Bibms(:,4))-mean(Bibms(:,3))*mean(Bibms(:,4))];
MVArow3 = [mean(Bibms(:,4).*Bibms(:,2))-mean(Bibms(:,4))*mean(Bibms(:,2)) mean(Bibms(:,4).*Bibms(:,3))-mean(Bibms(:,4))*mean(Bibms(:,3)) mean(Bibms(:,4).^2)-mean(Bibms(:,4))^2];    
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
Bibmmva(:,2) = Bibms(:,2)*vmax(1)+Bibms(:,3)*vmax(2)+Bibms(:,4)*vmax(3); 
Bibmmva(:,3) = Bibms(:,2)*vint(1)+Bibms(:,3)*vint(2)+Bibms(:,4)*vint(3);
Bibmmva(:,4) = Bibms(:,2)*vmin(1)+Bibms(:,3)*vmin(2)+Bibms(:,4)*vmin(3);


%Rotate into field aligned coordinates
Bfac = irf_convert_fac(Bibms,Bisr2,SCpos);

Bfacs = irf_tlim(Bfac,tint1);


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

plot(h(1),Bibmmva(:,4),Bibmmva(:,2));
ylabel(h(1),'B_{max} (nT)');
xlabel(h(1),'B_{min} (nT)');
irf_zoom(h(1),'x',[-0.3 0.3])
irf_zoom(h(1),'y',[-0.3 0.3])

plot(h(2),Bibmmva(:,4),Bibmmva(:,3));
ylabel(h(2),'B_{int} (nT)');
xlabel(h(2),'B_{min} (nT)');
irf_zoom(h(2),'x',[-0.3 0.3])
irf_zoom(h(2),'y',[-0.3 0.3])

plot(h(3),Bibmmva(:,3),Bibmmva(:,2));
ylabel(h(3),'B_{max} (nT)');
xlabel(h(3),'B_{int} (nT)');
irf_zoom(h(3),'x',[-0.3 0.3])
irf_zoom(h(3),'y',[-0.3 0.3])

 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters minvarhodograms.eps;

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
irf_plot(h(2),Bfac);
ylabel(h(2),'B_{FAC} (nT)');
irf_legend(h(2),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('Bfac2');
irf_plot(h(3),Bfac);
ylabel(h(3),'B_{FAC} (nT)');
irf_legend(h(3),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

xwidth = 0.86;
ywidth = 0.18;

tint2 = [irf_time([2003 03 01 04 02 20.8],'vector6>epoch') irf_time([2003 03 01 04 02 21.0],'vector6>epoch')];

irf_plot_axis_align(1,h(1:3))
irf_zoom(h(1:2),'x',tint1);
irf_zoom(h(3),'x',tint2);

set(h(1),'position',[0.12 0.99-ywidth xwidth ywidth]);
set(h(2),'position',[0.12 0.99-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.12 0.99-3*ywidth-0.05 xwidth ywidth]);

irf_pl_mark(h(1:2),tint2,[255 255 0]/255)
irf_plot_zoomin_lines_between_panels(h(2),h(3));

irf_timeaxis(h(1:2),'nodate');

 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters whistlerfac.eps;
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

plot(h(1),Bfacs(:,3),Bfacs(:,2),Bfacs(1,3),Bfacs(1,2),'ko',Bfacs(end,3),Bfacs(end,2),'kx');
ylabel(h(1),'B_x (nT)');
xlabel(h(1),'B_y (nT)');

 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters whistlerhodogram.eps;
end