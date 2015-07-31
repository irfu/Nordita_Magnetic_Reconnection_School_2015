ic = 2; %Spacecraft number 1--4

SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');
Bisr2=c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');

Bmag = sqrt(Bisr2(:,2).^2+Bisr2(:,3).^2+Bisr2(:,4).^2);
LHfreq = (1.6e-19)*Bmag*1e-9/sqrt(9.1e-31*1.673e-27)/(2*pi);
ecfreq = (1.6e-19)*Bmag*1e-9/(9.1e-31*2*pi);


Bibm=c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','mat');
Eibm=c_caa_var_get(irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L2_EB',ic),'caa','mat');

tint = [Eibm(1,1) Eibm(end,1)];


% Estimate Ez
[Eibm,angle] = irf_edb(Eibm,Bisr2,10,'E.B=0');

Bibm(:,2) = Bibm(:,2)-nanmean(Bibm(:,2));
Bibm(:,3) = Bibm(:,3)-nanmean(Bibm(:,3));
Bibm(:,4) = Bibm(:,4)-nanmean(Bibm(:,4));
Bibm(isnan(Bibm)) = 0;

Bfac = irf_convert_fac(Bibm,Bisr2,SCpos);


polarization = irf_ebsp(Eibm,Bibm,Bisr2,Bisr2,SCpos,[10 2200],'polarization','fac');

frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Esum = polarization.ee_xxyyzzss(:,:,4);
Esum2D = polarization.ee_ss;
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);


% Calculate phase speed v_ph = E/B.
vph = sqrt(Esum./Bsum)*1e6;

Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
vph(removepts) = NaN;
pfluxz(removepts) = NaN;


h=irf_plot(4,'newfigure'); 

h(1)=irf_panel('Esum');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=Esum2D;
    specrec.f_label='';
    specrec.p_label={'log_{10}E^{2}','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(h(1),specrec,'log','donotfitcolorbarlabel');
  irf_legend(h(1),'(a)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(1),'on');
irf_plot(h(1),[Bisr2(:,1),ecfreq],'linewidth',1.5,'color','w')
irf_plot(h(1),[Bisr2(:,1),ecfreq/2],'linewidth',1.5,'color','w')
irf_plot(h(1),[Bisr2(:,1),ecfreq/10],'linewidth',1.5,'color','w')
hold(h(1),'off');
set(h(1),'yscale','log');
set(h(1),'ytick',[1e1 1e2 1e3]);
caxis(h(1),[-6 -1])
ylabel(h(1),'f (Hz)','fontsize',12);

h(2)=irf_panel('Bsum');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=Bsum;
    specrec.f_label='';
    specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
    irf_spectrogram(h(2),specrec,'log','donotfitcolorbarlabel');
  irf_legend(h(2),'(b)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(2),'on');
irf_plot(h(2),[Bisr2(:,1),ecfreq],'linewidth',1.5,'color','w')
irf_plot(h(2),[Bisr2(:,1),ecfreq/2],'linewidth',1.5,'color','w')
irf_plot(h(2),[Bisr2(:,1),ecfreq/10],'linewidth',1.5,'color','w')
hold(h(2),'off');
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3]);
caxis(h(2),[-8 -2])
ylabel(h(2),'f (Hz)','fontsize',12);

h(3)=irf_panel('vph');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=vph;
    specrec.f_label='';
    specrec.p_label={'log_{10}E/B','m s^{-1}'};
    irf_spectrogram(h(3),specrec,'log','donotfitcolorbarlabel');
  irf_legend(h(3),'(c)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(3),'on');
irf_plot(h(3),[Bisr2(:,1),ecfreq],'linewidth',1.5,'color','w')
irf_plot(h(3),[Bisr2(:,1),ecfreq/2],'linewidth',1.5,'color','w')
irf_plot(h(3),[Bisr2(:,1),ecfreq/10],'linewidth',1.5,'color','w')
hold(h(3),'off');
set(h(3),'yscale','log');
set(h(3),'ytick',[1e1 1e2 1e3]);
caxis(h(3),[5 8])
ylabel(h(3),'f (Hz)','fontsize',12);

h(4)=irf_panel('thetak');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=pfluxz;
    specrec.f_label='';
    specrec.p_label={'S_{z}/|S|'};
    irf_spectrogram(h(4),specrec,'lin','donotfitcolorbarlabel');
  irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(4),'on');
irf_plot(h(4),[Bisr2(:,1),ecfreq],'linewidth',1.5,'color','w')
irf_plot(h(4),[Bisr2(:,1),ecfreq/2],'linewidth',1.5,'color','w')
irf_plot(h(4),[Bisr2(:,1),ecfreq/10],'linewidth',1.5,'color','w')
hold(h(4),'off');
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3]);
caxis(h(4),[-1 1])
ylabel(h(4),'f (Hz)','fontsize',12);

  
set(h(1:4),'xgrid','off','ygrid','off')
set(h(1:4),'Color',0.7*[1 1 1]);

irf_zoom(h(1:4),'x',tint);
set(h(1:4),'fontsize',12)

irf_plot_axis_align(h(1:4))

%save plot as .png
%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','polarizationeg.png');
