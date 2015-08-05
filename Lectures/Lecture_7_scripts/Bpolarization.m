ic = 4; %Spacecraft number 1--4

SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');
Bisr2=c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');

Bmag = sqrt(Bisr2.data(:,1).^2+Bisr2.data(:,2).^2+Bisr2.data(:,3).^2);
LHfreq = (1.6e-19)*Bmag*1e-9/sqrt(9.1e-31*1.673e-27)/(2*pi);
ecfreq = (1.6e-19)*Bmag*1e-9/(9.1e-31*2*pi);

ecfreqs = [ecfreq*0.1, ecfreq*0.5, ecfreq];
ecfreqs = TSeries(Bisr2.time,ecfreqs);


Bibm=c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','ts');
Eibm=c_caa_var_get(irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L2_EB',ic),'caa','ts');

tlimit = irf.tint(Eibm.time.start.utc,Eibm.time.stop.utc);


% Estimate Ez
[Eibm,angle] = irf_edb(Eibm,Bisr2,10,'E.B=0');

%Eibm(isnan(Eibm)) = 0;
Bibm.data(:,1) = Bibm.data(:,1)-nanmean(Bibm.data(:,1));
Bibm.data(:,2) = Bibm.data(:,2)-nanmean(Bibm.data(:,2));
Bibm.data(:,3) = Bibm.data(:,3)-nanmean(Bibm.data(:,3));
Bibm.data(isnan(Bibm.data)) = 0;

Bfac = irf_convert_fac(Bibm,Bisr2,SCpos);

%[~,idx]=irf_tlim(Bibm,[irf_time([2003 03 01 03 59 07.4]) irf_time([2003 03 01 03 59 08.2])]);
%Bibms = Bibm(idx,:);

polarization = irf_ebsp(Eibm,Bibm,Bisr2,Bisr2,SCpos,[10 2200],'polarization','fac');

frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Esum2D = polarization.ee_ss;
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;

Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;


h=irf_plot(6,'newfigure'); 

h(1)=irf_panel('Bfac');
irf_plot(h(1),Bfac);
ylabel(h(1),'B_{FAC} (nT)');
irf_legend(h(1),{'B_x','B_y','B_z'},[0.5 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Bsum');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=Bsum;
    specrec.f_label='';
    specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
    irf_spectrogram(h(2),specrec,'log','donotfitcolorbarlabel');
  irf_legend(h(2),'(b)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(2),'on');
irf_plot(h(2),ecfreqs,'linewidth',1.5,'color','w')
hold(h(2),'off');
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3]);
caxis(h(2),[-8 -2])
ylabel(h(2),'f (Hz)','fontsize',12);

h(3)=irf_panel('ellipt');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=ellipticity;
    specrec.f_label='';
    specrec.p_label={'Ellipticity'};
    irf_spectrogram(h(3),specrec,'lin','donotfitcolorbarlabel');
  irf_legend(h(3),'(c)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(3),'on');
irf_plot(h(3),ecfreqs,'linewidth',1.5,'color','w')
hold(h(3),'off');
set(h(3),'yscale','log');
set(h(3),'ytick',[1e1 1e2 1e3]);
caxis(h(3),[-1, 1])
ylabel(h(3),'f (Hz)','fontsize',12);

h(4)=irf_panel('thetak');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=thetak;
    specrec.f_label='';
    specrec.p_label={'\theta_{k}'};
    irf_spectrogram(h(4),specrec,'lin','donotfitcolorbarlabel');
  irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(4),'on');
irf_plot(h(4),ecfreqs,'linewidth',1.5,'color','w')
hold(h(4),'off');
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3]);
caxis(h(4),[0, 90])
ylabel(h(4),'f (Hz)','fontsize',12);

h(5)=irf_panel('dop');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=dop;
    specrec.f_label='';
    specrec.p_label={'DOP'};
    irf_spectrogram(h(5),specrec,'lin','donotfitcolorbarlabel');
  irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(5),'on');
irf_plot(h(5),ecfreqs,'linewidth',1.5,'color','w')
hold(h(5),'off');
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3]);
caxis(h(5),[0, 1])
ylabel(h(5),'f (Hz)','fontsize',12);

h(6)=irf_panel('planarity');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=planarity;
    specrec.f_label='';
    specrec.p_label={'planarity'};
    irf_spectrogram(h(6),specrec,'lin','donotfitcolorbarlabel');
  irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',12)
  hold(h(6),'on');
irf_plot(h(6),ecfreqs,'linewidth',1.5,'color','w')
hold(h(6),'off');
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3]);
caxis(h(6),[0, 1])
ylabel(h(6),'f (Hz)','fontsize',12);

xwidth = 0.73;
ywidth = 0.12;
%set(h(1),'position',[0.10 0.99-ywidth xwidth ywidth]);
%set(h(2),'position',[0.10 0.99-2*ywidth xwidth ywidth]);
%set(h(3),'position',[0.10 0.99-3*ywidth xwidth ywidth]);
  
set(h(2:6),'xgrid','off','ygrid','off')
set(h(2:6),'Color',0.7*[1 1 1]);

irf_zoom(h(1:6),'x',tlimit);
set(h(1:6),'fontsize',12)

irf_plot_axis_align(h(1:6))

set(gcf, 'InvertHardCopy', 'off');
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print('-dpng','-painters','-r600','polarizationeg.png');
