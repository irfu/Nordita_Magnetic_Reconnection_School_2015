ic = 2; %set spacecraft number 1--4.

Bisr2 = c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');
Eibm = c_caa_var_get(irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L2_EB',ic),'caa','mat');
Bibm = c_caa_var_get(irf_ssub('B_Vec_xyz_ISR2__C?_CP_EFW_L2_BB',ic),'caa','mat');
SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','mat');

tint = [Eibm(1,1) Eibm(end,1)];


%coordinate transformations 
[Epar,Eperp]=irf_dec_parperp(Bisr2,Eibm,1);
Eparperp = [Epar(:,1) Epar(:,2) Eperp(:,2)];

Bfac = irf_convert_fac(Bibm,Bisr2,SCpos);

%calculate wavelet transforms
Bibx = [Bibm(:,1), Bibm(:,2)];
Bxwavelet = irf_wavelet(Bibx);
Biby = [Bibm(:,1), Bibm(:,3)];
Bywavelet = irf_wavelet(Biby);
Bibz = [Bibm(:,1), Bibm(:,4)];
Bzwavelet = irf_wavelet(Bibz);

Eibx = [Eibm(:,1), Eibm(:,2)];
Exwavelet = irf_wavelet(Eibx);
Eiby = [Eibm(:,1), Eibm(:,3)];
Eywavelet = irf_wavelet(Eiby);

Bmag = sqrt(Bisr2(:,2).^2+Bisr2(:,3).^2+Bisr2(:,4).^2);
LHfreq = (1.6e-19)*Bmag*1e-9/sqrt(9.1e-31*1.673e-27)/(2*pi);
ecfreq = (1.6e-19)*Bmag*1e-9/(9.1e-31*2*pi);


%plot internal burst mode electric and magnetic fields
irf_plot(7)

h(1)=irf_panel('BISR2');
irf_plot(h(1),Bisr2);
ylabel(h(1),'B_{ISR2} (nT)');
irf_zoom(h(1),'y',[-100,100])
irf_legend(h(1),{'B_x','B_y','B_z'},[0.98 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Eibm');
irf_plot(h(2),Eibm);
ylabel(h(2),'E_{ISR2} (mV/m)');
irf_legend(h(2),{'E_x','E_y'},[0.98 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('Bibm');
irf_plot(h(3),Bibm);
ylabel(h(3),'B_{ISR2} (nT)');
irf_legend(h(3),{'B_x','B_y','B_z'},[0.98 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('Eparperp');
irf_plot(h(4),Eparperp);
ylabel(h(4),'E (mV/m)');
irf_legend(h(4),{'E_{||}','E_{\perp}'},[0.98 0.1])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

h(5)=irf_panel('Bfac');
irf_plot(h(5),Bfac);
ylabel(h(5),'B_{FAC} (nT)');
irf_legend(h(5),{'B_x','B_y','B_z'},[0.98 0.1])
irf_legend(h(5),'(e)',[0.99 0.98],'color','k')

h(6)=irf_panel('Bwavelet');
  specrec=struct('t',Bxwavelet.t,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=Bxwavelet.f;
    specrec.p=cell2mat(Bxwavelet.p)+cell2mat(Bywavelet.p)+cell2mat(Bzwavelet.p);
    specrec.f_label='';
    specrec.p_label={'log_{10} B^2','nT^{2} Hz^{-1}'};
    irf_spectrogram(h(6),specrec);
  irf_legend(h(6),'(f)',[0.99 0.98],'color','w')
  hold(h(6),'on');
irf_plot(h(6),[Bisr2(:,1),ecfreq],'linewidth',1.5,'color','w')
irf_plot(h(6),[Bisr2(:,1),ecfreq*0.5],'linewidth',1.5,'color','w')
irf_plot(h(6),[Bisr2(:,1),ecfreq*0.1],'linewidth',1.5,'color','w')
irf_plot(h(6),[Bisr2(:,1),LHfreq],'linewidth',1.5,'color','k')
hold(h(6),'off');
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3]);
caxis(h(6),[-9, -3])
ylabel(h(6),'f (Hz)');

h(7)=irf_panel('Ewavelet');
  specrec=struct('t',Exwavelet.t,'p_label',['mV^{2}m^{-2}Hz^{-1}']);
    specrec.f=Exwavelet.f;
    specrec.p=cell2mat(Exwavelet.p)+cell2mat(Eywavelet.p);
    specrec.f_label='';
    specrec.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(h(7),specrec,'donotfitcolorbarlabel');
  irf_legend(h(7),'(g)',[0.99 0.98],'color','w')
  hold(h(7),'on');
irf_plot(h(7),[Bisr2(:,1),ecfreq],'linewidth',1.5,'color','w')
irf_plot(h(7),[Bisr2(:,1),ecfreq*0.5],'linewidth',1.5,'color','w')
irf_plot(h(7),[Bisr2(:,1),ecfreq*0.1],'linewidth',1.5,'color','w')
irf_plot(h(7),[Bisr2(:,1),LHfreq],'linewidth',1.5,'color','k')
hold(h(7),'off');
set(h(7),'yscale','log');
set(h(7),'ytick',[1e1 1e2 1e3]);
caxis(h(7),[-6, 0])
ylabel(h(7),'f (Hz)');

irf_plot_axis_align(7)
irf_zoom(h,'x',tint);
irf_timeaxis(h);