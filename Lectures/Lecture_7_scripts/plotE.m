% Plot internal burst mode E fields and the associated wavelet
% spectrograms.

ic = 2; %set spacecraft number 1--4.

Bisr2 = c_caa_var_get(irf_ssub('B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');
Eibm = c_caa_var_get(irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L2_EB',ic),'caa','ts');
SCpos = c_caa_var_get(irf_ssub('sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2',ic),'caa','ts');

tint = irf.tint(Eibm.time.start.utc,Eibm.time.stop.utc);


%coordinate transformations 
[Epar,Eperp]=irf_dec_parperp(Bisr2,Eibm,1);
Eparperp = [Epar(:,1) Epar(:,2) Eperp(:,2)];


%calculate wavelet transforms
Ewavelet = irf_wavelet(Eparperp,'nf',50);

Bmag = sqrt(Bisr2.data(:,1).^2+Bisr2.data(:,2).^2+Bisr2.data(:,3).^2);
LHfreq = (1.6e-19)*Bmag*1e-9/sqrt(9.1e-31*1.673e-27)/(2*pi);
ecfreq = (1.6e-19)*Bmag*1e-9/(9.1e-31*2*pi);

ecfreqs = [ecfreq*0.1, ecfreq*0.5, ecfreq];
ecfreqs = TSeries(Bisr2.time,ecfreqs);
LHfreq = TSeries(Bisr2.time,LHfreq);

%plot internal burst mode electric fields
irf_plot(5,'newfigure')

h(1)=irf_panel('BISR2');
irf_plot(h(1),Bisr2);
ylabel(h(1),'B_{ISR2} (nT)','Interpreter','tex');
irf_zoom(h(1),'y',[-100,100])
irf_legend(h(1),{'B_x','B_y','B_z'},[0.98 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Eibm');
irf_plot(h(2),Eibm);
ylabel(h(2),'E_{ISR2} (mV/m)','Interpreter','tex');
irf_legend(h(2),{'E_x','E_y'},[0.98 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('Eparperp');
irf_plot(h(3),Eparperp);
ylabel(h(3),'E (mV/m)');
irf_legend(h(3),{'E_{||}','E_{\perp}'},[0.98 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('Ewaveletpar');
  specrec=struct('t',Ewavelet.t);
    specrec.f=Ewavelet.f;
    specrec.p=Ewavelet.p{1,1};
    specrec.f_label='';
    specrec.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(h(4),specrec,'donotfitcolorbarlabel');
  irf_legend(h(4),'(d)',[0.99 0.98],'color','w')
  hold(h(4),'on');
irf_plot(h(4),ecfreqs,'linewidth',1.5,'color','w')
irf_plot(h(4),LHfreq,'linewidth',1.5,'color','k')
hold(h(4),'off');
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3]);
caxis(h(4),[-6, 0])
ylabel(h(4),'f (Hz)');

h(5)=irf_panel('Ewaveletperp');
  specrec=struct('t',Ewavelet.t);
    specrec.f=Ewavelet.f;
    specrec.p=Ewavelet.p{1,2};
    specrec.f_label='';
    specrec.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(h(5),specrec,'donotfitcolorbarlabel');
  irf_legend(h(5),'(e)',[0.99 0.98],'color','w')
  hold(h(5),'on');
irf_plot(h(5),ecfreqs,'linewidth',1.5,'color','w')
irf_plot(h(5),LHfreq,'linewidth',1.5,'color','k')
hold(h(5),'off');
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3]);
caxis(h(5),[-6, 0])
ylabel(h(5),'f (Hz)');

irf_plot_axis_align(h(1:5))
irf_zoom(h(1:5),'x',tint); 
irf_timeaxis(h(1:5));