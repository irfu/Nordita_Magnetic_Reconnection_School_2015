% Plot electron and ion energy spectrograms using standard data
ic = 1;
  
h=irf_plot(5,'newfigure'); 

BGSE = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'caa','mat');

hca=irf_panel('BGSE');
irf_plot(hca,BGSE);
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.9 0.05]);
irf_legend(hca,'(a)',[0.99 0.98],'color','k')
irf_zoom(hca,'y',[-80 80]);
ylabel(hca,'B_{GSE}');

%Get electron differential energy flux data
idpflux = c_caa_var_get(irf_ssub('Differential_Particle_Flux__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
idpflux2 = idpflux.data;
idpflux2(idpflux2 < 0)=NaN;
idpflux2omni = squeeze(nanmean(idpflux2,2));
idpfluxt = c_caa_var_get(irf_ssub('time_tags__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
idpfluxt = idpfluxt.data;
idpfluxe = c_caa_var_get(irf_ssub('energy_table__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
idpfluxe = idpfluxe.data;
thdpflux = c_caa_var_get(irf_ssub('pitch_angle__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
thdpflux = thdpflux.data;


% select energy value
E1 = 5000;
[~,idx] = min(abs(E1-idpfluxe));
ipitch1 = idpflux2(:,:,idx);
E1 = idpfluxe(idx);

E2 = 800;
[~,idx] = min(abs(E2-idpfluxe));
ipitch2 = idpflux2(:,:,idx);
E2 = idpfluxe(idx);

E3 = 50;
[~,idx] = min(abs(E3-idpfluxe));
ipitch3 = idpflux2(:,:,idx);
E3 = idpfluxe(idx);

tint = [idpfluxt(1) idpfluxt(end)]; 


hca=irf_panel('CISomni');
  specrec=struct('t',idpfluxt,'dt',idpfluxt(2)-idpfluxt(1),'p_label',['(' idpflux.UNITS ')']);
    specrec.f=idpfluxe;
    specrec.p=idpflux2omni;
    specrec.f_label='';
    specrec.p_label={'log_{10} dPF','Particles/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(b)',[0.99 0.98],'color','w')
  hold(hca,'on')
  irf_plot(hca,[tint(1)-1000 tint(2)+1000; E1 E1]');
  irf_plot(hca,[tint(1)-1000 tint(2)+1000; E2 E2]');
  irf_plot(hca,[tint(1)-1000 tint(2)+1000; E3 E3]');
  hold(hca,'off')
  caxis(hca,[4 9]);
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[10,max(idpfluxe)])
  ylabel(hca,{'E_{i} (eV)'});


hca=irf_panel('CISpitch1');
  specrec=struct('t',idpfluxt,'dt',idpfluxt(2)-idpfluxt(1),'p_label',['(' idpflux.UNITS ')']);
    specrec.f=thdpflux;
    specrec.p=ipitch1;
    specrec.f_label='';
    specrec.p_label={'log_{10} dPF','Particles/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(c)',[0.99 0.98],'color','w')
  caxis(hca,[4 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E1(1)),'eV')});

hca=irf_panel('CISpitch2');
  specrec=struct('t',idpfluxt,'dt',idpfluxt(2)-idpfluxt(1),'p_label',['(' idpflux.UNITS ')']);
    specrec.f=thdpflux;
    specrec.p=ipitch2;
    specrec.f_label='';
    specrec.p_label={'log_{10} dPF','Particles/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(d)',[0.99 0.98],'color','w')
  caxis(hca,[4 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E2(1)),'eV')});

hca=irf_panel('CISpitch3');
  specrec=struct('t',idpfluxt,'dt',idpfluxt(2)-idpfluxt(1),'p_label',['(' idpflux.UNITS ')']);
    specrec.f=thdpflux;
    specrec.p=ipitch3;
    specrec.f_label='';
    specrec.p_label={'log_{10} dPF','Particles/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(e)',[0.99 0.98],'color','w')
  caxis(hca,[4 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E3(1)),'eV')});
 
irf_zoom(h,'x',tint);
set(h,'fontsize',12)

irf_plot_axis_align(h)

if 0, % set to 1 to save plot
set(gcf, 'InvertHardCopy', 'off');

set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print('-dpng','-painters','-r600','idefluxpitchangle.png');
end