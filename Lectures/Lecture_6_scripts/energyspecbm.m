% Plot electron and ion energy spectrograms using burst mode data
ic = 1;
  
h=irf_plot(6,'newfigure'); 


BGSE = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'caa','ts');

hca=irf_panel('BGSE');
irf_plot(hca,BGSE);
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.9 0.1]);
irf_legend(hca,'(a)',[0.99 0.98],'color','k')
irf_zoom(hca,'y',[-80 80]);
ylabel(hca,'B_{GSE}','interpreter','tex');

VGSE = c_caa_var_get(irf_ssub('velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic),'caa','ts');

hca=irf_panel('VGSE');
irf_plot(hca,VGSE);
irf_legend(hca,{'V_{x}','V_{y}','V_{z}'},[0.9 0.1]);
irf_legend(hca,'(b)',[0.99 0.98],'color','k')
ylabel(hca,'V_{GSE}','interpreter','tex');


%Get electron differential energy flux data at subspin resolution
edefluxL=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DXL_DEFlux',ic));
edeflux2omni = edefluxL.omni;
edefluxt = edefluxL.tt;
edefluxe = edefluxL.en;
thdeflux = edefluxL.theta;

hca=irf_panel('PEACEomni');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edefluxL.dataunits ')']);
    specrec.f=edefluxe;
    specrec.p=edeflux2omni;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(c)',[0.99 0.98],'color','k')
  caxis(hca,[5 9]);
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[20,max(edefluxe)])
  ylabel(hca,{'E_{e} (eV)'});
  
  
E1 = 200;
[~,idx] = min(abs(E1-edefluxe));
epitch1 = edefluxL.data(:,:,idx);
E1 = edefluxe(idx);

hca=irf_panel('PEACEpitch1');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edefluxL.dataunits ')']);
    specrec.f=thdeflux;
    specrec.p=epitch1;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(d)',[0.99 0.98],'color','k')
  caxis(hca,[5 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E1(1)),'eV')});
  
%subspin resolution ions
idefluxX3D=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_HIA_HS_MAG_IONS_PEF',ic));
ideflux2omni = idefluxX3D.omni;
idefluxt = idefluxX3D.tt;
idefluxe = idefluxX3D.en;
ithdeflux = idefluxX3D.theta;

hca=irf_panel('CISHAIomni');
  specrec=struct('t',idefluxt,'dt',idefluxt(2)-idefluxt(1),'p_label',['(' idefluxX3D.dataunits ')']);
    specrec.f=idefluxe;
    specrec.p=ideflux2omni;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(e)',[0.99 0.98],'color','k')
  caxis(hca,[4 8]);
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[20,max(idefluxe)])
  ylabel(hca,{'E_{i} (eV)'});

E1 = 100;
[~,idx] = min(abs(E1-idefluxe));
ipitch1 = idefluxX3D.data(:,:,idx);
E1 = idefluxe(idx);

hca=irf_panel('CISHIApitch1');
  specrec=struct('t',idefluxt,'dt',idefluxt(2)-idefluxt(1),'p_label',['(' idefluxX3D.dataunits ')']);
    specrec.f=ithdeflux;
    specrec.p=ipitch1;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(f)',[0.99 0.98],'color','k')
  caxis(hca,[4 8]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E1(1)),'eV')});

%select short time interval  
tint = [iso2epoch('2008-04-22T18:00:00.00Z') iso2epoch('2008-04-22T18:10:00.00Z')]
irf_zoom(h,'x',tint);
set(h,'fontsize',12)


irf_plot_axis_align(h)

if 1,
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print('-dpng','-painters','-r600','burstmodeexample.png');
end