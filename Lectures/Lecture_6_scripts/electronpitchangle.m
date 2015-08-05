% Plot electron and ion energy spectrograms using standard data
ic = 1;
  
h=irf_plot(5,'newfigure'); 

BGSE = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'caa','ts');

hca=irf_panel('BGSE');
irf_plot(hca,BGSE);
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.9 0.05]);
irf_legend(hca,'(a)',[0.99 0.98],'color','k')
irf_zoom(hca,'y',[-80 80]);
ylabel(hca,'B_{GSE}','interpreter','tex');

%Get electron differential energy flux data
edeflux = c_caa_var_get(irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edeflux2 = edeflux.data;
edeflux2(edeflux2 < 0)=NaN;
edeflux2omni = squeeze(nanmean(edeflux2,2));
edefluxt = c_caa_var_get(irf_ssub('time_tags__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxt = edefluxt.data;
edefluxe = c_caa_var_get(irf_ssub('Sweep_Energy__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxe = edefluxe.data(1,:);
thdeflux = c_caa_var_get(irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
thdeflux = thdeflux.data;


% select energy value
E1 = 2000;
[~,idx] = min(abs(E1-edefluxe));
epitch1 = edeflux2(:,:,idx);
E1 = edefluxe(idx);

E2 = 200;
[~,idx] = min(abs(E2-edefluxe));
epitch2 = edeflux2(:,:,idx);
E2 = edefluxe(idx);

E3 = 50;
[~,idx] = min(abs(E3-edefluxe));
epitch3 = edeflux2(:,:,idx);
E3 = edefluxe(idx);

tint = [edefluxt(1) edefluxt(end)]; 


hca=irf_panel('PEACEomni');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edeflux.UNITS ')']);
    specrec.f=edefluxe;
    specrec.p=edeflux2omni;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(a)',[0.99 0.98],'color','w')
  hold(hca,'on')
  irf_plot(hca,[tint(1)-1000 tint(2)+1000; E1 E1]');
  irf_plot(hca,[tint(1)-1000 tint(2)+1000; E2 E2]');
  irf_plot(hca,[tint(1)-1000 tint(2)+1000; E3 E3]');
  hold(hca,'off')
  caxis(hca,[5 9]);
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[10,max(edefluxe)])
  ylabel(hca,{'E_{e} (eV)'});

% Plot electron temperature over spectrogram
if 0,
Tepar = c_caa_var_get(irf_ssub('Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS',ic),'caa','mat');
Teperp = c_caa_var_get(irf_ssub('Data_Temperature_ComponentPerpendicularToMagField__C?__MOMENTS',ic),'caa','mat');
Te = [Tepar(:,1) (Tepar(:,2)+2*Teperp(:,2))/3*86.132]; %calculate total temperature and convert MK to eV
hold(hca,'on');
irf_plot(hca,Te);
hold(hca,'off');
end  

hca=irf_panel('PEACEpitch1');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edeflux.UNITS ')']);
    specrec.f=thdeflux;
    specrec.p=epitch1;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(b)',[0.99 0.98],'color','w')
  caxis(hca,[5 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E1(1)),'eV')});

hca=irf_panel('PEACEpitch2');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edeflux.UNITS ')']);
    specrec.f=thdeflux;
    specrec.p=epitch2;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(c)',[0.99 0.98],'color','w')
  caxis(hca,[5 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E1(1)),'eV')});

hca=irf_panel('PEACEpitch3');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edeflux.UNITS ')']);
    specrec.f=thdeflux;
    specrec.p=epitch3;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log');
  irf_legend(hca,'(d)',[0.99 0.98],'color','w')
  caxis(hca,[5 9]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[0, 180])
  ylabel(hca,{'\theta (deg)',strcat('E=',num2str(E1(1)),'eV')});
 
irf_zoom(h,'x',tint);
set(h,'fontsize',12)

irf_plot_axis_align(h)

if 0, % set to 1 to save plot
set(gcf, 'InvertHardCopy', 'off');

set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print('-dpng','-painters','-r600','edefluxpitchangle.png');
end