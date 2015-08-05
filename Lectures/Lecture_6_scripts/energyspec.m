% Plot electron and ion energy spectrograms using standard data
ic = 1;
  
h=irf_plot(4,'newfigure'); 


BGSE = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'caa','ts');

hca=irf_panel('BGSE');
irf_plot(hca,BGSE);
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.9 0.05]);
irf_legend(hca,'(a)',[0.99 0.98],'color','k')
irf_zoom(hca,'y',[-80 80]);
ylabel(hca,'B_{GSE}');

VGSE = c_caa_var_get(irf_ssub('velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic),'caa','ts');

hca=irf_panel('VGSE');
irf_plot(hca,VGSE);
irf_legend(hca,{'V_{x}','V_{y}','V_{z}'},[0.9 0.05]);
irf_legend(hca,'(b)',[0.99 0.98],'color','k')
ylabel(hca,'V_{GSE}');


%Get electron differential energy flux data
edeflux = c_caa_var_get(irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edeflux2 = edeflux.data;
edeflux2(edeflux2 < 0)=NaN;
edeflux2omni = squeeze(nanmean(edeflux2,2));
edefluxt = c_caa_var_get(irf_ssub('time_tags__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxt = edefluxt.data;
edefluxe = c_caa_var_get(irf_ssub('Sweep_Energy__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxe = edefluxe.data(1,:);


hca=irf_panel('PEACEomni');
  specrec=struct('t',edefluxt,'dt',edefluxt(2)-edefluxt(1),'p_label',['(' edeflux.UNITS ')']);
    specrec.f=edefluxe;
    specrec.p=edeflux2omni;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(c)',[0.99 0.98],'color','w')
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
hold(h(3),'on');
irf_plot(h(3),Te);
hold(h(3),'off');
end  
  
  
ideflux = c_caa_var_get(irf_ssub('flux__C?_CP_CIS_HIA_HS_1D_PEF',ic),'caa');  
ideflux2 = ideflux.data;
ideflux2(ideflux2 < 0)=NaN;
idefluxt = c_caa_var_get(irf_ssub('time_tags__C?_CP_CIS_HIA_HS_1D_PEF',ic),'caa');
idefluxt = idefluxt.data;
idefluxe = c_caa_var_get(irf_ssub('energy_table__C?_CP_CIS_HIA_HS_1D_PEF',ic),'caa');
idefluxe = idefluxe.data;

hca=irf_panel('CISHIAomni');
  specrec=struct('t',idefluxt,'dt',idefluxt(2)-idefluxt(1),'p_label',['(' ideflux.UNITS ')']);
    specrec.f=idefluxe;
    specrec.p=ideflux2;
    specrec.f_label='';
    specrec.p_label={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(d)',[0.99 0.98],'color','w')
  caxis(hca,[4 8]);
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  set(hca,'xgrid','off')
  irf_zoom(hca,'y',[min(idefluxe),max(idefluxe)])
  ylabel(hca,{'E_{i} (eV)'});
  
% Plot ion temperature over spectrogram
if 0,
Ti = c_caa_var_get(irf_ssub('temperature__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic),'caa','mat');
Ti(:,2) = Ti(:,2)*86.132; %convert from MK to eV
hold(h(4),'on');
irf_plot(h(4),Ti);
hold(h(4),'off');
end

tint = [edefluxt(1) edefluxt(end)];  
irf_zoom(h,'x',tint);
set(h,'fontsize',12)

irf_plot_axis_align(h)

if 0, % set to 1 to save plot
set(gcf, 'InvertHardCopy', 'off');

set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print('-dpng','-painters','-r600','electroniondefluxes.png');
end
