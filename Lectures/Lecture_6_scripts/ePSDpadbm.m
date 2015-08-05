%Plot pitch angle distribution of phase-space density using 3DX, 3DR
%data. 

ic = 1;

tint = [iso2epoch('2008-04-22T18:06:59.90Z') iso2epoch('2008-04-22T18:07:02.20Z')];

LEEA3DPSD=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DXL_PSD',ic));
HEEA3DPSD=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DRH_PSD',ic)); %Use 3DXH if available otherwise 3DRH
n=c_caa_var_get(irf_ssub('Spacecraft_potential__C?_CP_EFW_L2_P',ic),'caa','ts');
ePSDbg = c_caa_var_get(irf_ssub('BackgroundLevel__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSDbg = squeeze(nanmean(ePSDbg.data,2));
ePSDtbg = c_caa_var_get(irf_ssub('time_tags__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSDtbg = ePSDtbg.data;

[~,idx] = min(abs(tint(1)-ePSDtbg));
bg = ePSDbg(idx,:);
energybg = c_caa_var_get(irf_ssub('Sweep_Energy__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
energybg = energybg.data(idx,:);
    
energy = LEEA3DPSD.en;
energyH = HEEA3DPSD.en;
theta = LEEA3DPSD.theta;

%LEEA data
n = n.tlim(irf_time(tint,'epoch>epochTT'));
offset1 = 1.2*nanmean(n.data);
if isnan(offset1),
    offset1 = 0;
end
energy = energy+offset1;
energybg = energybg+offset1;
[~,idx]=irf_tlim(LEEA3DPSD.tt,tint);
disttemp = squeeze(LEEA3DPSD.data(idx,:,:));
dst = squeeze(nanmean(disttemp,1));

%HEEA data 
[~,idx]=irf_tlim(HEEA3DPSD.tt,tint);
disttemp = squeeze(HEEA3DPSD.data(idx,:,:));
dsth = squeeze(nanmean(disttemp,1));

%Combine LEEA and HEEA data
elstemp = find(energy < energyH(1));
et = cat(2,energy(elstemp),energyH);
dst = cat(2,dst(:,elstemp),dsth);
et(find(et < 2)) = NaN;

h=irf_plot(1,'newfigure');
h(1)=irf_panel('EdistE1');
plot(h(1),et,dst(1,:),'ko',et,mean(dst([6 7],:)),'ro',et,dst(12,:),'bo',energybg,bg,'k--');
ylabel(h(1),'f_e (s^3 km^{-6})','fontsize',20);
xlabel(h(1),'E (eV)','fontsize',20);
set(h(1),'yscale','log');
set(h(1),'xscale','log');
set(h(1),'xtick',[1 1e1 1e2 1e3 1e4 1e5]);
irf_legend(h(1),{'0 deg'},[0.94 0.94],'color','k','fontsize',20)
irf_legend(h(1),{'90 deg'},[0.94 0.88],'color','r','fontsize',20)
irf_legend(h(1),{'180 deg'},[0.94 0.82],'color','b','fontsize',20)
irf_zoom(h(1),'x',[2 3e4])
%irf_zoom(h(1),'y',[0.00001 1e4])
set(h(1),'position',[0.15 0.15 0.76 0.78]);
set(h(1),'fontsize',20)

if 0, % set to 1 to save plot
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters ePSDpadtestburst.eps;
end