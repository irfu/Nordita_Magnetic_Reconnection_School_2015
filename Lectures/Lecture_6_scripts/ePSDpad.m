%Plot phase-space electron densities, at pitch angle 0, 90, and 180 deg.
%Spin resolution data.

ic = 1;

tint = iso2epoch('2008-04-22T18:06:55.0Z');

ePSD = c_caa_var_get(irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSD = ePSD.data;
ePSDt = c_caa_var_get(irf_ssub('time_tags__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSDt = ePSDt.data;
ePSDtmin = c_caa_var_get(irf_ssub('time_tags_DeltaLower__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSDtmin = double(ePSDtmin.data);
ePSDtmax = c_caa_var_get(irf_ssub('time_tags_DeltaUpper__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSDtmax = double(ePSDtmax.data);
theta = c_caa_var_get(irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
theta = theta.data';
ePSDbg = c_caa_var_get(irf_ssub('BackgroundLevel__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
ePSDbg = squeeze(nanmean(ePSDbg.data,2));

n=c_caa_var_get(irf_ssub('Spacecraft_potential__C?_CP_EFW_L2_P',ic),'caa','ts');

%Find closed time to tint
[~,idx] = min(abs(tint-ePSDt));
energy = c_caa_var_get(irf_ssub('Sweep_Energy__C?_CP_PEA_PITCH_SPIN_PSD',ic),'caa');
energy = energy.data(idx,:);

trange = [ePSDt(idx)-ePSDtmin(idx) ePSDt(idx)+ePSDtmax(idx)];
fprintf(strcat('Electron distribution obtained over interval: \n',epoch2iso(trange(1)),'/',epoch2iso(trange(2)),'\n'));

dst = squeeze(ePSD(idx,:,:));
bg = ePSDbg(idx,:);

n = n.tlim(irf_time(trange,'epoch>epochTT'));
offset1 = 1.2*nanmean(n.data);
if isnan(offset1),
    offset1 = 0;
end
energy = energy+offset1;

energy(find(energy < 1)) = NaN;

h=irf_plot(1,'newfigure');
h(1)=irf_panel('EdistE1');
plot(h(1),energy,dst(1,:),'ko',energy,mean(dst([6 7],:)),'ro',energy,dst(12,:),'bo',energy,bg,'k--');
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
 print -depsc -painters ePSDpadteststandard.eps;
end