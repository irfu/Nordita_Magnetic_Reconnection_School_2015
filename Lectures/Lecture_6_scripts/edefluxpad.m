%Plot pitch angle distribution of differential energy flux

ic = 1;

tint = iso2epoch('2008-04-22T18:05:00.0Z');

edeflux = c_caa_var_get(irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edeflux = edeflux.data;
edefluxt = c_caa_var_get(irf_ssub('time_tags__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxt = edefluxt.data;
edefluxtmin = c_caa_var_get(irf_ssub('time_tags_DeltaLower__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxtmin = double(edefluxtmin.data);
edefluxtmax = c_caa_var_get(irf_ssub('time_tags_DeltaUpper__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
edefluxtmax = double(edefluxtmax.data);
theta = c_caa_var_get(irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
theta = theta.data';

n=c_caa_var_get(irf_ssub('Spacecraft_potential__C?_CP_EFW_L2_P',ic),'caa','ts');

[~,idx] = min(abs(tint-edefluxt));
energy = c_caa_var_get(irf_ssub('Sweep_Energy__C?_CP_PEA_PITCH_SPIN_DEFlux',ic),'caa');
energy = energy.data(idx,:);

trange = [edefluxt(idx)-edefluxtmin(idx) edefluxt(idx)+edefluxtmax(idx)];
fprintf(strcat('Electron distribution obtained over interval: \n',epoch2iso(trange(1)),'/',epoch2iso(trange(2)),'\n'));


dsf = squeeze(edeflux(idx,:,:));

%LEEA data
n = n.tlim(irf_time(trange,'epoch>epochTT'));
offset1 = 1.2*nanmean(n.data);
if isnan(offset1),
    offset1 = 0;
end
energy = energy+offset1;

energy(find(energy < 1)) = NaN;

energy = log10(energy);
theta2 = [theta, 180+theta, theta(1)]+90-theta(1);
dsf = [dsf; flipdim(dsf,1); dsf(1,:)];

[TH,R] = meshgrid(theta2,energy);
[X,Y] = pol2cart(TH*pi/180.0,R);


h=irf_plot(1,'newfigure');
% Plot DIFFERENTIAL ENERGY FLUX 2D plot
h(1)=irf_panel('defpad');
pcolor(h(1),double(X),double(Y),real(double(log10(dsf'))))
hold(h(1),'on')
plot(h(1),[0 0],[-5 5],'k');
plot(h(1),[-5 5],[0 0],'k');
hold(h(1),'off')
irf_legend(h(1),{'\theta = 0^{o}'},[0.65 0.98],'color','k','fontsize',20)
irf_legend(h(1),{'\theta = 180^{o}'},[0.70 0.02],'color','k','fontsize',20)
%view(h(3),2); 
%caxis(h(3),[2, 10])
axis(h(1),'equal','tight'); shading(h(1),'flat'); grid(h(1),'off');
ylabel(h(1),'log_{10}E [eV]','fontsize',20);
xlabel(h(1),'log_{10}E [eV]','fontsize',20);
axis(h(1),[-max(energy) max(energy) -max(energy) max(energy)]);
caxis(h(1),[5 8]);
c=colorbar;
ylabel(c,'log_{10} DEFlux','fontsize',20);
set(c,'fontsize',20)
set(h(1),'position',[0.09 0.2 0.76 0.78]);
set(h(1),'fontsize',20)

if 0, % set to 1 to save plot
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters edefluxpadteststandard.eps;
end