%Plot pitch angle distribution of differential energy flux. Using 2D pitch
%angle distributions. 

ic = 1;

%select time of interest
tint = [iso2epoch('2008-04-22T18:06:26.0Z')];

idpflux = c_caa_var_get(irf_ssub('Differential_Particle_Flux__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
idpflux = idpflux.data;
idpfluxt = c_caa_var_get(irf_ssub('time_tags__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
idpfluxt = idpfluxt.data;
idpfluxta = c_caa_var_get(irf_ssub('duration__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
idpfluxta = double(idpfluxta.data);
idpfluxtmid = idpfluxt+double(idpfluxta/2.0);
energy = c_caa_var_get(irf_ssub('energy_table__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
energy = energy.data;
theta = c_caa_var_get(irf_ssub('pitch_angle__C?_CP_CIS_HIA_PAD_HS_MAG_IONS_PF',ic),'caa');
theta = theta.data';

[~,idx] = min(abs(tint-idpfluxtmid));
trange = [idpfluxt(idx) idpfluxt(idx)+idpfluxta(idx)];
fprintf(strcat('Ion distribution obtained over interval: \n',epoch2iso(trange(1)),'/',epoch2iso(trange(2)),'\n'));

dsf = squeeze(idpflux(idx,:,:));

energy(find(energy < 1)) = NaN;

%energy = log10(energy);
energy = sqrt(2*1.6e-19*energy/1.67e-27)/1000; %convert to km/s

theta2 = [theta, 180+theta, theta(1)]+90-theta(1);
dsf = [dsf; flipdim(dsf,1); dsf(1,:)];

[TH,R] = meshgrid(theta2,energy);
[X,Y] = pol2cart(TH*pi/180.0,R);


h=irf_plot(1,'newfigure');
% Plot DIFFERENTIAL ENERGY FLUX 2D plot
h(1)=irf_panel('dpfpad');
pcolor(h(1),double(X),double(Y),real(double(log10(dsf'))))
hold(h(1),'on')
plot(h(1),[0 0],[-max(energy) max(energy)],'k');
plot(h(1),[-max(energy) max(energy)],[0 0],'k');
hold(h(1),'off')
irf_legend(h(1),{'\theta = 0^{o}'},[0.65 0.98],'color','k','fontsize',20)
irf_legend(h(1),{'\theta = 180^{o}'},[0.70 0.02],'color','k','fontsize',20)
%view(h(3),2); 
%caxis(h(3),[2, 10])
axis(h(1),'equal','tight'); shading(h(1),'flat'); grid(h(1),'off');
ylabel(h(1),'v_{||} (km s^{-1})','fontsize',20);
xlabel(h(1),'v_{\perp} (km s^{-1})','fontsize',20);
maxv = 800;
axis(h(1),[-maxv maxv -maxv maxv]);
%caxis(h(1),[5 8]);
c=colorbar;
ylabel(c,'log_{10} DPFlux','fontsize',20);
set(c,'fontsize',20)
set(h(1),'position',[0.10 0.2 0.76 0.78]);
set(h(1),'fontsize',20)

if 0, % set to 1 to save plot
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters idefluxpadsheath.eps;
end