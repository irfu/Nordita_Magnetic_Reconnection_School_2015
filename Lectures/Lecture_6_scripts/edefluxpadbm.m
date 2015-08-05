%Plot pitch angle distribution of differential energy flux using 3DX, 3DR
%data. 

ic = 1;

tint = [iso2epoch('2008-04-22T18:01:56.271000Z') iso2epoch('2008-04-22T18:02:00.451000Z')];

LEEA3Ddef=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DXL_DEFlux',ic));
HEEA3Ddef=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DRH_DEFlux',ic)); %Use 3DXH if available otherwise 3DRH
n=c_caa_var_get(irf_ssub('Spacecraft_potential__C?_CP_EFW_L2_P',ic),'caa','ts');
   
    
energy = LEEA3Ddef.en;
energyH = HEEA3Ddef.en;
theta = LEEA3Ddef.theta;

%LEEA data
n = n.tlim(irf_time(tint,'epoch>epochTT'));
offset1 = 1.2*nanmean(n.data);
if isnan(offset1),
    offset1 = 0;
end
energytemp = energy+offset1;
[~,idx]=irf_tlim(LEEA3Ddef.tt,tint);
disttempf = squeeze(LEEA3Ddef.data(idx,:,:));
dsf = squeeze(nanmean(disttempf,1));

%HEEA data 
[~,idx]=irf_tlim(HEEA3Ddef.tt,tint);
disttempf = squeeze(HEEA3Ddef.data(idx,:,:));
dshf = squeeze(nanmean(disttempf,1));

%Combine LEEA and HEEA data
elstemp = find(energytemp < energyH(1));
et = cat(2,energytemp(elstemp),energyH);
dstf = cat(2,dsf(:,elstemp),dshf);
et(find(et < 2)) = NaN;

Energy = log10(et);
theta2 = [theta, 180+theta, theta(1)]+90-theta(1);
dstf = [dstf; flipdim(dstf,1); dstf(1,:)];

[TH,R] = meshgrid(theta2,Energy);
[X,Y] = pol2cart(TH*pi/180.0,R);


h=irf_plot(1,'newfigure');
% Plot DIFFERENTIAL ENERGY FLUX 2D plot - Pitch Angle
h(1)=irf_panel('defpad');
pcolor(h(1),double(X),double(Y),double(log10(dstf')))
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
axis(h(1),[-max(Energy) max(Energy) -max(Energy) max(Energy)]);
c=colorbar;
ylabel(c,'log_{10} DEFlux','fontsize',20);
set(c,'fontsize',20)
set(h(1),'position',[0.09 0.2 0.76 0.78]);
set(h(1),'fontsize',20)

if 0, % set to 1 to save plot
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters edefluxpadtest.eps;
end