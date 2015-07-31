ic = 2; %set spacecraft number 1--4.

tint = [irf_time('2007-04-10T11:20:25.998125Z','iso>epoch') irf_time('2007-04-10T11:20:26.576791Z','iso>epoch')];

VL1 = c_caa_var_get(irf_ssub('Data__C?_CP_EFW_L1_IB',ic),'caa','mat'); 
%N.B. Different L1 data format because of B field
VL1 = [VL1(:,1) VL1(:,6) VL1(:,7) VL1(:,8) VL1(:,9)];
B = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'caa','mat'); 
[power,freq,wavenumber] = wkpowerspec(VL1, B, tint, ic, 3);

freq = freq*1e-3; % covert from Hz to kHz

h=irf_plot(1);
h(1)=irf_panel('disprel');
pcolor(h(1),wavenumber,freq,log10(power))
shading(h(1),'flat')
ylabel(h(1),'f (kHz)','fontsize',20);
xlabel(h(1),'k_{||} (m^{-1})','fontsize',20);
set(h(1),'fontsize',20)
c=colorbar;
ylabel(c,'log_{10} P(f,k)/P_{max}','fontsize',20);
set(c,'fontsize',20)
set(h(1),'position',[0.09 0.2 0.76 0.78]);

if 0, 
% Overplot linear dispersion relation fit to the data
hold(h(1),'on')
kfit = [0.0001:0.0001:0.1];
v1 = 120; % speed in km s^{-1}
ffit1 = v1/(2*pi)*kfit*1e3;
plot(h(1),kfit,ffit1/1000.0,'linewidth',3,'color','r')
irf_legend(h(1),strcat('v = ',num2str(v1(1)),'km s^{-1}'),[0.6 0.9],'fontsize',20,'color','r')
% if wave number of power spectrum are negative use -kfit in plot
hold(h(1),'off');
end

if 0,
 set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
 print -depsc -painters wkpowerspec.eps;
end