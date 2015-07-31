ic = 4;

Eibm = c_caa_var_get(irf_ssub('Data__C?_CP_EFW_L1_IB',ic),'caa','mat');
Bgse = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'caa','mat');
Ebxy = c_caa_var_get(irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L2_EB',ic),'caa','mat');
BDSC=c_coord_trans('GSE','DSC',Bgse,'cl_id',ic);
BDSI=c_coord_trans('GSE','DSI',Bgse,'cl_id',ic);
tint = [Eibm(1,1) Eibm(end,1)];

[tt,phase_data] = irf_isdat_get([irf_ssub('Cluster/?/ephemeris/phase_2',ic)], Eibm(1,1)-10, Eibm(end,1)-Eibm(1,1)+20);
idx = find(phase_data == 0);
phasedata = [tt(idx) [0:1:length(idx)-1]'*360];


probett = Eibm(:,1);
BDSC = irf_resamp(BDSC,probett);
phasedata = irf_resamp(phasedata,probett);
phasedata = [phasedata(:,1) mod(phasedata(:,2),360)];

%probe phases from c_pl_sc_orient
phase_p1=phasedata(:,2)/180*pi + 3*pi/4 ;
phase_p3=phase_p1     - pi/2   ;
phase_p2=phase_p1     + pi     ;
phase_p4=phase_p1     + pi/2 ;
rp1=[44*cos(phase_p1) 44*sin(phase_p1)]; % in DSC reference frame
rp2=[44*cos(phase_p2) 44*sin(phase_p2)];
rp3=[44*cos(phase_p3) 44*sin(phase_p3)];
rp4=[44*cos(phase_p4) 44*sin(phase_p4)];

thetap1b = (rp1(:,1).*BDSC(:,2)+rp1(:,2).*BDSC(:,3))./(sqrt(rp1(:,1).^2+rp1(:,2).^2).*sqrt(BDSC(:,2).^2+BDSC(:,3).^2));
thetap1b = acosd(abs(thetap1b));

thetap3b = (rp3(:,1).*BDSC(:,2)+rp3(:,2).*BDSC(:,3))./(sqrt(rp3(:,1).^2+rp3(:,2).^2).*sqrt(BDSC(:,2).^2+BDSC(:,3).^2));
thetap3b = acosd(abs(thetap3b));

resc = 0.002165;
VL1 = [Eibm(:,1) Eibm(:,2)*resc Eibm(:,3)*resc Eibm(:,4)*resc Eibm(:,5)*resc];
SCV12 = [VL1(:,1) (VL1(:,2)+VL1(:,3))/2];
SCV34 = [VL1(:,1) (VL1(:,4)+VL1(:,5))/2];
E1 = (VL1(:,2)-SCV34(:,2))*1e3/44;
E2 = (SCV34(:,2)-VL1(:,3))*1e3/44;
E3 = (VL1(:,4)-SCV12(:,2))*1e3/44;
E4 = (SCV12(:,2)-VL1(:,5))*1e3/44;

idxB = find(sqrt(BDSC(:,2).^2+BDSC(:,3).^2) < abs(BDSC(:,4)));
thresang = 25.0;

E1(find(thetap1b > thresang)) = NaN;
E2(find(thetap1b > thresang)) = NaN;
E3(find(thetap3b > thresang)) = NaN;
E4(find(thetap3b > thresang)) = NaN;
SCV12(find(thetap3b > thresang),2) = NaN;
SCV34(find(thetap1b > thresang),2) = NaN;
E1(idxB) = NaN;
E2(idxB) = NaN;
E3(idxB) = NaN;
E4(idxB) = NaN;
SCV12(idxB,2) = NaN;
SCV34(idxB,2) = NaN;

%irf_plot([Eibm(:,1) thetap1b thetap3b])
h=irf_plot(6);

h(1)=irf_panel('BGSE');
irf_plot(h(1),Bgse);
ylabel(h(1),'B (GSE) (nT)');
irf_zoom(h(1),'y',[-100,100])
irf_legend(h(1),{'B_x','B_y','B_z'},[0.98 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('BDSC');
irf_plot(h(2),[BDSC(:,1) sqrt(BDSC(:,2).^2+BDSC(:,3).^2) abs(BDSC(:,4))]);
ylabel(h(2),'B (DSC) (nT)');
irf_zoom(h(2),'y',[0,100])
irf_legend(h(2),{'|B_{plane}|','|B_z|'},[0.98 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('theta');
irf_plot(h(3),[Eibm(:,1) thetap1b thetap3b])
ylabel(h(3),'\theta (deg)');
%irf_zoom(h(3),'y',[-100,100])
irf_legend(h(3),{'\theta_{1}','\theta_{3}'},[0.98 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('Vp');
irf_plot(h(4),[SCV12(:,1) SCV12(:,2) SCV34(:,2)])
ylabel(h(4),'V (V/m)');
irf_legend(h(4),{'V_{SC12}','V_{SC34}'},[0.98 0.1])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

h(5)=irf_panel('Ep');
irf_plot(h(5),[Eibm(:,1) E1 E2 E3 E4])
ylabel(h(5),'E (mV/m)');
irf_zoom(h(5),'y',[-30,30])
irf_legend(h(5),{'E_{SC-p1}','E_{p2-SC}','E_{SC-p3}','E_{p4-SC}'},[0.98 0.1])
irf_legend(h(5),'(d)',[0.99 0.98],'color','k')

[Apar,Aperp]=irf_dec_parperp(BDSI,Ebxy,1);

h(6)=irf_panel('Eparperp');
irf_plot(h(6),[Apar(:,1) Apar(:,2) Aperp(:,2)])
ylabel(h(6),'E (mV/m)');
irf_legend(h(6),{'E_{||}','E_{\perp}'},[0.98 0.1])
irf_legend(h(6),'(e)',[0.99 0.98],'color','k')

irf_plot_axis_align(6)
irf_zoom(h,'x',tint);
irf_timeaxis(h);


%Find and print times when the probes satisfy alignment conditions
E1temp = isnan(E1)-[-1; isnan(E1(3:length(E1))); 1];
starttime = Eibm(find(E1temp == 1),1);
starttime = irf_time(starttime,'epoch>iso');

endtime = Eibm(find(E1temp == -1),1);
endtime = irf_time(endtime,'epoch>iso');

fprintf('Intervals when Probe 1 and 2 satisfy alignment conditions \n')
for q = 1:length(starttime(:,1))
    fprintf(strcat('No.',num2str(q),' : ',starttime(q,:),' -- ',endtime(q,:),'\n'))
end

E3temp = isnan(E3)-[-1; isnan(E3(3:length(E3))); 1];
starttime = Eibm(find(E3temp == 1),1);
starttime = irf_time(starttime,'epoch>iso');

endtime = Eibm(find(E3temp == -1),1);
endtime = irf_time(endtime,'epoch>iso');

fprintf('Intervals when Probe 3 and 4 satisfy alignment conditions \n')
for q = 1:length(starttime(:,1))
    fprintf(strcat('No.',num2str(q),' : ',starttime(q,:),' -- ',endtime(q,:),'\n'))
end
    

