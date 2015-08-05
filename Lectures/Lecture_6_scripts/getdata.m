%Example script to download Particle data

%%%%%%%%%%%%%%%%%%%%%%%
% specify time interval

tint = [iso2epoch('2007-03-15T08:01:00.000000Z') iso2epoch('2007-03-15T08:10:00.000000Z')];


%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1, % put to 0 if data already downloaded !!!!
    
    % Magnetic fields
    caa_download(tint,'C?_CP_FGM_FULL');
    
    % Spacecraft potential
    caa_download(tint,'C?_CP_EFW_L2_P');
    
    % Ion data
    caa_download(tint,'C?_CP_CIS-HIA_HS_1D_PEF');
    caa_download(tint,'C?_CP_CIS-HIA_HS_MAG_IONS_PEF');
    caa_download(tint,'C?_CP_CIS-HIA_ONBOARD_MOMENTS');
    caa_download(tint,'C?_CP_CIS-HIA_PAD_HS_MAG_IONS_PF');    
    
    % Electron data
    caa_download(tint,'C?_CP_PEA_MOMENTS');
    caa_download(tint,'C?_CP_PEA_PITCH_SPIN_DEFlux');
    caa_download(tint,'C?_CP_PEA_PITCH_SPIN_PSD');
    caa_download(tint,'C?_CP_PEA_PITCH_3DXL_PSD');
    caa_download(tint,'C?_CP_PEA_PITCH_3DXL_DEFlux');
    caa_download(tint,'C?_CP_PEA_PITCH_3DXH_PSD');
    caa_download(tint,'C?_CP_PEA_PITCH_3DXH_DEFlux');
    caa_download(tint,'C?_CP_PEA_PITCH_3DRH_PSD');
    caa_download(tint,'C?_CP_PEA_PITCH_3DRH_DEFlux');
    
    
    download_status=caa_download; % repeat until all data are downloaded
    if download_status==0, % some data are still in queue
      disp('___________!!!!_____________')
      disp('Some data where put in queue!')
      disp('To see when they are ready and to download execute "caa_download".');
      return
    end
end





