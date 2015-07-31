
%%%%%%%%%%%%%%%%%%%%%%%
% specify time interval
tint = [iso2epoch('2002-02-06T08:54:43.000000Z') iso2epoch('2002-02-06T08:54:56.000000Z')]; % time interval
tint2 = [iso2epoch('2002-02-06T08:45:00.000000Z') iso2epoch('2002-02-06T09:05:00.000000Z')]

%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1, % put to 0 if data already downloaded !!!!


    caa_download(tint,'C2_CP_EFW_L2_EB')
    caa_download(tint,'C2_CP_EFW_L2_PB')
    caa_download(tint,'C2_CP_EFW_L1_IB')
    caa_download(tint,'C3_CP_EFW_L2_EB')
    caa_download(tint,'C3_CP_EFW_L2_PB')
    caa_download(tint,'C3_CP_EFW_L1_IB')
    caa_download(tint2,'C2_CP_FGM_FULL')
    caa_download(tint2,'C2_CP_FGM_FULL_ISR2')
    caa_download(tint2,'C3_CP_FGM_FULL')
    caa_download(tint2,'C3_CP_FGM_FULL_ISR2')
    
    download_status=caa_download; % repeat until all data are downloaded
    if download_status==0, % some data are still in queue
      disp('___________!!!!_____________')
      disp('Some data where put in queue!')
      disp('To see when they are ready and to download execute "caa_download".');
      return
    end
end





