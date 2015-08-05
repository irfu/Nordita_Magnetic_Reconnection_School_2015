

%%%%%%%%%%%%%%%%%%%%%%%
% specify time interval

tintC2bm = [iso2epoch('2004-03-08T16:03:17.000000Z') iso2epoch('2004-03-08T16:03:30.000000Z')];
tintC4bm = [iso2epoch('2004-03-08T15:55:24.000000Z') iso2epoch('2004-03-08T15:55:36.000000Z')];

tintB = [iso2epoch('2004-03-08T15:45:0.000000Z') iso2epoch('2004-03-08T16:15:0.000000Z')];


%%%%%%%%%%%%%%%%%%%%%%%%
% download data from CAA (needed only once!!!!!)
if 1, % put to 0 if data already downloaded !!!!


    
    caa_download(tintC2bm,'C2_CP_EFW_L2_EB');
    caa_download(tintC2bm,'C2_CP_EFW_L2_BB');
    
    caa_download(tintC4bm,'C4_CP_EFW_L2_EB');
    caa_download(tintC4bm,'C4_CP_EFW_L2_BB');
    
    
    caa_download(tintB,'C?_CP_FGM_FULL_ISR2');    

    
    
    download_status=caa_download; % repeat until all data are downloaded
    if download_status==0, % some data are still in queue
      disp('___________!!!!_____________')
      disp('Some data where put in queue!')
      disp('To see when they are ready and to download execute "caa_download".');
      return
    end
end





