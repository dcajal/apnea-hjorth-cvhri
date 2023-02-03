function [tdvol, Me, me, artifacts, Mpeaks, Mlocs, mpeaks, mlocs] = tdvol_est(resp,fs)
    
    minDist = 1.5*fs;
    zerocross = diff(sign(resp));
    downcross = find(zerocross==-2)+1;
    upcross = find(zerocross==2)+1; 
    
    downcross(diff(upcross)<minDist) = [];
    upcross(diff(upcross)<minDist) = [];

    % Clean artifacts
    Setup.plotflag = false;
    Setup.seg = 40;
    Setup.step = 20;
    Setup.thresholdH0Low = 0.05;
    Setup.thresholdH0Up = 30;
    Setup.thresholdH1Low = 0;
    Setup.thresholdH1Up = 2;
    Setup.thresholdH2Low = 0;
    Setup.thresholdH2Up = 5;
    Setup.minSegmentSeparation = 5;
    artifacts = hjorthArtifacts(resp, fs, Setup);
    artifacts = round(artifacts*fs)+1;
    
    
    % Upcross first
    if downcross(1)<upcross(1)
        downcross(1) = [];
    end    
%     for kk=1:size(artifacts,1)       
%         resp(artifacts(kk,1):artifacts(kk,2)) = nan;
%         downcross(downcross>=artifacts(kk,1) & downcross<=artifacts(kk,2)) = [];
%         upcross(upcross>=artifacts(kk,1) & upcross<=artifacts(kk,2)) = [];
% 
%         % Delete last upcross if > last downcross
%         if upcross(find(upcross<=artifacts(kk,1),1,'last')) > downcross(find(downcross<=artifacts(kk,1),1,'last'))
%             upcross(find(upcross<=artifacts(kk,1),1,'last')) = [];
%         end
% 
%         % Delete first downcross if < first upcross
%         if downcross(find(downcross>=artifacts(kk,2),1)) < upcross(find(upcross>=artifacts(kk,2),1))
%             downcross(find(downcross>=artifacts(kk,2),1)) = [];
%         end
%     end
    if downcross(end)<upcross(end)
        upcross(end) = [];
    end   
        
    
%     figure
%     plot(resp); hold on
%     plot(downcross,resp(downcross),'vg')
%     plot(upcross,resp(upcross),'^b')

    
    Mpeaks = nan(size(downcross));
    Mlocs = Mpeaks;
    mpeaks = Mpeaks;
    mlocs = Mpeaks;
    if downcross(1)<upcross(1)
        for kk=2:length(downcross)
            indexes = upcross(kk-1):downcross(kk);
            [Mpeaks(kk), Mlocs(kk)] = max(resp(indexes));
            Mlocs(kk) = Mlocs(kk) + upcross(kk-1) - 1;
        end
        for kk=1:length(upcross)
            [mpeaks(kk), mlocs(kk)] = min(resp(downcross(kk):upcross(kk)));
            mlocs(kk) = mlocs(kk) + downcross(kk) - 1;
        end
    else %downcross(1)>upcross(1)
        for kk=1:length(downcross)
            [Mpeaks(kk), Mlocs(kk)] = max(resp(upcross(kk):downcross(kk)));
            Mlocs(kk) = Mlocs(kk) + upcross(kk) - 1;
        end
        for kk=2:length(upcross)
            [mpeaks(kk), mlocs(kk)] = min(resp(downcross(kk-1):upcross(kk)));
            mlocs(kk) = mlocs(kk) + downcross(kk-1) - 1;
        end
    end
    
    Mpeaks = Mpeaks(~isnan(Mlocs));
    Mlocs = Mlocs(~isnan(Mlocs));
    mpeaks = mpeaks(~isnan(mlocs));
    mlocs = mlocs(~isnan(mlocs));
    
    
%     Me = spline(Mlocs, Mpeaks, 1:length(resp));
%     me = spline(mlocs, mpeaks, 1:length(resp));
    Me = interp1(Mlocs, Mpeaks, 1:length(resp),'pchip');
    me = interp1(mlocs, mpeaks, 1:length(resp),'pchip');
    
    Me(1:upcross(1)) = nan;
    Me(downcross(end):end) = nan;
    me(1:downcross(1)) = nan;
    me(upcross(end):end) = nan;
    
    mpeaks = mpeaks(mlocs>Mlocs(1));
    mlocs = mlocs(mlocs>Mlocs(1));
    
    aux_tdvol = Mpeaks(1:length(mpeaks)) - mpeaks;
    aux_t = (Mlocs(1:length(mpeaks))+mlocs)/2;
    
%     tdvol = spline(aux_t, aux_tdvol, 1:length(resp));
    tdvol = interp1(aux_t, aux_tdvol, 1:length(resp),'pchip');
    tdvol(isnan(Me-me)) = nan;
    
    tdvol(isnan(resp)) = nan;
%     tdvol = Me - me;
    
end