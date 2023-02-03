function [ basal ] = computeBasalRespiration( tidalVolume, spo2, tidalVolumeFs, spo2Fs, artifacts )

Setup.plotflag = false;
Setup.seg = 3*60;
Setup.step = 1*60;
Setup.adaptive = false;
Setup.thresholdH0Low = 0;
Setup.thresholdH0Up = 7;
Setup.thresholdH1Low = 0;
Setup.thresholdH1Up = 10;
Setup.thresholdH2Low = 0;
Setup.thresholdH2Up = 10;
iAbnormalBreathing = hjorthArtifacts(spo2, spo2Fs, Setup);
iAbnormalBreathing = round(iAbnormalBreathing*tidalVolumeFs)+1;

normalBreathing = tidalVolume;
for kk=1:size(iAbnormalBreathing,1)
    normalBreathing(iAbnormalBreathing(kk,1):iAbnormalBreathing(kk,2)) = nan;
end
for kk=1:size(artifacts,1)
    normalBreathing(artifacts(kk,1):artifacts(kk,2)) = nan;
end

% Adaptive update
momentum = 0.6;
updatePeriod = 60*tidalVolumeFs;
basal = nan(size(normalBreathing));
basal(1:updatePeriod) = median(normalBreathing,'omitnan');
for kk = updatePeriod:updatePeriod:numel(normalBreathing)
    newMed = median(normalBreathing(kk+1:min(kk+updatePeriod,end)),'omitnan');
    if isnan(newMed)
        % Full abnormal respiration within segment. Reuse previous value
        newMed = basal(kk);
    end
    basal(kk+1:min(kk+updatePeriod,end)) = basal(kk)*momentum + newMed*(1-momentum);
end

% First part abnormal
if isnan(basal(1))
    basal(1:find(isnan(basal),1,'last')) =  basal(find(~isnan(basal),1));
end

% Adaptive update backwards
basalBackwards = nan(size(normalBreathing));
basalBackwards(1:updatePeriod) = basal(end-updatePeriod+1:end);
cleanRespirationBackwards = flip(normalBreathing);
for kk = updatePeriod:updatePeriod:numel(normalBreathing)
    newMed = median(cleanRespirationBackwards(kk+1:min(kk+updatePeriod,end)),'omitnan');
    if isnan(newMed)
        % Full abnormal respiration within segment. Reuse previous value
        newMed = basalBackwards(kk);
    end
    basalBackwards(kk+1:min(kk+updatePeriod,end)) = basalBackwards(kk)*momentum + newMed*(1-momentum);
end

basal = 0.5*basal + 0.5*flip(basalBackwards);

% figure
% plot(normalBreathing); hold on
% plot(basal);

end

