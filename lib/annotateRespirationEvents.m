function [ apneas, hypopneas, doubts ] = annotateRespirationEvents( tidalVolume, fs, basalRespiration, spo2, tSpo2, hypno, computeDoubts )
%APNEALABELS Onsets and ensets of apnea events based on respiration signal
%reductions. Events vector contains all type of events while labels contain
%their type (1 for apnoea, 2 for hypo, 3 for doubt). Units are in seconds

    apneaReduction = 0.1;
    hypopneaReduction = 0.7;
    doubtReduction = 0.3;

    apneaThreshold = basalRespiration*apneaReduction;
    hypopneaThreshold = basalRespiration*hypopneaReduction;
    doubtThreshold = basalRespiration*doubtReduction;

    apneas = computeEvents(tidalVolume, apneaThreshold, fs, hypno);
    hypopneas = computeEvents(tidalVolume, hypopneaThreshold, fs, hypno);
    if computeDoubts
        doubts = computeEvents(tidalVolume, doubtThreshold, fs, hypno);

        noHypo = zeros(size(hypopneas,1),1);
        for kk = 1:numel(noHypo)
            for jj = 1:size(doubts,1)
                noHypo(kk) = hypopneas(kk,1)<=doubts(jj,1) & hypopneas(kk,2)>=doubts(jj,2);
                if noHypo(kk), break; end
            end
        end
        hypopneas = hypopneas(~noHypo,:); % Remove doubts (at this stage doubts still include apneas)
        
        noDoubt = zeros(size(doubts,1),1);
        for kk = 1:numel(noDoubt)
            for jj = 1:size(apneas,1)
                noDoubt(kk) = doubts(kk,1)<=apneas(jj,1) && doubts(kk,2)>=apneas(jj,2);
                if noDoubt(kk), break; end
            end
        end
        doubts = doubts(~noDoubt,:);
    else       
        noHypo = zeros(size(hypopneas,1),1);
        for kk = 1:numel(noHypo)
            for jj = 1:size(apneas,1)
                noHypo(kk) = hypopneas(kk,1)<=apneas(jj,1) && hypopneas(kk,2)>=apneas(jj,2);
                if noHypo(kk), break; end
            end
        end
        hypopneas = hypopneas(~noHypo,:);
        doubts = [];
    end


    % Hypopneas must be related to arousals or desaturations >= 3%
    if isempty(hypopneas), return; end
    windowMargin = 5*fs; % Desaturation may appear few seconds after restored breathing
    hypoNoDesat = false(1,size(hypopneas,1));
    for kk = 1:size(hypopneas,1)
        indexes = find(tSpo2>hypopneas(kk,1),1):(find(tSpo2<hypopneas(kk,2),1,'last')+windowMargin);
        indexes = indexes(indexes<=length(spo2)); % indexes must be constrained to spo2 length
        desaturation = max(spo2(indexes))-min(spo2(indexes));
        if desaturation < 3
            hypoNoDesat(kk) = true;
        end
    end; clear kk indexes desaturation

%     figure
%     plot(tSpo2,spo2); hold on;
%     for ll=1:size(hypopneas,1)
%         if hypoNoDesat(ll)
%             p(1) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[80 80 100 100],[1 0 0], ...
%                 'FaceAlpha',.3,'EdgeColor','none');
%         else
%             p(2) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[80 80 100 100],[0 1 0], ...
%                 'FaceAlpha',.3,'EdgeColor','none');
%         end
%     end; clear ll

    hypopneas(logical(hypoNoDesat),:) = [];


end

function [events] = computeEvents( tidalVolume, threshold, fs, hypno )
    minTime = 10; % seconds under reduction
    maxTime = 90; % seconds under reduction (consider artifact)
    
    respirationBelowThreshold = tidalVolume<threshold; 
    respirationBelowThreshold(1) = 0; % First sample above threshold
    reductions = find(diff(respirationBelowThreshold)>0);
    raises = find(diff(respirationBelowThreshold)<0);
    raises = raises(1:numel(reductions));
    reductionsDuration = raises-reductions;
    events = [reductions(reductionsDuration>=minTime*fs & reductionsDuration<=maxTime*fs)'...
        raises(reductionsDuration>=minTime*fs & reductionsDuration<=maxTime*fs)'];
    events = (events-1)/fs; 
    
    % No apnea events
    if isempty(events), return; end
    
    % Last event not allowed within the last 10 seconds
    if events(end) >= (numel(tidalVolume)-1)/fs - 10
        events(end,:) = [];
    end

    % No apnea events
    if isempty(events), return; end

    % Remove events during Wake (if any)
    tHypno = linspace(0,(length(tidalVolume)-1)/fs,numel(hypno));
    wakeEvents = false(1,size(events,1));
    for kk = 1:size(events,1)
        iHypno = find(tHypno>events(kk,1),1):find(tHypno<events(kk,2),1,'last');
        if sum(hypno(iHypno)==5) > 0 % Wake during event
            wakeEvents(kk) = true;
        end
    end
    events(logical(wakeEvents),:) = [];
    
end