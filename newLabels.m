% Annotate database. processPSG must be executed fist to obtain
% preprocessed PSG signals.

clear
addpath(genpath('lib'));
plotflag = false;

minNEvent = 2; % Minimum number of events to consider a burst
maxNEvent = 8; % Maximum number of events in a burst. At maxNEvent+1 may start a new burst that can be labeled as hypo or apnea independently
maxDistEvent = 180; % Maximum distance between events to consider a burst

% Load all files in database directory
dirlist = dir('dataset/nasal_pressure_signals/');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end


for kk = 1:length(files)

    % Load PSG signals
    subject = split(files{kk},'_');
    subject = subject(1);
    fprintf('Computing subject: %s...',string(subject));
    load(strcat('results/signals/',string(subject),'_psg.mat'),'nasalPressureArtifacts','nasalPressureFs', ...
        'nasalPressureProcessed','tNasalPressure','spo2Processed','spo2Fs','tSpo2','hypno','tHypno','tidalVolume');

    % Labeling    
    basalRespiration = computeBasalRespiration(tidalVolume,spo2Processed,nasalPressureFs,spo2Fs, nasalPressureArtifacts);
    
    [apneas, hypopneas, doubts] = ...
        annotateRespirationEvents(tidalVolume,nasalPressureFs,basalRespiration,spo2Processed,tSpo2,hypno,true);

    % Remove labels during nasal pressure artifacts
    nasalPressureArtifactsSeconds = (nasalPressureArtifacts-1)/nasalPressureFs;
    for ll=1:size(nasalPressureArtifactsSeconds,1)
        for jj = 1:size(apneas,1)
            if ((apneas(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (apneas(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (apneas(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (apneas(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                apneas(jj,:) = nan;
            end   
        end
        for jj = 1:size(hypopneas,1)
            if ((hypopneas(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (hypopneas(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (hypopneas(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (hypopneas(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                hypopneas(jj,:) = nan;
            end   
        end
        for jj = 1:size(doubts,1)
            if ((doubts(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (doubts(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (doubts(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (doubts(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                doubts(jj,:) = nan;
            end   
        end
    end
    apneas(isnan(apneas(:,1)),:) = [];
    hypopneas(isnan(hypopneas(:,1)),:) = [];
    doubts(isnan(doubts(:,1)),:) = [];

    
    % labels: 1=apnea, 2=hypopneas, 3=doubts
    events = [apneas; hypopneas; doubts];
    labels = [ones(size(apneas,1),1); 2*ones(size(hypopneas,1),1); 3*ones(size(doubts,1),1)];
    [events,indexes] = sort(events);
    indexes = indexes(:,1);
    labels = labels(indexes);
    
    counter = 0;
    position = 1;
    abnormalSegments = [];
    nApnea = [];
    nHypo = [];
    nDoubt = [];
    continues = false;
    for ee=1:size(events,1)-1
        if events(ee+1,1)-events(ee,2) <= maxDistEvent && ee<size(events,1)-1
            counter = counter+1;
            if counter == minNEvent-1
                init = ee-minNEvent+2;
                abnormalSegments(position,1) = events(init,1); %#ok<*SAGROW> % burst onset
            elseif counter == maxNEvent
                abnormalSegments(position,2) = events(ee,2); % burst endset
                nApnea(position) = sum(labels(init:ee)==1);
                nHypo(position) = sum(labels(init:ee)==2);
                nDoubt(position) = sum(labels(init:ee)==3);
                position = position+1;
                
                init = ee+1;
                abnormalSegments(position,1) = events(init,1); %#ok<*SAGROW> % burst onset
                counter = 0;
                continues = true;
            end                
        else
            if counter >= minNEvent-1 || continues
                abnormalSegments(position,2) = events(ee,2); % burst endset
                nApnea(position) = sum(labels(init:ee)==1);
                nHypo(position) = sum(labels(init:ee)==2);
                nDoubt(position) = sum(labels(init:ee)==3);
                position = position+1;
            end
            counter = 0;
            continues = false;
        end
    end 
    clear ee counter position
    if numel(abnormalSegments) == 1
        abnormalSegments(1,2) = events(end); % All events in a burst (rare)
    end
    
    
    % Burst labeling: bursts are labeled as apneic if contain at least one
    % apnea or more than a half of the events are "doubts"
    apneaSegments = [];
    hypoSegments = [];
    if ~isempty(abnormalSegments)
        for aa=1:size(abnormalSegments,1)
            doubtRatio = nDoubt(aa)/(nApnea(aa)+nHypo(aa)+nDoubt(aa));
            if nApnea(aa)>0 || doubtRatio>0.5
               apneaSegments = [apneaSegments; abnormalSegments(aa,:)]; %#ok<*AGROW>
            else
               hypoSegments = [hypoSegments; abnormalSegments(aa,:)];
            end
        end
    end
    
    % Remove labels during nasal pressure artifacts
    for ll=1:size(nasalPressureArtifactsSeconds,1)
        for jj = 1:size(apneaSegments,1)
            if ((apneaSegments(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (apneaSegments(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (apneaSegments(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (apneaSegments(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                apneaSegments(jj,:) = nan;
            end   
        end
        for jj = 1:size(hypoSegments,1)
            if ((hypoSegments(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (hypoSegments(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (hypoSegments(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (hypoSegments(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                hypoSegments(jj,:) = nan;
            end   
        end
    end
    if ~isempty(apneaSegments)
        apneaSegments(isnan(apneaSegments(:,1)),:) = [];
    end
    if ~isempty(hypoSegments)
        hypoSegments(isnan(hypoSegments(:,1)),:) = [];
    end

    
    % Plots
    if plotflag
    
        nasalPressureWithoutArtifacts = nasalPressureProcessed; %#ok<UNRCH> 
        for ll=1:size(nasalPressureArtifacts,1)
            nasalPressureWithoutArtifacts(nasalPressureArtifacts(ll,1):nasalPressureArtifacts(ll,2)) = nan;
        end
    
        figure;
        ax(1) = subplot(511);
        hold on
        for ll=1:size(apneas,1)
            p(1) = patch([apneas(ll, 1) apneas(ll, 2) apneas(ll, 2) apneas(ll, 1)],[-6 -6 15 15],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(hypopneas,1)
            p(2) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[-6 -6 15 15],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(doubts,1)
            p(3) = patch([doubts(ll, 1) doubts(ll, 2) doubts(ll, 2) doubts(ll, 1)],[-6 -6 15 15],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        plot(tNasalPressure,nasalPressureProcessed,'r')
        plot(tNasalPressure,nasalPressureWithoutArtifacts,'b')
        plot(tNasalPressure,basalRespiration+5,'k--');
        plot(tNasalPressure,tidalVolume+5);
        yline(5)
        ylabel('Nasal pressure')
        axis tight;
        
        ax(2) = subplot(512);
        plot(tSpo2,spo2Processed,'b'); axis tight;
        hold on
        for ll=1:size(apneas,1)
            p(1) = patch([apneas(ll, 1) apneas(ll, 2) apneas(ll, 2) apneas(ll, 1)],[80 80 100 100],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(hypopneas,1)
            p(2) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[80 80 100 100],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(doubts,1)
            p(3) = patch([doubts(ll, 1) doubts(ll, 2) doubts(ll, 2) doubts(ll, 1)],[80 80 100 100],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        ylabel('SpO2')
        
        ax(3) = subplot(513);
        plot(tAbdBelt,ripSum,'b'); axis tight;
        hold on
        for ll=1:size(apneas,1)
            p(1) = patch([apneas(ll, 1) apneas(ll, 2) apneas(ll, 2) apneas(ll, 1)],[-10 -10 10 10],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(hypopneas,1)
            p(2) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[-10 -10 10 10],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(doubts,1)
            p(3) = patch([doubts(ll, 1) doubts(ll, 2) doubts(ll, 2) doubts(ll, 1)],[-10 -10 10 10],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        ylabel('RIPsum')
        
        ax(4) = subplot(514);
        p = plot(tHypno,hypno,'b'); axis tight;
        hold on
        for ll=1:size(apneas,1)
            p(1) = patch([apneas(ll, 1) apneas(ll, 2) apneas(ll, 2) apneas(ll, 1)],[0 0 6 6],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(hypopneas,1)
            p(2) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[0 0 6 6],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(doubts,1)
            p(3) = patch([doubts(ll, 1) doubts(ll, 2) doubts(ll, 2) doubts(ll, 1)],[0 0 6 6],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        xlabel('Time (seconds)');
        ylabel('Hypnogram')
        yticks(1:5);
        yticklabels({'NREM3','NREM2','NREM1','REM','WAKE'})

        ax(5) = subplot(515);
        hold on
        for ll = 1:size(apneaSegments,1)
            p(1) = patch([apneaSegments(ll, 1) apneaSegments(ll, 2) apneaSegments(ll, 2) apneaSegments(ll, 1)],[0 0 2 2],[1 0 0],'EdgeColor','none');
        end; clear ll
        for ll = 1:size(hypoSegments,1)
            p(2) = patch([hypoSegments(ll, 1) hypoSegments(ll, 2) hypoSegments(ll, 2) hypoSegments(ll, 1)],[0 0 2 2],[0 0 1],'EdgeColor','none');
        end; clear ll
        ylabel('Segment labels')
        axis tight;
    
        linkaxes(ax,'x')
        pause
    end

       
    % AHI
    for ll=1:size(nasalPressureArtifactsSeconds,1) % Exclude nasal pressure artifacts from sleeptime
        indexes = find(tHypno>=nasalPressureArtifactsSeconds(ll,1),1):...
           find(tHypno<=nasalPressureArtifactsSeconds(ll,2),1,'last');
        hypno(indexes) = nan; 
    end
    sleeptime = tHypno(end)*mean(hypno<5);
    ah = size(apneas,1)+size(hypopneas,1)+size(doubts,1);
    ahi = ah/sleeptime*3600;


    save(strcat('results/labels/',string(subject),'_newlabels.mat'),...
        'apneas','hypopneas','doubts','apneaSegments','hypoSegments','abnormalSegments', ...
        'sleeptime','ah','ahi');
    
    fprintf('Done\n');

end
