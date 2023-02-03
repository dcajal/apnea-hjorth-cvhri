clear
addpath(genpath('lib'));

minNEvent = 3; % Minimum number of events to consider a burst
maxNEvent = 8; % Maximum number of events in a burst. At maxNEvent+1 may start a new burst that can be labeled as hypo or apnea independently
maxDistEvent = 90; % Maximum distance between events to consider a burst

% Load all files in database directory
dirlist = dir('dataset/annotations/pneumo/');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end
load("dataset\AHI.mat"); AHI = AHI.AHI;

% For dataset statistics
nAobsDataset = 0;
nAcenDataset = 0;
nAmixDataset = 0;
nHpopDataset = 0;
nHobsDataset = 0;

for kk = 1:length(files)
    subject = split(files{kk},'.');
    subject = subject(1);
    load(strcat('dataset/annotations/pneumo/',string(files{kk})));
    fprintf('Computing subject %s. ',string(subject));

    starts = (pneumo{1,1}-1)/500;
    ends = starts+pneumo{1,3};
    events = [starts ends];
    types = pneumo{1,2};
    apneas = events(types<5,:);
    hypopneas = events(types>4,:);
    doubts = [];

    % For dataset statistics
    nAobsDataset = nAobsDataset+sum(types==2);
    nAcenDataset = nAcenDataset+sum(types==3);
    nAmixDataset = nAmixDataset+sum(types==4);
    nHpopDataset = nHpopDataset+sum(types==5);
    nHobsDataset = nHobsDataset+sum(types==6);
    fprintf('AHI: %.2f. ',AHI(kk));
    fprintf('nAobs: %i, nAcen: %i, nAmix: %i, nHpop: %i, nHobs: %i... ' ...
        ,sum(types==2),sum(types==3),sum(types==4),sum(types==5),sum(types==6));
    
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
            if counter == minNEvent
                init = ee-minNEvent+1;
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
            if counter >= minNEvent || continues
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
    
    
    save(strcat('results/labels/',string(subject),'_labels.mat'),...
        'apneas','hypopneas','doubts','apneaSegments','hypoSegments','abnormalSegments');
    
    fprintf('Done\n');

end

nEventsDataset = nAobsDataset+nAcenDataset+nAmixDataset+nHpopDataset+nHobsDataset;
fprintf('\n'); fprintf('Database statistics.\n');
fprintf('nAobs: %i (%.1f%%), nAcen: %i (%.1f%%), nAmix: %i (%.1f%%), nHpop: %i (%.1f%%), nHobs: %i (%.1f%%)\n', ...
        nAobsDataset,nAobsDataset/nEventsDataset*100, ...
        nAcenDataset,nAcenDataset/nEventsDataset*100, ...
        nAmixDataset,nAmixDataset/nEventsDataset*100, ...
        nHpopDataset,nHpopDataset/nEventsDataset*100, ...
        nHobsDataset,nHobsDataset/nEventsDataset*100);
fprintf('AHI<5: %i (%.1f%%), 5<=AHI<15: %i (%.1f%%), 15<=AHI<30: %i (%.1f%%), AHI>=30: %i (%.1f%%)\n', ...
        sum(AHI<5),sum(AHI<5)/numel(AHI)*100, ...
        sum(AHI>=5 & AHI<15),sum(AHI>=5 & AHI<15)/numel(AHI)*100, ...
        sum(AHI>=15 & AHI<30),sum(AHI>=15 & AHI<30)/numel(AHI)*100, ...
        sum(AHI>30),sum(AHI>30)/numel(AHI)*100);
