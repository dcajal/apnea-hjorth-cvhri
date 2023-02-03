%% Model selection

clear
addpath(genpath('lib'));
rng('default');

tic
warning('off')

% Load all files in database directory
dirlist = dir('results/ppi/');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

% Compute Hj�rth parameters and labels for every subject
fprintf('Computing Hj�rth parameters and labels for every subject\n');
labels = [];
hjorthParameters = [];
subjects = [];
for kk = 1:length(files)
    subject = split(files{kk},'_');
    subject = subject(1);
    fprintf('Computing subject: %s...',string(subject));
    load(strcat('results/labels/',string(subject),'_newlabels.mat'));
    load(strcat('results/signals/',string(subject),'_psg.mat'),'nasalPressureArtifacts','nasalPressureFs', ...
        'spo2Processed','spo2Fs','hypno','tHypno');
    load(strcat('results/ppi/',string(subject),'_ppi.mat'));
     

    % Compute Hjorth parameters
    [hjorthParametersAux, segments] = computeHjorthParameters(spo2Processed,spo2Fs,ppi,ppiFs);
    
    % Compute labels
    labelsAux = computeSegmentLabels(segments, abnormalSegments, apneaSegments, hypoSegments);

    % Remove labels during nasal pressure artifacts
    nasalPressureArtifactsSeconds = (nasalPressureArtifacts-1)/nasalPressureFs;
    for ll=1:size(nasalPressureArtifactsSeconds,1)
        for jj = 1:size(segments,1)
            if ((segments(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (segments(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (segments(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (segments(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                hjorthParametersAux(jj,:) = nan;
            end   
        end
    end
    
    labelsAux(logical(sum(isnan(hjorthParametersAux),2)),:) = []; % Delete NaNs
    hjorthParametersAux(logical(sum(isnan(hjorthParametersAux),2)),:) = []; % Delete NaNs
            
    labels = [labels; labelsAux]; %#ok<*AGROW> 
    hjorthParameters = [hjorthParameters; hjorthParametersAux];
    subjects = [subjects; kk*ones(size(labelsAux))];
    fprintf('Done\n');

end

% Random class balancing
fprintf('Random class balancing... ');
[hjorthParameters,labels] = balanceRandom(hjorthParameters,labels,subjects);
fprintf('Done\n');

learnerMatrix = [hjorthParameters labels];

toc

% Select the best model using Classification Learner app with 5-folds based
% on AUC

%% Validation preprocessing
% Leave one out strategy

clear
addpath(genpath('lib'));
addpath(genpath('models'));
rng('default');
plotflag = false;
tic

warning('off')

% Load all files in database directory
dirlist = dir('results/ppi/');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

% Compute parameters and labels for every subject
fprintf('Computing Hj�rth parameters and labels for every subject...\n');
hjorthParameters = [];
labels = [];
subjects = [];
hypnoLabels = [];
for kk = 1:length(files)
    subject = split(files{kk},'_');
    subject = subject(1);
    fprintf('Computing subject: %s...',string(subject));
    load(strcat('results/labels/',string(subject),'_newlabels.mat'));
    load(strcat('results/signals/',string(subject),'_psg.mat'),'nasalPressureArtifacts','nasalPressureFs', ...
        'spo2Processed','spo2Fs','hypno','tHypno');
    load(strcat('results/ppi/',string(subject),'_ppi.mat'));


    % Compute Hjorth parameters
    [hjorthParametersAux, segments] = computeHjorthParameters(spo2Processed,spo2Fs,ppi,ppiFs);

    % Compute labels
    labelsAux = computeSegmentLabels(segments, abnormalSegments, apneaSegments, hypoSegments);
    hypnoLabelsAux = computeHypnoLabels(segments,hypno,tHypno);

    % Remove labels during nasal pressure artifacts
    nasalPressureArtifactsSeconds = (nasalPressureArtifacts-1)/nasalPressureFs;
    for ll=1:size(nasalPressureArtifactsSeconds,1)
        for jj = 1:size(segments,1)
            if ((segments(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (segments(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (segments(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (segments(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                hjorthParametersAux(jj,:) = nan;
            end   
        end
    end

    labelsAux(logical(sum(isnan(hjorthParametersAux),2)),:) = []; % Delete NaNs
    hypnoLabelsAux(logical(sum(isnan(hjorthParametersAux),2)),:) = []; % Delete NaNs
    hjorthParametersAux(logical(sum(isnan(hjorthParametersAux),2)),:) = []; % Delete NaNs
    
    hjorthParameters = [hjorthParameters; hjorthParametersAux];
    labels = [labels; labelsAux];
    hypnoLabels = [hypnoLabels; hypnoLabelsAux];
    subjects = [subjects; kk*ones(size(labelsAux))];

    fprintf('Done\n');

end; clear labelsAux hjorthParametersAux hypnoLabelsAux kk
fprintf('\n');

% Random class balancing
fprintf('Random class balancing... ');
[hjorthParameters,labels,subjects] = balanceRandom(hjorthParameters,labels,subjects);
fprintf('Done\n');

toc

%% Validation using LOO strategy

tic

disp('Train and test using LOO strategy');
predictionsTest = [];
labelsTest = [];
cvhriTest = []; % new index to check correlation
subjectTest = [];
corrects = []; % rate of correct detections
for kk = 1:length(files)

    % Train
    trainedClassifier = trainClassifier([hjorthParameters(subjects~=kk,:) labels(subjects~=kk)]);
       
    % Test subject out
    subject = split(files{kk},'_');
    subject = subject(1);
    fprintf('Computing subject: %s...',string(subject));
    load(strcat('results/labels/',string(subject),'_newlabels.mat'));
    if plotflag
        load(strcat('results/signals/',string(subject),'_psg.mat'),'nasalPressureArtifacts','nasalPressureFs', ...
            'nasalPressureProcessed','tNasalPressure','spo2Processed','spo2Fs','tSpo2','hypno','tHypno','tidalVolume');
    else
        load(strcat('results/signals/',string(subject),'_psg.mat'),'nasalPressureArtifacts','nasalPressureFs', ...
            'spo2Processed','spo2Fs');
    end
    load(strcat('results/ppi/',string(subject),'_ppi.mat'));

    % Compute Hjorth parameters
    [hjorthParametersTestAux, segments] = computeHjorthParameters(spo2Processed,spo2Fs,ppi,ppiFs);

    % Compute labels
    labelsTestAux = computeSegmentLabels(segments, abnormalSegments, apneaSegments, hypoSegments);

    % Remove labels during nasal pressure artifacts
    nasalPressureArtifactsSeconds = (nasalPressureArtifacts-1)/nasalPressureFs;
    for ll=1:size(nasalPressureArtifactsSeconds,1)
        for jj = 1:size(segments,1)
            if ((segments(jj,1) > nasalPressureArtifactsSeconds(ll,1)) && (segments(jj,1) < nasalPressureArtifactsSeconds(ll,2))) || ...
                    (segments(jj,2) > nasalPressureArtifactsSeconds(ll,1)) && (segments(jj,2) < nasalPressureArtifactsSeconds(ll,2))
                hjorthParametersTestAux(jj,:) = nan;
            end   
        end
    end
    
    labelsTestAux(logical(sum(isnan(hjorthParametersTestAux),2)),:) = []; % Delete NaNs
    segments(logical(sum(isnan(hjorthParametersTestAux),2)),:) = []; % Delete NaNs
    hjorthParametersTestAux(logical(sum(isnan(hjorthParametersTestAux),2)),:) = []; % Delete NaNs
    
    % Predict
    predictionsTestAux = trainedClassifier.predictFcn(hjorthParametersTestAux);
    predictionsTest = [predictionsTest; predictionsTestAux];
    labelsTest = [labelsTest; labelsTestAux];
    subjectTest = [subjectTest; ones(size(labelsTestAux))*kk];
    
    correctsAux = mean(~xor(labelsTestAux,predictionsTestAux));
    corrects = [corrects; correctsAux];

    % CVHRI
    segmentIndexes = floor(segments*ppiFs)+1;
    cvhriTestAux = nan(1,size(segments,1));
    for ll = 1:size(segmentIndexes,1)
        if predictionsTestAux(ll) > 0
            segmentFFT = abs(fft(detrend(ppi(segmentIndexes(ll,1):segmentIndexes(ll,2))),2^15));
            segmentFFT = segmentFFT(1:length(segmentFFT)/2);
            f = (linspace(0,0.5,numel(segmentFFT)))*ppiFs;
            [peak,cvhriTestAux(ll)] = findpeaks(segmentFFT(f<0.1),f(f<0.1),'NPeaks',1,'SortStr','descend');
        end
    end; clear ll aux
    cvhriTestAux = sum(cvhriTestAux,'omitnan')/size(segmentIndexes,1);
    cvhriTest = [cvhriTest; cvhriTestAux];

    % Plots    
    if plotflag

        nasalPressureWithoutArtifacts = nasalPressureProcessed;
        for ll=1:size(nasalPressureArtifacts,1)
            nasalPressureWithoutArtifacts(nasalPressureArtifacts(ll,1):nasalPressureArtifacts(ll,2)) = nan;
        end

        figure
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
        plot(tNasalPressure,medfilt1(tidalVolume,180*nasalPressureFs)+5,'k--');
        plot(tNasalPressure,tidalVolume+5);
        yline(5)
        ylabel('Nasal pressure')
        axis tight;

        ax(2) = subplot(512);
        plot(tSpo2,spo2Processed,'b'); axis tight;
        hold on
        for ll=1:size(apneas,1)
            p(1) = patch([apneas(ll, 1) apneas(ll, 2) apneas(ll, 2) apneas(ll, 1)],[90 90 100 100],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(hypopneas,1)
            p(2) = patch([hypopneas(ll, 1) hypopneas(ll, 2) hypopneas(ll, 2) hypopneas(ll, 1)],[90 90 100 100],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        for ll=1:size(doubts,1)
            p(3) = patch([doubts(ll, 1) doubts(ll, 2) doubts(ll, 2) doubts(ll, 1)],[90 90 100 100],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
        end; clear ll
        ylabel('SpO2')

        ax(3) = subplot(513);
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
        ylabel('Hypnogram')
        yticks(1:5);
        yticklabels({'NREM3','NREM2','NREM1','REM','WAKE'})

        ax(4) = subplot(514);
        hold on
        for ll = 1:size(segments,1)
            if labelsTestAux(ll)==1 % apnea
                p(4) = patch([segments(ll, 1) segments(ll, 2) segments(ll, 2) segments(ll, 1)],[0 0 2 2],[1 0 0],'EdgeColor','none');
            end
            if labelsTestAux(ll)==2 % hypo
                p(4) = patch([segments(ll, 1) segments(ll, 2) segments(ll, 2) segments(ll, 1)],[0 0 2 2],[0 0 1],'EdgeColor','none');
            end
        end; clear ll
        plot(tPpi,ppi)
        ylabel('Labels')
        axis tight;

        ax(5) = subplot(515);
        hold on
        for ll = 1:size(segments,1)
            if predictionsTestAux(ll)==1 % apnea
                p(2) = patch([segments(ll, 1) segments(ll, 2) segments(ll, 2) segments(ll, 1)],[0 0 2 2],[1 0 0],'EdgeColor','none');
            end
            if predictionsTestAux(ll)==2 % hypo
                p(2) = patch([segments(ll, 1) segments(ll, 2) segments(ll, 2) segments(ll, 1)],[0 0 2 2],[0 0 1],'EdgeColor','none');
            end
        end; clear ll
        plot(tPpi,ppi)
        ylabel('Predicted')

        axis tight;
        xlabel('Time (seconds)');
        linkaxes(ax,'x');
        set(gcf, 'Position', get(0, 'Screensize'));

        pause;
    end

    fprintf('Done\n');
end; clear kk
    
toc

c = confusionmat(labelsTest,predictionsTest);
plotconfusion(labelsTest',predictionsTest');


targets = zeros(3,length(labelsTest));
outputs = zeros(3,length(predictionsTest));
targetsIdx = sub2ind(size(targets), labelsTest'+1, 1:length(labelsTest));
outputsIdx = sub2ind(size(outputs), predictionsTest'+1, 1:length(predictionsTest));
targets(targetsIdx) = 1;
outputs(outputsIdx) = 1;
plotconfusion(targets,outputs);

%% Exclude test errors

% Load AHI
load(strcat('results/AHI_new.mat'));

% 33 was wake most of the time
exclude = 33; 
labelsTest(subjectTest~=kk,:) = nan;
predictionsTest
ahiDataset

%% Segment detection

figure;
stem(corrects*100); ylim([0 100])
xlabel('Subject'); ylabel('Correct detections (%)')
% hold on; yline(mean(corrects)*100,'--'); 
hold on; plot(xlim,[mean(corrects)*100 mean(corrects)*100],'--'); 
mean(corrects)

% figure;
% cm = confusionchart(labelsTest,predictionsTest);
cm = confusionmat(labelsTest,predictionsTest)
% cmValues = cm.NormalizedValues;
cmValues = cm;
precision = 100*cmValues(2,2)/(cmValues(2,2)+cmValues(2,1)) % TP/(TP+FP)
recall = 100*cmValues(2,2)/(cmValues(2,2)+cmValues(1,2)) % TP/(TP+FN)
accuracy = 100*(cmValues(1,1)+cmValues(2,2))/(cmValues(1,1)+cmValues(2,2)+cmValues(1,2)+cmValues(2,1)) % nCorrect/nTotal

%% Correlation

[rhoPearson,pvalPearson] = corr(ahiDataset,cvhriTest,'Type','Pearson') %#ok<*ASGLU> 
[rhoKendall,pvalKendall] = corr(ahiDataset,cvhriTest,'Type','Kendall')
[rhoSpearman,pvalSpearman] = corr(ahiDataset,cvhriTest,'Type','Spearman')

figure;
plot(ahi,cvhriTest,'o'); 
xlabel('AHI'); ylabel('CVHRI')
hold on; xline(15,'--k');% yline(1/15,'--k')
hold on; xline(30,'--k');% yline(1/30,'--k')
xlim([0 70])
xticks(0:10:70)
% hold on; plot([15 15],ylim,'--');% plot(xlim,[1/15 1/15],'--')

%%

[rhoPearson,pvalPearson] = corr(ahiDataset(ahiDataset<15),cvhriTest(ahiDataset<15),'Type','Pearson')
[rhoKendall,pvalKendall] = corr(ahiDataset(ahiDataset<15),cvhriTest(ahiDataset<15),'Type','Kendall')
[rhoSpearman,pvalSpearman] = corr(ahiDataset(ahiDataset<15),cvhriTest(ahiDataset<15),'Type','Spearman')

[rhoPearson,pvalPearson] = corr(ahiDataset(ahiDataset>15),cvhriTest(ahiDataset>15),'Type','Pearson')
[rhoKendall,pvalKendall] = corr(ahiDataset(ahiDataset>15),cvhriTest(ahiDataset>15),'Type','Kendall')
[rhoSpearman,pvalSpearman] = corr(ahiDataset(ahiDataset>15),cvhriTest(ahiDataset>15),'Type','Spearman')

