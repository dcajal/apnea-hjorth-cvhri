% Visual inspection for new annotations

clear
addpath(genpath('lib'));
addpath(genpath('dataset'));

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaulttextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize', 20);
set(0,'defaultTextFontSize', 20);
set(0, 'defaultFigureRenderer', 'painters')


% Load signals
subjectNumber = '02'; % CVHR events: 2, 5, 6, 8, 10...

load(['results/labels/UZLeuven0' subjectNumber '_newlabels.mat'])
load(['results/signals/UZLeuven0' subjectNumber '_psg.mat'])
load(['results/ppi/UZLeuven0' subjectNumber '_ppi.mat'])


% Plots
figure;
ax(1) = subplot(511);
hold on
for kk=1:size(apneas,1)
    p(1) = patch([apneas(kk, 1) apneas(kk, 2) apneas(kk, 2) apneas(kk, 1)],[-6 -6 15 15],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(hypopneas,1)
    p(2) = patch([hypopneas(kk, 1) hypopneas(kk, 2) hypopneas(kk, 2) hypopneas(kk, 1)],[-6 -6 15 15],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(doubts,1)
    p(3) = patch([doubts(kk, 1) doubts(kk, 2) doubts(kk, 2) doubts(kk, 1)],[-6 -6 15 15],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
plot(tNasalPressure,nasalPressureProcessed)
plot(tNasalPressure,medfilt1(tidalVolume,180*nasalPressureFs)+5,'k--');
plot(tNasalPressure,tidalVolume+5);
yline(5)
ylabel('Nasal pressure')
axis tight;

ax(2) = subplot(512);
plot(tSpo2,spo2Processed,'b'); axis tight;
hold on
for kk=1:size(apneas,1)
    p(1) = patch([apneas(kk, 1) apneas(kk, 2) apneas(kk, 2) apneas(kk, 1)],[80 80 100 100],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(hypopneas,1)
    p(2) = patch([hypopneas(kk, 1) hypopneas(kk, 2) hypopneas(kk, 2) hypopneas(kk, 1)],[80 80 100 100],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(doubts,1)
    p(3) = patch([doubts(kk, 1) doubts(kk, 2) doubts(kk, 2) doubts(kk, 1)],[80 80 100 100],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
ylabel('SpO2')
% legend(p,'Apnea','Hypopnea','Doubts','Location','southeast'); clear p

ax(3) = subplot(513);
plot(tPpi,ppi,'b'); axis tight;
hold on
for kk=1:size(apneas,1)
    p(1) = patch([apneas(kk, 1) apneas(kk, 2) apneas(kk, 2) apneas(kk, 1)],[0.6 0.6 1.4 1.4],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(hypopneas,1)
    p(2) = patch([hypopneas(kk, 1) hypopneas(kk, 2) hypopneas(kk, 2) hypopneas(kk, 1)],[0.6 0.6 1.4 1.4],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(doubts,1)
    p(3) = patch([doubts(kk, 1) doubts(kk, 2) doubts(kk, 2) doubts(kk, 1)],[0.6 0.6 1.4 1.4],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
ylabel('PPI')

ax(4) = subplot(514);
plot(tAbdBelt,ripSum,'b'); axis tight;
hold on
for kk=1:size(apneas,1)
    p(1) = patch([apneas(kk, 1) apneas(kk, 2) apneas(kk, 2) apneas(kk, 1)],[-10 -10 10 10],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(hypopneas,1)
    p(2) = patch([hypopneas(kk, 1) hypopneas(kk, 2) hypopneas(kk, 2) hypopneas(kk, 1)],[-10 -10 10 10],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(doubts,1)
    p(3) = patch([doubts(kk, 1) doubts(kk, 2) doubts(kk, 2) doubts(kk, 1)],[-10 -10 10 10],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
ylabel('RIPsum')

ax(5) = subplot(515);
p = plot(tHypno,hypno,'b'); axis tight;
hold on
for kk=1:size(apneas,1)
    p(1) = patch([apneas(kk, 1) apneas(kk, 2) apneas(kk, 2) apneas(kk, 1)],[0 0 6 6],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(hypopneas,1)
    p(2) = patch([hypopneas(kk, 1) hypopneas(kk, 2) hypopneas(kk, 2) hypopneas(kk, 1)],[0 0 6 6],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(doubts,1)
    p(3) = patch([doubts(kk, 1) doubts(kk, 2) doubts(kk, 2) doubts(kk, 1)],[0 0 6 6],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
xlabel('Time (seconds)');
ylabel('Hypnogram')
yticks(1:5);
yticklabels({'NREM3','NREM2','NREM1','REM','WAKE'})

linkaxes(ax,'x'); clear ax



figure;
ax(1) = subplot(211);
hold on
for kk=1:size(apneas,1)
    p(1) = patch([apneas(kk, 1) apneas(kk, 2) apneas(kk, 2) apneas(kk, 1)],[-6 -6 15 15],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(hypopneas,1)
    p(2) = patch([hypopneas(kk, 1) hypopneas(kk, 2) hypopneas(kk, 2) hypopneas(kk, 1)],[-6 -6 15 15],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(doubts,1)
    p(3) = patch([doubts(kk, 1) doubts(kk, 2) doubts(kk, 2) doubts(kk, 1)],[-6 -6 15 15],[0 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
plot(tNasalPressure,nasalPressureProcessed)
plot(tNasalPressure,medfilt1(tidalVolume,180*nasalPressureFs)+5,'k--');
plot(tNasalPressure,tidalVolume+5);
yline(5)
ylabel('Event-based labels')
axis tight;

ax(2) = subplot(212);
hold on
for kk=1:size(hypoSegments,1)
    p(1) = patch([hypoSegments(kk, 1) hypoSegments(kk, 2) hypoSegments(kk, 2) hypoSegments(kk, 1)], ...
        [-6 -6 15 15],[0 0 1],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
for kk=1:size(apneaSegments,1)
    p(2) = patch([apneaSegments(kk, 1) apneaSegments(kk, 2) apneaSegments(kk, 2) apneaSegments(kk, 1)], ...
        [-6 -6 15 15],[1 0 0],'FaceAlpha',.3,'EdgeColor','none');
end; clear kk
plot(tNasalPressure,nasalPressureProcessed)
plot(tNasalPressure,medfilt1(tidalVolume,180*nasalPressureFs)+5,'k--');
plot(tNasalPressure,tidalVolume+5);
yline(5)
ylabel('Segment-based labels')
axis tight;

linkaxes(ax,'x'); clear ax
