% Process PSG signals. Resampling, filtering, tidal volume..

clear
addpath(genpath('lib'));

% Load all files in database directory
dirlist = dir('dataset/nasal_pressure_signals/');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

fso = 500; % Original sample rate

for kk = 1:length(files)

    subject = split(files{kk},'_');
    subject = subject(1);
    fprintf('Computing subject: %s...',string(subject));

    % SpO2
    load(strcat('dataset/spo2_signals/',string(subject),'_sao2.mat'),'sig'); spo2 = sig; clear sig
    spo2Fs = 25; % Fs after decimate (AASM Recommended: 25 Hz)
    r = fso/spo2Fs; spo2 = decimate(spo2,r);
    tSpo2 = 0:1/spo2Fs:(length(spo2)-1)/spo2Fs;
    spo2Processed = spo2Processing(spo2,spo2Fs);

    % Nasal pressure
    load(strcat('dataset/nasal_pressure_signals/',string(subject),'_oronas_pres.mat'),'sig'); nasalPressure = sig; clear sig
    nasalPressureFs = 100; % Fs after decimate (AASM Recommended: 100 Hz)
    r = fso/nasalPressureFs; nasalPressure = decimate(nasalPressure,r);
    tNasalPressure = 0:1/nasalPressureFs:(length(nasalPressure)-1)/nasalPressureFs;
    [nasalPressureProcessed, tidalVolume, nasalPressureUpperEnvelope, nasalPressureLowerEnvelope, nasalPressureArtifacts] = ...
        nasalPressureProcessing(nasalPressure,nasalPressureFs);
    
    % Abdominal belt
    load(strcat('dataset/abd_belt_signals/',string(subject),'_resp_abd.mat'),'sig'); abdBelt = sig; clear sig
    abdBeltFs = 100; % Fs after decimate (AASM Recommended: 100 Hz)
    r = fso/abdBeltFs; abdBelt = decimate(abdBelt,r);
    tAbdBelt = 0:1/abdBeltFs:(length(abdBelt)-1)/abdBeltFs;
    
    % Thoracic belt
    load(strcat('dataset/thor_belt_signals/',string(subject),'_resp_thor.mat'),'sig'); thorBelt = sig; clear sig
    thorBeltFs = 100; % Fs after decimate (AASM Recommended: 100 Hz)
    r = fso/thorBeltFs; thorBelt = decimate(thorBelt,r);
    tThorBelt = 0:1/thorBeltFs:(length(thorBelt)-1)/thorBeltFs;

    % RIPsum
    ripSum = normalize(abdBelt+thorBelt);
    
    % Hypno
    load(strcat('dataset/annotations/hypno/',string(subject),'.mat')); hypno = cell2mat(hypno);
    tHypno = linspace(0,tNasalPressure(end),numel(hypno));


    save(strcat('results/signals/',string(subject),'_psg.mat'),...
        'abdBelt','abdBeltFs','thorBelt','thorBeltFs','hypno','nasalPressureProcessed','nasalPressureFs','nasalPressureLowerEnvelope', ...
        'nasalPressureUpperEnvelope','tidalVolume','ripSum','spo2Processed','spo2Fs','tAbdBelt','tThorBelt','tNasalPressure',"tHypno" ...
        ,"tSpo2",'nasalPressureArtifacts');

    fprintf('Done\n');

end