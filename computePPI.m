clear
addpath(genpath('lib'));

fs = 500;
ppiFs = 4;

% Load all files in database directory
dirlist = dir('dataset/ppg_signals');
files = cell([1 length(dirlist)-2]);
for kk = 3:length(dirlist)
    files{kk-2} = dirlist(kk).name;
end

for kk = 1:length(files)

    subject = split(files{kk},'_');
    subject = subject(1);
    load(strcat('dataset/ppg_signals/',string(subject),'_pulse.mat'),'sig'); ppg = sig; clear sig
    load(strcat('dataset/ecg_signals/',string(subject),'_ecg.mat'),'sig'); ecg = -sig; clear sig
    fprintf('Computing PPI from %s...',files{kk});

    % Preprocessing
    [bb, aa] = butter(3, 0.3*2/fs, 'high');
    ppgFiltered = filtfilt(bb, aa, ppg);
    ecgFiltered = filtfilt(bb, aa, ecg);
    ecgFiltered = normalize(ecgFiltered);
    [bb, aa] = butter(3, 30*2/fs, 'low');
    ppgFiltered = filtfilt(bb, aa, ppgFiltered);
    ppgFiltered = normalize(ppgFiltered);
    clear bb aa

    % Delineation
    % Maximum upslope detection
    DelineatorConfig.alfa = 0.1; % Percentage of threshold max value used to its low limit
    DelineatorConfig.refract = 250e-3; % Refractory interval in wich threshold remains constant
    DelineatorConfig.taoRR = 1; % Fraction of RR estimated interval where threshold falls to its min malue
    DelineatorConfig.w_nA = 700e-3; % Window size for searching the maximum (time from the derivate max)
    DelineatorConfig.plotflag = false;
    [pulses] = pulseDelineation(ppgFiltered',fs,DelineatorConfig); %set(gcf,'position',[0,0,1000,500])

    tn = gapcorrectorNonLinear(pulses,false); tn(isnan(tn)) = []; tn = unique(tn);
    tPpi = tn(2):1/ppiFs:tn(end);
    ppi = interp1(tn(2:end),diff(tn),tPpi,'linear');
    ppi = smooth(tPpi,ppi,20*ppiFs,'loess'); % 2-order polynomial fitting
    
    subject = split(files{kk},'_');
    subject = subject(1);
    save(strcat('results/ppi/',string(subject),'_ppi.mat'),'ppi','ppiFs','tPpi');

    clear pulses ppgFiltered tn tPpg
    fprintf('Done\n')
end