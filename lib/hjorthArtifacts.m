function [artifact, Messages] = hjorthArtifacts(inputSignal, fs, Setup, inputArtifacts)
% This function detects non valid signal segments and those positions are
% the outputs
%
% INPUTS:
%        inputSignal = signal vector (only one lead)
%        inputArtifacts = vector with artifact signal indexes
%        fs = sampling frequency
%        Setup.seg = time window search in seconds (Default: 5 s)
%        Setup.step = step in seconds to shift window (Default: 1 s)
%        Setup.thresholdH0Low = Lower thershold for Hjorth H0 parameter (Default: 0 Hz)
%        Setup.thresholdH0Up = Upper thershold for Hjorth H0 parameter (Default: 1 Hz)
%        Setup.thresholdH1Low = Lower threshold for Hjorth H1 parameter (Default: 1 Hz)
%        Setup.thresholdH1Up = Upper threshold for Hjorth H1 parameter (Default: 1.4 Hz)
%        Setup.thresholdH2Low = Lower thershold for Hjorth H2 parameter (Default: 0 Hz)
%        Setup.thresholdH2Up = Upper thershold for Hjorth H2 parameter (Default: 3 Hz)
%        Setup.minSegmentSeparation = Min length between noisy segments, measured in segments (Default: 5 segments)
%        Setup.plotflag = Enable, intermediate graphic results (Default false)
%        Setup.adaptive = Adaptive thresholds (Default false)
%
% OUTPUTS:
%        atifact = Position vector with non valid segments in seconds.
%                   dimensions 2 x number of non valid segment:
%                   [pos_inic_s1    pos_inic_s2 ...
%                    pos_fin_s1     pos_fin_s2  ...]

% Clean artifacts
inputSignal = inputSignal(:);
if nargin > 3
    inputSignal(logical(inputArtifacts(:))) = NaN;
end

% Initialization
Messages.errors = [];
Messages.errorsDesc = [];
Messages.warnings = [];
Messages.status = 1;

if nargin < 2
    error('Not enough arguments');
end
if ~isfield(Setup,'seg')
    Setup.seg = 5;
end
if ~isfield(Setup,'step')
    Setup.step = 1;
end
if ~isfield(Setup,'thresholdH0Low')
    Setup.thresholdH0Low = 0;
end
if ~isfield(Setup,'thresholdH0Up')
    Setup.thresholdH0Up = 1;
end
if ~isfield(Setup,'thresholdH1Low')
    Setup.thresholdH1Low = 1;
end
if ~isfield(Setup,'thresholdH1Up')
    Setup.thresholdH1Up = 1.4;
end
if ~isfield(Setup,'thresholdH2Low')
    Setup.thresholdH2Low = 0;
end
if ~isfield(Setup,'thresholdH2Up')
    Setup.thresholdH2Up = 3;
end
if ~isfield(Setup,'minSegmentSeparation')
    Setup.minSegmentSeparation = 3;
end
if ~isfield(Setup,'plotflag')
    Setup.plotflag = false;
end
if ~isfield(Setup,'adaptive')
    Setup.adaptive = false;
end


try
    % Compute derivatives
    sfilt = inputSignal-mean(inputSignal,'omitnan');
    sfilt = sfilt(:)';
    sd = [NaN diff(sfilt)];
    sdd = [NaN diff(sfilt,2) NaN];

    % Initialize Hjorth parameters
    nWindow = floor(Setup.seg*fs); % Segmentation window [samples]
    nStep = floor(Setup.step*fs);  % Sliding step [samples]
    nSegments = floor((length(inputSignal)-nWindow)/nStep);  
    H0 = zeros(1,nSegments+1);
    H1 = zeros(1,nSegments+1);
    H2 = zeros(1,nSegments+1);

    % Compute Hjorth parameters
    for i=0:nSegments
        indexes = (i*nStep)+1:(i*nStep+nWindow);
        [H0(i+1),H1(i+1),H2(i+1)] = hjorth3(sfilt(indexes),sd(indexes),sdd(indexes),fs);
    end

    % Compute thresholds
    if Setup.adaptive
        medfiltOrder = 60;
        thresholdH0Low  = medfilt1(H0,medfiltOrder,'truncate','omitnan') - Setup.thresholdH0Low;
        thresholdH0Up  = medfilt1(H0,medfiltOrder,'truncate','omitnan') + Setup.thresholdH0Up;
        thresholdH1Low = medfilt1(H1,medfiltOrder,'truncate','omitnan') - Setup.thresholdH1Low;
        thresholdH1Up = medfilt1(H1,medfiltOrder,'truncate','omitnan') + Setup.thresholdH1Up;
        thresholdH2Low  = medfilt1(H2,medfiltOrder,'truncate','omitnan') - Setup.thresholdH2Low;
        thresholdH2Up  = medfilt1(H2,medfiltOrder,'truncate','omitnan') + Setup.thresholdH2Up;
    else
        thresholdH0Low = ones(size(H0))*Setup.thresholdH0Low;
        thresholdH0Up = ones(size(H0))*Setup.thresholdH0Up;
        thresholdH1Low = ones(size(H1))*Setup.thresholdH1Low;
        thresholdH1Up = ones(size(H1))*Setup.thresholdH1Up;
        thresholdH2Low = ones(size(H2))*Setup.thresholdH2Low;
        thresholdH2Up = ones(size(H2))*Setup.thresholdH2Up;
    end
    thresholdH0Low(thresholdH0Low<0) = 0;
    thresholdH1Low(thresholdH1Low<0) = 0;
    thresholdH2Low(thresholdH2Low<0) = 0;

    % Look for low quality segments   
    artifactSegments = find((H2 > thresholdH2Up) | (H2 < thresholdH2Low) | (H1 > thresholdH1Up) | (H1 < thresholdH1Low)...
        | (H0 > thresholdH0Up) | (H0 < thresholdH0Low)); % First segment is always false due to nans in H1 and H2
    if (H0(1) > thresholdH0Up(1)) || (H0(1) < thresholdH0Low(1))
        artifactSegments = [1 artifactSegments]; 
    end

    if isempty(artifactSegments)
        artifact = [];
    else
        artifactSegments = [artifactSegments artifactSegments(end)+Setup.minSegmentSeparation]; % To use last one too
        newSegmentPosition = find(diff(artifactSegments)>= Setup.minSegmentSeparation); % New bad segment indexes
        iartifact = nan(length(newSegmentPosition),2);
        artifact = nan(length(newSegmentPosition),2);
        k = 1;
        for i = 1:length(newSegmentPosition)
           iartifact(i,1) = artifactSegments(k);
           iartifact(i,2) = artifactSegments(newSegmentPosition(i));
           k = newSegmentPosition(i)+1;
        end

        % Samples -> Seconds
        artifact(:,1) = (iartifact(:,1)-1)*nStep + 1;
        artifact(:,2) = (iartifact(:,2)-1)*nStep + nWindow;
        artifact = (artifact-1)/fs;
    end

catch me
    Messages.errors = [Messages.errors {'Error on hjorthArtifacts.m identifying non valid segments;'}];
    Messages.errorsDesc = [Messages.errorsDesc {me.message}];
    warning(char(Messages.errors(end)));
    Messages.status = 0;
end


try
    if Setup.plotflag     

        t = linspace(0,(length(inputSignal)-1)/(fs),length(H1))/60;
        tSignal = (0:1/fs:(length(inputSignal)-1)/fs)/60;
        artifact = artifact./60;
        
        figure
        ax(1) = subplot(4,1,1);
        plot(tSignal, inputSignal); hold on; grid on
        for kk = 1:size(artifact,1)
            plot(tSignal(tSignal>=artifact(kk,1) & tSignal<=artifact(kk,2)), ...
                inputSignal(tSignal>=artifact(kk,1) & tSignal<=artifact(kk,2)),'r');
        end
        ylabel('Input');
        
        ax(2) = subplot(4,1,2);
        plot(t,H0); hold on;
        plot(t,thresholdH0Up,'k','LineWidth',2);
        plot(t,thresholdH0Low,'k','LineWidth',2);
        grid on; ylabel('Activity');
        
        ax(3) = subplot(4,1,3);
        plot(t,H1); hold on;
        plot(t,thresholdH1Up,'k','LineWidth',2);
        plot(t,thresholdH1Low,'k','LineWidth',2);
        grid on; ylabel('Mobility');
        
        ax(4) = subplot(4,1,4);
        plot(t,H2); hold on;
        plot(t,thresholdH2Up,'k','LineWidth',2)
        plot(t,thresholdH2Low,'k','LineWidth',2)
        grid on; ylabel('Complexity'); xlabel('t (minutes)');
        
        linkaxes(ax,'x'); axis tight;
    end
    
catch me
    Messages.errors = [Messages.errors {'Error on hjorthArtifacts.m drawing graphical results;'}];
    Messages.errorsDesc = [Messages.errorsDesc {me.message}];
    warning(char(Messages.errors(end)));
    Messages.status = 0;
end
end


function [h0,h1,h2] = hjorth3(s,sd,sdd,fm)
    w0=(2*pi/length(s))*sum(s.^2);
    w2=(2*pi/length(s))*sum(sd.^2);
    w4=(2*pi/length(s))*sum(sdd.^2);

    h0=w0;
    h1=abs(sqrt(w2/w0))*fm/(2*pi);
    h2=abs(sqrt((w4/w2)-(w2/w0)))*fm/(2*pi);
end