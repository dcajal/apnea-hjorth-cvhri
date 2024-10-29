function [ parameters, segments ] = computeHjorthParametersAndHRV( spo2,spo2Fs,ppi,ppiFs,ppiNoFilter )
    Setup.hrv = false;
    [ppiHjorth,segments] = hjorth(ppi, ppiFs, Setup);
    [spo2Hjorth,~] = hjorth(spo2, spo2Fs, Setup);
    
    Setup.hrv = true;
    [prv,~] = hjorth(ppiNoFilter, ppiFs, Setup);
    
    spo2Hjorth = spo2Hjorth(1:size(segments,1),:);
    
    % Output
    parameters = [spo2Hjorth ppiHjorth prv];
    %parameters = [spo2Hjorth prv];
    %parameters = [ppiHjorth prv];
    end
    
    function [parameters, segments] = hjorth(inputSignal, fs, Setup)
    % Compute Hjorth parameters
    %
    % INPUTS:
    %        inputSignal = signal vector
    %        fs = sampling frequency
    %        Setup.seg = time window search in seconds (Default: 240 s)
    %        Setup.step = step in seconds to shift window (Default: 30 s)
    %
    % OUTPUTS:
    %        parameters = Activity, mobility, complexity in columns
    
    
    % Check inputs and initialize defaults
    
    if nargin < 3
        Setup =  [];
    end
    if ~isfield(Setup,'seg')
        Setup.seg = 180;
    end
    if ~isfield(Setup,'step')
        Setup.step = 30;
    end
    if ~isfield(Setup,'hrv')
        Setup.hrv = false;
    end
    
    nWindow = floor(Setup.seg*fs); % Segmentation window [samples]
    nStep = floor(Setup.step*fs);  % Sliding step [samples]
    nSegments = floor((length(inputSignal)-nWindow)/nStep);  
    segments = zeros(nSegments+1,2);
    
    if Setup.hrv
        % Initialize HRV arrays
        lf = zeros(1,nSegments+1);
        hf = zeros(1,nSegments+1);
        lfhf = zeros(1,nSegments+1);
    
        % Compute HRV
        for kk=0:nSegments
            indexes = (kk*nStep)+1:(kk*nStep+nWindow);
            [lf(kk+1),hf(kk+1),lfhf(kk+1)] = computeHRV(inputSignal(indexes),fs);
            segments(kk+1,:) = ([indexes(1) indexes(end)]-1)/fs;
        end
    
        parameters = [lf; hf; lfhf]';
    else
        % Compute derivatives
        sfilt = detrend(inputSignal,0,1:fs*Setup.seg:length(inputSignal)); % Removes mean value segmentwise
        sfilt = sfilt(:)';
        sd = [NaN diff(sfilt)];
        sdd = [NaN diff(sfilt,2) NaN];
        
        % Initialize Hjorth parameters
        h0 = zeros(1,nSegments+1);
        h1 = zeros(1,nSegments+1);
        h2 = zeros(1,nSegments+1);
    
        % Compute Hjorth parameters
        for kk=0:nSegments
            indexes = (kk*nStep)+1:(kk*nStep+nWindow);
            [h0(kk+1),h1(kk+1),h2(kk+1)] = computeHjorth(sfilt(indexes),sd(indexes),sdd(indexes),fs);
            segments(kk+1,:) = ([indexes(1) indexes(end)]-1)/fs;
        end
        
        parameters = [h0; h1; h2]';
    end
end

function [lf,hf,lfhf] = computeHRV(ppi,fs)
    f = linspace(0,0.4,2^10);
    pxx = periodogram(ppi,hamming(length(ppi)),f,fs);
    [lf,hf,lfhf] = freqind(pxx,f);
end

function [h0,h1,h2] = computeHjorth(s,sd,sdd,fm)
    w0=(2*pi/length(s))*sum(s.^2);
    w2=(2*pi/length(s))*sum(sd(~isnan(sd)).^2);
    w4=(2*pi/length(s))*sum(sdd(~isnan(sdd)).^2);

    h0=w0;
    h1=abs(sqrt(w2/w0))*fm/(2*pi);
    h2=abs(sqrt((w4/w2)-(w2/w0)))*fm/(2*pi);
end