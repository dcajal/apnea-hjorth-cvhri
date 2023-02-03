function [ nD, threshold ] = adaptiveThresholdDelineator( signal, fs, Setup )
% adaptiveThresholdDelineator
%   find local maxima using adaptive thresholding
%
% <new>
%   For faster computation, this function segmentizes the signal in 100000 sps (independent to fs)
%   In those segments, the general adaptive thresholding is computed.
%   Border effects are corrected
%
% --------------------------------------------------------------------------------
%   Sintax: [ nD, threshold ] = adaptiveThresholdDelineator( signal, fs, Setup )
%   In:   signal	=	input signal (PPG or BP)
%         fs        =	sampling rate [Hz]
%         alfa      =	multiplies previus amplitude of detected maximum in
%                           filtered signal for updating the threshold (Default 0.2)
%         refractPeriod	= refractory period for threshold [seconds] (Default: 150e-03)
%         tauRR     =	fraction of estimated RR where threshold reaches its minimum value
%                           (alfa*amplitude of previous SSF peak) (Default: 1)
%         thrIncidences= threshold for incidences
%
%   Out:  nD        = location of peaks detected in filtered signal [seconds]
%         threshold	= computed time varying theshold
%
% Created   by Jesús Lázaro  <jlazarop@unizar.es> in 2018
% Optimized by Pablo Armañac <parmanac@unizar.es> in 2019
%

%% Check Inputs
if nargin <= 2,   Setup = struct();                       end
if nargin <  2,   error('Not enough input arguments.');   end

% Default Values
if ~isfield(Setup,'alfa'),                      Setup.alfa = 0.2;               end
if ~isfield(Setup,'tauRR'),                     Setup.tauRR = 1;                end
if ~isfield(Setup,'refractPeriod'),             Setup.refractPeriod = 150e-03;	end
if ~isfield(Setup,'thrIncidences'),             Setup.thrIncidences = 1.5;      end

% Get, assign and store the variable names
data_names = fieldnames(Setup);
for ii = 1:length(data_names), eval([ data_names{ii} ' = Setup.' data_names{ii} ';']); end
clear Setup data_names ii

refractPeriod	=	round(refractPeriod*fs); %#ok

%% Segmentize the signals (this is JUST for FASTER COMPUTATION)

signal_length = length(signal);

% Set the segment size to 100000 SPS and compute the amount of minutes in the signal
segmentLength = 100000; % in SPS.	If fs = 1000 Hz -> 1.5  mins
%                                   If fs = 200  Hz -> 8.33 mins

nSegmentsOriginal = floor(signal_length/segmentLength);

% Add zeros if the signal does not contain a round number of minutes
if nSegmentsOriginal > 1
    newSignal = [signal; NaN(segmentLength-(signal_length-nSegmentsOriginal*segmentLength),1)];
    % Get the new length and amount of segments
    newSignalLength = length(newSignal);
    
    nSegments = floor(newSignalLength/segmentLength);
    segmentedSignal = (double(reshape(newSignal,segmentLength,nSegments)));
else
    nSegments = 1;
    segmentedSignal = signal;
end

nD          =	[];
threshold	=	[];

% Go through each segment
for ii = 1:nSegments
    % Make the segments tAdd seconds longer on each side
    
    tAdd = 5*fs;
    
    if nSegments > 1
        if ii == 1
            % End the segment 10 seconds later
            segment = [ segmentedSignal(:,1); segmentedSignal(1:tAdd,2)];
        elseif ii==nSegments
            % Start the segment 10 seconds earlier
            segment = [ segmentedSignal(end-(tAdd)+1:end,end-1); segmentedSignal(:,end)];
        else
            % Start the segement 10 seconds earlier and end 10 seconds later
            segment = [ segmentedSignal(end-(tAdd)+1:end,ii-1);segmentedSignal(:,ii);segmentedSignal(1:tAdd,ii+1)];
        end
    else
        % Take the whole segment
        segment = segmentedSignal;
    end
    
    time = (0:1:length(segment)-1)./fs;
    
    %% Detect peaks in the LPD signal by adaptive thresholding
    
    [ nD_segment , thresholdSegment ] = adaptiveThresholding( segment(:) , fs , alfa , refractPeriod , tauRR , thrIncidences );
    
%     figure;
%     plot(time,segment);  hold on;
%     plot(time(nD_segment),segment(nD_segment),'o');
%     plot(time,thresholdSegment);
%     plot(time(nD_segment),thresholdSegment(nD_segment),'o');
    
    %% Remove added signal on both sides
    %     if length(nD_segment) > time(end)/2
    if nSegments > 1
        if ii == 1
            % Remove the last five seconds
            nD_segment = nD_segment(nD_segment<=segmentLength);
            thresholdSegment = thresholdSegment(time*fs < segmentLength);
        elseif ii == nSegments
            % Remove the first five seconds
            nD_segment = nD_segment(nD_segment>tAdd)-(tAdd);
            thresholdSegment = thresholdSegment(time*fs >= tAdd);
        else
            % Remove the first and last five seconds
            nD_segment = nD_segment(nD_segment>=tAdd)-tAdd;
            nD_segment = nD_segment(nD_segment<segmentLength);
            
            thresholdSegment = thresholdSegment(time*fs >= tAdd);
            thresholdSegment = thresholdSegment(time*fs < segmentLength);
        end
    end
    
    % Store the signals in a cell
    nD          =	[ nD;           nD_segment(:)+segmentLength*(ii-1)	]; %#ok
    threshold	=	[ threshold;	thresholdSegment(:)                    ]; %#ok
    
end

%% Arrange Outputs
t = (0:1:length(signal)-1)./fs;
nD = ( unique(nD)-1 ) ./ fs;
threshold (length(t)+1:end) = [];

end





function [ nD , thres ] = adaptiveThresholding( sig_filt , fs , alfa , refractPeriod , tauRR , thrIncidences )

nD = [];
peaks_filt_orig = [];
peaks_added = [];
cond_vec = [];
thres_ini_w_ini = find(~isnan(sig_filt), 1, 'first');
thres_ini_w_end = thres_ini_w_ini + round(10*fs); thres_ini_w_end(thres_ini_w_end>=length(sig_filt)) = [];
aux = sig_filt(thres_ini_w_ini:thres_ini_w_end);
thres_ini = 3*nanmean(aux(aux>=0));
thres = nan(size(sig_filt));
t = 1:length(sig_filt);
RR = round(60/80*fs);

if (1+RR)<length(sig_filt)
    thres(1:1+RR) = thres_ini - (thres_ini*(1-alfa)/RR)*(t(1:RR+1)-1);
    thres(1+RR:end) = alfa*thres_ini;
else
    thres(1:end) = thres_ini - (thres_ini*(1-alfa)/RR)*(t(1:end)-1);
end

kk=1;
while true
    cross_u = kk-1 + find(sig_filt(kk:end)>thres(kk:end), 1, 'first'); %Next point to cross the actual threshold (down->up)
    if isempty(cross_u)
        % No more pulses -> end
        break;
    end
    
    cross_d = cross_u-1 + find(sig_filt(cross_u:end)<thres(cross_u:end), 1, 'first'); %Next point to cross the actual threshold (up->down)
    
    if isempty(cross_d)
        % No more pulses -> end
        break;
    end
    
    % Pulse detected:
    [vmax, imax] = max(sig_filt(cross_u:cross_d));
    p = cross_u-1+imax;
    nD = [nD, p];
    peaks_filt_orig = nD;
    Npeaks = length(nD);
    if Npeaks>3
        tk_c = peaks_filt_orig(end);
        tk1_c = peaks_filt_orig(end-1);
        tk2_c = peaks_filt_orig(end-2);
        tk3_c = peaks_filt_orig(end-3);
        cond = abs((2*tk1_c - tk2_c - tk_c) / ((tk1_c-tk2_c)*(tk_c-tk1_c)*(tk_c-tk2_c)));
        cond_vec = [cond_vec cond];
        if cond >=thrIncidences/(fs*fs)
            tk = nD(end);
            tk1 = nD(end-1);
            tk2 = nD(end-2);
            tk3 = nD(end-3);
            
            %Inserting a beat between tk2 and tk1:
            aux_15 = sig_filt(tk2:tk1);
            [aux_peaks, aux_locs] = findpeaks(aux_15);
            if ~isempty(aux_locs)
                aux_locs = aux_locs(aux_peaks>=0.5*max(aux_peaks));
                [~, aux_loc] = min(abs(aux_locs-(length(aux_15))/2));
                tk15 = tk2 - 1 + aux_locs(aux_loc);
            else
                tk15 = nan;
            end
            %Inserting a beat between tk1 and tk:
            aux_05 = sig_filt(tk1:tk);
            [aux_peaks, aux_locs] = findpeaks(aux_05);
            if ~isempty(aux_locs)
                aux_locs = aux_locs(aux_peaks>=0.5*max(aux_peaks));
                [~, aux_loc] = min(abs(aux_locs-(length(aux_05))/2));
                tk05 = tk1 - 1 + aux_locs(aux_loc);
            else
                tk05 = nan;
            end
            
            %Condition removing previous detection (cond2)
            cond1 = abs((2*tk2 - tk3 - tk1) / ((tk2-tk3)*(tk1-tk2)*(tk1-tk3)));
            
            %Condition removing previous detection (cond2)
            cond2 = abs((2*tk2 - tk3 - tk) / ((tk2-tk3)*(tk-tk2)*(tk-tk3)));
            
            %Condition adding a new detection between tk2 and tk1 (cond3)
            if ~isnan(tk15)
                cond3 = abs((2*tk1 - tk15 - tk) / ((tk1-tk15)*(tk-tk1)*(tk-tk15)));
            else
                cond3 = inf;
            end
            
            %Condition adding a new detection between tk1 and tk (cond4)
            if ~isnan(tk05)
                cond4 = abs((2*tk05 - tk1 - tk) / ((tk05-tk1)*(tk-tk05)*(tk-tk1)));
            else
                cond4 = inf;
            end
            
            
            [~, high_cond] = min([cond1, cond2, cond3, cond4]);
            
            switch high_cond
                case 1 %Best is to remove current detection
                    nD = nD(1:end-1);
                    cond_vec = cond_vec(1:end-1);
                    Npeaks = Npeaks-1;
                    kk = cross_d;
                    continue;
                case 2 %Best is to remove previous detection
                    [vmax, imax] = max(sig_filt(cross_u-refractPeriod:cross_d));
                    if imax~=1
                        p = cross_u-refractPeriod-1+imax;
                    end
                    nD = [nD(1:end-2) p];
                    cond_vec = cond_vec(1:end-1);
                    Npeaks = Npeaks-1;
                case 3 % Best is to add a detection between tk2 and tk1
                    %                         peaks_filt = [peaks_filt(1:end-2) tk15 peaks_filt(end-1:end)];
                    %                         cond_vec = [cond_vec(1:end-2) cond3 cond_vec(end-1:end)];
                    %                         Npeaks = Npeaks+1;
                    peaks_added = [peaks_added tk15];
                case 4 % Best is to add a detection between tk1 and tk
                    %                         peaks_filt = [peaks_filt(1:end-1) tk05 peaks_filt(end)];
                    %                         cond_vec = [cond_vec(1:end-1) cond4 cond_vec(end)];
                    %                         Npeaks = Npeaks+1;
                    peaks_added = [peaks_added tk05];
            end
        end
    end
    
    % Update threshold
    N_RR_estimation = 3;
    N_ampli_est = 3;
    if Npeaks>=N_RR_estimation+1
        RR = round(median(diff(nD(end-N_RR_estimation:end))));
    elseif Npeaks>=2
        RR = round(mean(diff(nD)));
    end
    kk = min(p+refractPeriod, length(sig_filt));
    thres(p:kk) = vmax;
    
    vfall = vmax*alfa;
    if Npeaks>=(N_ampli_est+1)
        ampli_est = median(sig_filt(nD(end-N_ampli_est:end-1)));
        if vmax>=(2*ampli_est)
            vfall = alfa*ampli_est;
            vmax = ampli_est;
        end
    end
    
    fall_end = round(tauRR*RR);
    if (kk+fall_end)<length(sig_filt)
        thres(kk:kk+fall_end) = vmax - (vmax-vfall)/fall_end*(t(kk:kk+fall_end)-kk);
        thres(kk+fall_end:end) = vfall;
    else
        thres(kk:end) = vmax - (vmax-vfall)/fall_end*(t(kk:end)-kk);
    end
    
end

nD = unique([nD peaks_added]);

end