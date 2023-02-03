function [ nD , nA , nB , nM , signalLPD , threshold ] = pulseDelineation ( signal , fs , Setup )
% Plethysmography signals delineation, based on low-pass differentiator filter.
% [ nD , nA , nB , nM , signalLPD , threshold ] = pulseDelineation ( signal , fs , Setup )
%
%
% <new>
%       For faster computation in the delineation, this function segmentizes the signal in 100000 sps
%       (independent to fs). In those segments, the general adaptive thresholding is computed.
%       Border effects are corrected
%
%       With the previous verson, the Computation Time grew exponentially with the length of the signal.
%           24 hours of recordings were imposible to be delineated.
%
%       With this version, the length of the signal is almost not influencing the Computation Time
%
%
%--------------------------------------------------------
%   In:
%         ppg           = PPG signal
%         fs            = sampling rate (Hz)
%         isArtifact	= PPG artifacts [Default: []]
%         fpLPD        = last frequency of the pass-band for the
%                           low-pass-differentiator filter (Hz) [Default: 7.8]
%         fcLPD        = cut-off frequency for the low-pass-differentiator filter
%                           (Hz) [Default: 8]
%                           minimum value (alfa*amplitude of previous SSF peak)
%                           [Default: 1]
%                           If tauRR increases, steeper slope
%         thrIncidences= threshold for incidences
%         plotflag      = if true, plots a figure with PPG and SSF [Default: false]
%
%--------------------------------------------------------
%
%         orderLPD     = filter order [Default: 3*fsppg]
%         alfa          = multiplies previus amplitude of detected maximum in
%                           filtered signal for updating the threshold [Default: 0.2]
%         refract       = refractary period for threshold (s) [Default: 150e-3]
%         tauRR         = fraction of estimated RR where threshold reaches its
%   Out:
%         nD            =	location of peaks detected in filtered signal (samples)
%         signalLPD = band-pass filtered signal
%         thres = computed time varying theshold
%
%--------------------------------------------------------
%   Usage of isArtifact
%           PPG         =	[ 0    1   2  3   2   1   0 ]
%           isArtifact	=	[ 0	   1   1  0   0   0   0 ]
%
%       result:
%           PPG         =	[ 0  NaN NaN  3   2   1   0 ]
%
%--------------------------------------------------------
%
% % % % % % % % % % % % % % % % % % % % % % % % %
% LO QUE MAS LE CUESTA CALCULAR A LA FUNCION ES EL LPD_FILTER (75% del t_computo)
% SI VAMOS A DELINEAR UNA BdD CON LOS MISMOS PARAMETROS, MEJOR METER EL FILTRO COMO INPUT
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Created                by Jesús Lázaro <jlazarop@unizar.es> in 2018, 
% Optimized and extended by Pablo Armañac in 2019-2020
%
%

%% Check Inputs
if nargin <= 2,   Setup = struct();                       end
if nargin <  2,   error('Not enough input arguments.');   end

% Default Values
if ~isfield(Setup,'isArtifact'),                Setup.isArtifact = false(length(signal),1);     end

if ~isfield(Setup,'filterLPD'),                 Setup.filterLPD = [];                           end
if ~isfield(Setup,'orderLPD'),                  Setup.orderLPD = 3*fs;                          end
if ~isfield(Setup,'fpLPD'),                     Setup.fpLPD = 7.8;                              end
if ~isfield(Setup,'fcLPD'),                     Setup.fcLPD = 8;                                end

if ~isfield(Setup,'nD'),                        Setup.nD =	[];                                 end

if ~isfield(Setup,'Lenvelope'),                 Setup.Lenvelope = (fs*300)/1000;                end

if ~isfield(Setup,'alfa'),                      Setup.alfa = 0.2;                               end
if ~isfield(Setup,'tauRR'),                     Setup.tauRR = 1;                                end
if ~isfield(Setup,'refractPeriod'),             Setup.refractPeriod = 150e-03;                  end
if ~isfield(Setup,'thrIncidences'),             Setup.thrIncidences = 1.5;                      end

if ~isfield(Setup,'wdw_nA'),                    Setup.wdw_nA = 250e-3;                          end
if ~isfield(Setup,'wdw_nB'),                    Setup.wdw_nB = 150e-3;                          end
if ~isfield(Setup,'fsi'),                       Setup.fsi = 2*fs;                                 end

if ~isfield(Setup,'computeLPDFiltering'),       Setup.computeLPDFiltering     =	true;           end
if ~isfield(Setup,'computeAdaptiveThreshold'),  Setup.computeAdaptiveThreshold = true;          end
if ~isfield(Setup,'computeEnvelopesThreshold'), Setup.computeEnvelopesThreshold = false;       end
if ~isfield(Setup,'computeEyeCorrection'),      Setup.computeEyeCorrection        =	false;      end
if ~isfield(Setup,'computePeakDelineation'),	Setup.computePeakDelineation	=	true;       end

if ~isfield(Setup,'plotflag'),                  Setup.plotflag                  =	false;      end

% Get, assign and store the variable names
data_names = fieldnames(Setup);
for ii = 1:length(data_names), eval([ data_names{ii} ' = Setup.' data_names{ii} ';']); end
clear Setup data_names ii


%% Remove artifacts
isArtifact          =	logical(isArtifact); %#ok
PPG_orig            =	signal - nanmean(signal(~isArtifact));
signal(isArtifact)	=	NaN;
signal              =	signal(:) - nanmean(signal(:)); %Ensure "ppg" is a column vector without trend


%% Compute LPD filtering
signalLPD = signal; % in case you want to detect directly in the raw signal
if computeLPDFiltering
    
    LPDFilterSetup.filterLPD	=	filterLPD;
    LPDFilterSetup.orderLPD     =	orderLPD;
    LPDFilterSetup.fpLPD        =	fpLPD;
    LPDFilterSetup.fcLPD        =	fcLPD;
    
    signalLPD                   =	LPDFiltering( signal , fs , LPDFilterSetup );
    %     signalLPD_orig              =	LPDFiltering( PPG_orig , fs , LPDFilterSetup );
end


%% Compute threshold and nD detection
threshold = NaN(length(signal),1);
if isempty (nD) %#ok
    
    if computeEnvelopesThreshold
        
        delineatorSetup.Lenvelope	=	Lenvelope;
        
        [ nD , threshold ]                      =	envelopeThresholdDelineator ( signalLPD, fs, delineatorSetup );
        
        
    elseif computeAdaptiveThreshold
        
        delineatorSetup.alfa             =	alfa;
        delineatorSetup.tauRR            =	tauRR;
        delineatorSetup.refractPeriod	=	refractPeriod;
        delineatorSetup.thrIncidences	=	thrIncidences;
        
        [ nD , threshold ]                      =	adaptiveThresholdDelineator ( signalLPD, fs, delineatorSetup );
        
    end
    
end


%% Correction by Eye
if computeEyeCorrection
    nD = 1+round(nD*fs);
    [nD, nD_orig] = tk_inspect_PPG (  PPG_orig, nD, fs, signalLPD_orig ); %#ok
    nD = (nD-1)./fs;
end


%% Compute delineation
nA = []; nB = []; nM = [];
if computePeakDelineation
    
    peakSetup.wdw_nB        =   wdw_nB;
    peakSetup.wdw_nA        =   wdw_nA;
    peakSetup.fsi           =   fsi;
    peakSetup.diffSignal    =   signalLPD;
    
    [ nA , nB , nM , nD ] = fastPeakDelineation ( signal , fs , nD , peakSetup ) ;
%     [ nA , nB , nM , nD ] = peakDelineation ( signal , fs , nD , peakSetup ) ;
end


%% Plot
if plotflag
    
    t = 0:1/fs:(length(signal)-1)/fs;
    
    figure;
    
    ax(1) = subplot(2,1,1); hold on;
    plot(t,  PPG_orig, 'Color',[0.65 0.65 0.65] );
    plot(t,  signal, 'k');
    plot(nD(~isnan(nD)), signal(1+round(nD(~isnan(nD))*fs)), 'ro');
    plot(nA(~isnan(nA)), signal(1+round(nA(~isnan(nA))*fs)), 'b*');
    plot(nM(~isnan(nM)), signal(1+round(nM(~isnan(nM))*fs)), 'k*');
    plot(nB(~isnan(nB)), signal(1+round(nB(~isnan(nB))*fs)), 'b*');
    title('PPG');
    
    ax(2) = subplot(2,1,2); hold on;
%     plot(t, signalLPD_orig, 'Color',[0.65 0.65 0.65]);
    plot(t, signalLPD, 'k');
    plot(t, threshold, ':');
    plot(nD(~isnan(nD)), signalLPD(1+round(nD(~isnan(nD))*fs)), 'ro');
    title('LPD Filtered PPG');
    
    linkaxes(ax, 'x');
    zoom on
    
end


end

