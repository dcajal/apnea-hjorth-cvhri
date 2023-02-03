function varargout = fastPeakDelineation ( signal , fs , i_nD , Setup )
%
% Peak delineation for plethysmography signals, given nD_in (=nD) as anchor point
%
% Inputs
%         ppg           signal
%         fs            sampling frequency [Hertz]
%         peaks_filt	detections in the maximum of the first derivative of the PPG [seconds] (detected at fs)
%         wdw_nB        window width for searching the minimum before nD [seconds] (default = 150e-03)
%         wdw_nA        window width for searching the maximum after  nD [seconds] (default = 200e-03)
%         fsi           sampling frequency for interpolation. Fine search of the peaks [Hertz]
%         diffPPG     derivative of the PPG signal for searching nD
%
%
% Outputs all (detected at fsi) [seconds]
%         nA            Maximum of the PPG pulse
%         nB            Minimum of the PPG pulse
%         nM            Medium  of the PPG pulse
%         nD            Maximum of the first derivative of the PPG
%
%
% Created               by Jes�s L�zaro  <jlazarop@unizar.es> in 2014
% Fixed and Optimized   by Pablo Arma�ac <parmanac@unizar.es> in 2019
%
% Esto se tiene que optimizar.
%       1) findpeaks.
%       2) interp1 lo hace or vectores, pues se hace directamente en vez de en el for
%
%

%% Check Inputs
if nargin <= 3,   Setup = struct();                       end
if nargin <  3,   error('Not enough input arguments.');   end

% Default Values
if ~isfield(Setup,'wdw_nA'),	Setup.wdw_nA    = 250e-3;	end
if ~isfield(Setup,'wdw_nB'),	Setup.wdw_nB    = 150e-3;	end
if ~isfield(Setup,'fsi'),       Setup.fsi       = fs;       end
if ~isfield(Setup,'diffSignal'),Setup.diffSignal = diff(fillmissing(signal(:),'linear'));	end
if ~isfield(Setup,'plotflag'),  Setup.plotflag  = false;	end

% Get, assign and store the variable names
data_names = fieldnames(Setup);
for ii = 1:length(data_names), eval([ data_names{ii} ' = Setup.' data_names{ii} ';']); end
clear Setup data_names ii bb aa

%% Check inputs
if isempty(i_nD)
    varargout{1} = NaN;
    varargout{2} = NaN;
    varargout{3} = NaN;
    varargout{4} = NaN;
    return;
end


%% Delineation
warning off

i_nD = i_nD( ~isnan(i_nD(:)) );
i_nD = 1 + round ( i_nD*fsi );

t           =	0:1/fs:  (length(signal)-1)/fs;
t_i         =	0:1/fsi:((length(signal)*(fsi/fs)-1)/fsi);
signal_i	=   interp1( t , signal , t_i , 'spline' );


%%
if nargout>=1
    
    % nA
    mtx_nA = repmat( 0:round(wdw_nA*fsi) ,	length(i_nD) , 1 ) + i_nD;
    mtx_nA(mtx_nA<1)=1; mtx_nA(mtx_nA>length(signal_i))=length(signal_i);
    [~,i_nA] = max ( signal_i(mtx_nA),[],2 ); i_nA = i_nA + i_nD;
    i_nA(i_nA<1 | i_nA>length(signal_i)) = NaN;
    
    nA = NaN(length(i_nA),1);
    nA(~isnan(i_nA)) = t_i(i_nA(~isnan(i_nA)));
    
    varargout{1} = nA;
    
end


%%
if nargout>=2
    
    % nB
    mtx_nB = repmat( -round(wdw_nB*fsi):0 ,	length(i_nD) , 1 ) + i_nD;
    mtx_nB(mtx_nB<1)=1; mtx_nB(mtx_nB>length(signal_i))=length(signal_i);
    [~,i_nB] = min ( signal_i(mtx_nB),[],2 ); i_nB = i_nB +(i_nD-round(wdw_nB*fsi));
    i_nB(i_nB<1 | i_nB>length(signal_i)) = NaN;
    
    nB = NaN(length(i_nB),1);
    nB(~isnan(i_nB)) = t_i(i_nB(~isnan(i_nB)));
    
    varargout{2} = nB;
    
end


%%
if nargout>=3
    % nM
    
    nM = NaN(length(i_nD),1);
    for ii = 1:length(i_nD)
        
        if (isnan(i_nB(ii)) || isnan(i_nA(ii))), continue; end
        pulseAmplitude = (signal_i(i_nB(ii))+signal_i(i_nA(ii)))/2;
        mtx_nM = i_nB(ii):i_nA(ii); mtx_nM(mtx_nM<1)=1; mtx_nM(mtx_nM>length(signal_i))=length(signal_i);
        [~,i_nM] = max ( - abs( signal_i(mtx_nM) - pulseAmplitude' ) ,[],2 ); i_nM = i_nM + i_nB(ii) ;
        i_nM(i_nM<1 | i_nM>length(signal_i)) = NaN;
        if ~isnan(i_nM) || ~isempty(i_nM), nM(ii) = t_i(i_nM); end
        
    end
    
    %     %Estas dos i_xx_fix son solo para corregir cuando los i_nB o i_nA contienen algun NaN
    %     i_nB_fix = isnan(i_nB); i_nB(i_nB_fix) = i_nD(i_nB_fix) - round(wdw_nB/2*fsi); i_nB(i_nB<1)=1;
    %     i_nA_fix = isnan(i_nA); i_nA(i_nA_fix) = i_nD(i_nA_fix) + round(wdw_nA/2*fsi); i_nA(i_nA>length(signal_i))=length(signal_i);
    %
    %     pulseAmplitude = (signal_i(i_nB)+signal_i(i_nA))/2;
    %     mtx_nM = repmat( i_nB:i_nA ,	length(i_nD) , 1 );
    %     mtx_nM(mtx_nM<1)=1; mtx_nM(mtx_nM>length(signal_i))=length(signal_i);
    %     [~,i_nM] = max ( - abs( signal_i(mtx_nM) - pulseAmplitude' ) ,[],2 ); i_nM = i_nM + round(wdw_nB*fsi) ;
    %     i_nM(i_nM<1 | i_nM>length(signal_i)) = NaN;
    %
    %     nM = NaN(length(i_nM),1);
    %     nM(~isnan(i_nM)) = t_i(i_nM(~isnan(i_nM)));
    %
    %     nM ( i_nA_fix | i_nB_fix ) = NaN;
    %
    %     i_nB ( i_nB_fix ) = NaN;
    %     i_nA ( i_nA_fix ) = NaN;
    
    varargout{3} = nM;
    %     varargout{3} = NaN(length(i_nB),1);
end


%%
if nargout>=4
    
    % nD
    [bb,aa] = butter(5,25*2/fs,'low');
    diffSignal = nanfiltfilt(bb,aa,diffSignal); %#ok
    
    wdw_nD	=   round(0.030*fsi);
    mtx_nD = repmat( -wdw_nD:wdw_nD , length(i_nD) , 1 ) + i_nD;
    mtx_nD(mtx_nD<1)=1; mtx_nD(mtx_nD>length(signal_i))=length(signal_i);
    
    diffSignal_i	=   interp1( t , diffSignal , t_i , 'spline' );
    [~,i_nDi] = max ( diffSignal_i(mtx_nD),[],2 ); i_nDi = i_nDi + (i_nD-wdw_nD);
    i_nDi(i_nDi<1 | i_nDi>length(signal_i)) = NaN;
    
    nD = NaN(length(i_nDi),1);
    nD(~isnan(i_nDi)) = t_i(i_nDi(~isnan(i_nDi)));
    
    varargout{4} = nD;
    
end


%%
if plotflag
    
    i_nM = round(nM*fsi);
    
    figure
    plot(t_i,signal_i); hold on
    plot(t_i(i_nDi),signal_i(i_nDi), 'r*');
    plot(t_i(i_nA),signal_i(i_nA),'ro');
    plot(t_i(i_nB),signal_i(i_nB),'*');
    plot(t_i(i_nM),signal_i(i_nM),'x');
    plot(t_i(i_nD),signal_i(i_nD),'*');
    title('Pulse Delineation')
    
end

warning on

end



% %--------------------------------------------------------------------------
% function [iPk, iInflect] = findLocalMaxima(yTemp)
% % bookend Y by NaN and make index vector
% yTemp = [NaN; yTemp; NaN];
% iTemp = (1:length(yTemp)).';
%
% % keep only the first of any adjacent pairs of equal values (including NaN).
% yFinite = ~isnan(yTemp);
% iNeq = [1; 1 + find((yTemp(1:end-1) ~= yTemp(2:end)) & ...
%                     (yFinite(1:end-1) | yFinite(2:end)))];
% iTemp = iTemp(iNeq);
%
% % take the sign of the first sample derivative
% s = sign(diff(yTemp(iTemp)));
%
% % find local maxima
% iMax = 1 + find(diff(s)<0);
%
% % find all transitions from rising to falling or to NaN
% iAny = 1 + find(s(1:end-1)~=s(2:end));
%
% % index into the original index vector without the NaN bookend.
% iInflect = iTemp(iAny)-1;
% iPk = iTemp(iMax)-1;