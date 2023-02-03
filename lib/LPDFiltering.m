function [ signalLPD , filterLPD ] = LPDFiltering( signal , fs , Setup )
% Low Pass Derivative Filter
%
% --------------------------------------------------------------------------------
%   Sintax: [ signalLPD , filterLPD ] = LPDFilter( signal , fs , Setup )
%   In:   signal	=	input signal (PPG or BP)
%         fs        =	sampling rate [Hz]
%         filterLPD =   LPD filter (default: empty)
%         orderLPD	=	filter order [Default: 3*fs]
%         fpLPD	=	last frequency of the pass-band for the LPD filter (Hz) [Default: 7.8]
%         fcLPD	=	cut-off frequency for the low-pass-differentiator filter (Hz) [Default: 8]
%
%   Out:  signalLPD = 	signal LPD-filtered
%         filterLPD =	LPD filter impulse response(fir1s)
%
% Created by Jesús Lázaro  <jlazarop@unizar.es> in 2014
% Fixed   by Pablo Armañac <parmanac@unizar.es> in 2019
%

%% Check Inputs
if nargin <= 2,   Setup = struct();                       end
if nargin <  2,   error('Not enough input arguments.');   end

% Default Values
if ~isfield(Setup,'filterLPD'),	Setup.filterLPD = [];	end

if ~isfield(Setup,'orderLPD'),	Setup.orderLPD = 3*fs;	end
if ~isfield(Setup,'fpLPD'),     Setup.fpLPD = 7.8;      end
if ~isfield(Setup,'fcLPD'),     Setup.fcLPD = 8;        end

if ~isfield(Setup,'plotflag'),	Setup.plotflag = false;	end

% Get, assign and store the variable names
data_names = fieldnames(Setup);
for ii = 1:length(data_names), eval([ data_names{ii} ' = Setup.' data_names{ii} ';']); end
clear Setup data_names ii


%% LPD filter

if isempty(filterLPD)
    
    if mod(orderLPD,2)~=0, orderLPD = orderLPD+1; end %Force even order
    
    filterSpecs     =   fdesign.differentiator('n,fp,fst', orderLPD, fpLPD*2/fs, fcLPD*2/fs);
    % filterLPD   	=   design(filterSpecs, 'equiripple');
    filterLPD       =   design(filterSpecs, 'firls');
end

bb          =   filterLPD.Numerator*fs/(2*pi);
delay       =   round((numel(bb)-1)/2);

signalLPD	=   filter(bb, 1, signal);
signalLPD	=   [signalLPD(1+delay:end); zeros(delay, 1)];


if plotflag
    
    figure
    
    subplot(211);
    impz(bb,1,[],fs);
    title('Impulse response of FIR Differentiator FILTER');
    
    subplot(212);
    [H,W] = freqz(bb,1,2^15); f = W.*(fs/(2*pi)); H(f>2*fcLPD)=[]; f(f>2*fcLPD)=[];
    plot(f,abs(H)/max(abs((H)))); 
    xlabel('f[Hz]');    title('Frequency Response');
    
end

end