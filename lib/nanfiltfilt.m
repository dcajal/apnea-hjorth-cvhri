function SIGNAL_out = nanfiltfilt( bb, aa , SIGNAL_in , Setup )
% filtfilt supporting nan data. 
% 
% 
% Created by Pablo Armañac <parmanac@unizar.es> and Spyridon Kontaxis <sikontax@unizar.es> in 2018
%
% 
%% Parse Inputs

% Check Inputs
if nargin < 4,   Setup = struct();                       end
if nargin < 3,   error('Not enough input arguments.');   end

% Default Values
if ~isfield(Setup,'plotflag'),         Setup.plotflag = false;	end

% Get the variable names
data_names          = fieldnames(Setup);

% Assign and store the variables
for ii = 1:length(data_names)
    eval([ data_names{ii} ' = Setup.' data_names{ii} ';'])
end
clear Setup data_names ii

SIGNAL_in      =	SIGNAL_in (:);


%% Algorithm

% Find NaNs (>T)
idx_nan = isnan(SIGNAL_in);

if ~any(~idx_nan), 
    SIGNAL_out = SIGNAL_in;
    return
end
    
% Fill missing
SIGNAL_filled = fillmissing(SIGNAL_in,'linear');

% Filter
SIGNAL_filtered = filtfilt(bb,aa,SIGNAL_filled);

% Delete NaN segments
SIGNAL_out              = SIGNAL_filtered;
SIGNAL_out (idx_nan)	= NaN;


%% plot
if plotflag
    figure
    plot(SIGNAL_filled,'DisplayName','filled'); hold on
    plot(SIGNAL_in,'DisplayName','original');    
    plot(SIGNAL_filtered,'DisplayName','filled & filtered'); 
    plot(SIGNAL_out,'DisplayName','final');   
    legend
end


end