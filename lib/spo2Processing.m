function [spo2Processed] = spo2Processing(spo2,spo2Fs)

spo2Processed = round(spo2); % Remove Gibbs
spo2Processed(spo2Processed<60 | spo2Processed>100) = nan; % Remove impossible values
spo2Processed = medfilt1(spo2Processed,3*spo2Fs); % Remove spikes

end