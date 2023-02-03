function [npProcessed, tidalVolume, npUpperEnvelope, npLowerEnvelope, artifacts] = nasalPressureProcessing(np,npFs)

[bb, aa] = butter(3, 0.1*2/npFs, 'high'); % AASM Recommended
npProcessed = filtfilt(bb, aa, np); clear bb aa
[bb, aa] = butter(3, 15*2/npFs, 'low'); % AASM Recommended
npProcessed = filtfilt(bb, aa, npProcessed); clear bb aa
npProcessed = normalize(npProcessed);
[tidalVolume, npUpperEnvelope, npLowerEnvelope, artifacts] = tdvol_est(npProcessed,npFs);

end