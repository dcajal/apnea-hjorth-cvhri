function [ labels ] = computeSegmentLabels( segments, abnormalSegments, apneaSegments, hypoSegments )
labels = zeros(size(segments,1),1);
containAbnormal = zeros(size(segments,1),1);
containApnea = zeros(size(segments,1),1);
containHypo = zeros(size(segments,1),1);

for aa = 1:size(abnormalSegments,1)
    onsetSegment = find(segments(:,2)>=abnormalSegments(aa,1),1); % First segment including abnormal burst
    endsetSegment = find(segments(:,1)<=abnormalSegments(aa,2),1,'last'); % Last segment abnormal burst
    containAbnormal(onsetSegment:endsetSegment) = true;
end; clear aa
labels(logical(containAbnormal)) = 1; % 1 for abnormal (later 1 for apnea, 2 for hypo), 0 for normal.

for aa = 1:size(apneaSegments,1)
    onsetSegment = find(segments(:,2)>=apneaSegments(aa,1),1); % First segment including apnea burst
    endsetSegment = find(segments(:,1)<=apneaSegments(aa,2),1,'last'); % Last segment including apnea burst
    containApnea(onsetSegment:endsetSegment) = true;
end; clear aa
for hh = 1:size(hypoSegments,1)
    onsetSegment = find(segments(:,2)>=hypoSegments(hh,1),1); % First segment including apnea burst
    endsetSegment = find(segments(:,1)<=hypoSegments(hh,2),1,'last'); % Last segment including apnea burst
    containHypo(onsetSegment:endsetSegment) = true;
end; clear hh
labels(logical(containHypo)) = 2; % hypo
labels(logical(containApnea)) = 1; % apnea (if contain apnea and hypo -> apnea)
end

