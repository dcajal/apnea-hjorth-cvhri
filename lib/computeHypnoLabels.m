function [ hypnoLabels ] = computeHypnoLabels( segments, hypno, tHypno )

hypnoLabels = zeros(size(segments,1),1);

for kk = 1:size(segments,1)
    indexes = find(tHypno>=segments(kk,1),1):find(tHypno<=segments(kk,2),1,'last');
    hypnoLabels(kk) = mode(hypno(indexes));
end; clear kk

end

