function [ hjorthParameters,labels,subjects ] = balanceRandomBinary( hjorthParameters,labels,subjects )

nNB = sum(labels==0);
nApnea = sum(labels==1);
nMinorityClass = min([nNB nApnea]);
[~,majorityClass] = max([nNB nApnea]);
majorityClass = majorityClass-1;

% Select all the elements of the class to reduce
hpAux = hjorthParameters(labels==majorityClass,:);
subjectsAux = subjects(labels==majorityClass);

% Remove cases from final matrix
hjorthParameters(labels==majorityClass,:) = [];
subjects(labels==majorityClass) = [];
labels(labels==majorityClass) = [];

% Select nMinorityClass random cases from majority class cases
idx = randperm(size(hpAux,1));
hjorthParameters = [hjorthParameters; hpAux(idx(1:nMinorityClass),:)];
labels = [labels; ones(nMinorityClass,1)*majorityClass];
subjects = [subjects; subjectsAux(idx(1:nMinorityClass))];


end

