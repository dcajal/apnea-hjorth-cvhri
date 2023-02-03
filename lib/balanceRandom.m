function [ hjorthParameters,labels,subjects ] = balanceRandom( hjorthParameters,labels,subjects )

nNB = sum(labels==0);
nApnea = sum(labels==1);
nHypo = sum(labels==2);
[nMinorityClass,minorityClass] = min([nNB nApnea nHypo]);
minorityClass = minorityClass-1;

classesToReduce = 0:2;
classesToReduce(classesToReduce==minorityClass) = [];

for kk = 1:2
    % Select all the elements of the class to reduce
    hpAux = hjorthParameters(labels==classesToReduce(kk),:);
    subjectsAux = subjects(labels==classesToReduce(kk));
    
    % Remove cases from final matrix
    hjorthParameters(labels==classesToReduce(kk),:) = [];
    subjects(labels==classesToReduce(kk)) = [];
    labels(labels==classesToReduce(kk)) = [];
    
    % Select nMinorityClass random cases
    idx = randperm(size(hpAux,1));
    hjorthParameters = [hjorthParameters; hpAux(idx(1:nMinorityClass),:)]; %#ok<*AGROW> 
    labels = [labels; ones(nMinorityClass,1)*classesToReduce(kk)];
    subjects = [subjects; subjectsAux(idx(1:nMinorityClass))];
end

end

