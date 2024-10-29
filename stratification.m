clear all;

load(strcat('results/classification/binaryClassification_noseverehypo.mat'), ...
    'corrects','cvhriTest','labelsTest','predictionsTest','subjectTest');
load(strcat('results/AHI_new.mat'));
addpath(genpath('models'));

% Exclude test errors
% 33 was wake most of the time
% 57 problems with nasal pressure signal. Not reliable annotations
exclude = [33 57];

for kk = 1:numel(exclude)
    labelsTest(subjectTest==exclude(kk)) = nan;
    predictionsTest(subjectTest==exclude(kk)) = nan;
    cvhriTest(exclude(kk)) = nan;
    corrects(exclude(kk)) = nan;
    ahiDataset(exclude(kk)) = nan;
end; clear kk exclude


%% Support Vector Machine
%  Stratification AHI >= 5

rng('default');

isAhiHigherThan5 = ahiDataset>=5;
cvhriAhiLowerThan5 = cvhriTest(~isAhiHigherThan5 & ~isnan(cvhriTest));
cvhriAhiHigherThan5 = cvhriTest(isAhiHigherThan5 & ~isnan(cvhriTest));
nObservations = numel(cvhriAhiHigherThan5)+numel(cvhriAhiLowerThan5);

cvhriAhiForHigherThan5 = [cvhriAhiHigherThan5; cvhriAhiLowerThan5;];
labelsForHigherThan5 = [true(numel(cvhriAhiHigherThan5),1); false(numel(cvhriAhiLowerThan5),1)];
costsForHigherThan5 = [0 nObservations/(2*numel(cvhriAhiLowerThan5)); nObservations/(2*numel(cvhriAhiHigherThan5)) 0];


% Cross varidation (train: 50%, test: 50%)
cv = cvpartition(size(cvhriAhiForHigherThan5,1),'HoldOut',0.5);
idx = cv.test;
% Separate to training and test data
dataTrain = cvhriAhiForHigherThan5(~idx,:);
dataTest  = cvhriAhiForHigherThan5(idx,:);
labelsTrain = labelsForHigherThan5(~idx,:);
labelsTest = labelsForHigherThan5(idx,:);

trainedClassifier = trainStratificationSVM(dataTrain, labelsTrain, costsForHigherThan5);
predictionsTest = trainedClassifier.predictFcn(dataTest);
c5 = confusionmat(labelsTest,predictionsTest) %#ok<*NASGU> 
isCorrect = ~xor(labelsTest,predictionsTest);

threshold5 = (min(dataTest(predictionsTest==1))+max(dataTest(predictionsTest==0)))/2

figure
plot(dataTest(isCorrect & predictionsTest==0),0,'bo'); hold on;
plot(dataTest(isCorrect & predictionsTest==1),0,'ro');
plot(dataTest(~isCorrect & predictionsTest==0),0,'bx');
plot(dataTest(~isCorrect & predictionsTest==1),0,'rx');
xlabel('CVHRI')
xline(threshold5,'--')

%%
% Stratification AHI >= 15

rng('default');

isAhiHigherThan15 = ahiDataset>=15;
cvhriAhiLowerThan15 = cvhriTest(~isAhiHigherThan15 & ~isnan(cvhriTest));
cvhriAhiHigherThan15 = cvhriTest(isAhiHigherThan15 & ~isnan(cvhriTest));
nObservations = numel(cvhriAhiHigherThan15)+numel(cvhriAhiLowerThan15);

cvhriAhiForHigherThan15 = [cvhriAhiHigherThan15; cvhriAhiLowerThan15;];
labelsForHigherThan15 = [true(numel(cvhriAhiHigherThan15),1); false(numel(cvhriAhiLowerThan15),1)];
costsForHigherThan15 = [0 nObservations/2/numel(cvhriAhiLowerThan15); nObservations/2/numel(cvhriAhiHigherThan15) 0];


% Cross varidation (train: 50%, test: 50%)
cv = cvpartition(size(cvhriAhiForHigherThan15,1),'HoldOut',0.5);
idx = cv.test;
% Separate to training and test data
dataTrain = cvhriAhiForHigherThan15(~idx,:);
dataTest  = cvhriAhiForHigherThan15(idx,:);
labelsTrain = labelsForHigherThan15(~idx,:);
labelsTest = labelsForHigherThan15(idx,:);

trainedClassifier = trainStratificationSVM(dataTrain, labelsTrain, costsForHigherThan15);
predictionsTest = trainedClassifier.predictFcn(dataTest);
c15 = confusionmat(labelsTest,predictionsTest)
isCorrect = ~xor(labelsTest,predictionsTest);

threshold15 = (min(dataTest(predictionsTest==1))+max(dataTest(predictionsTest==0)))/2

figure
plot(dataTest(isCorrect & predictionsTest==0),0,'bo'); hold on;
plot(dataTest(isCorrect & predictionsTest==1),0,'ro');
plot(dataTest(~isCorrect & predictionsTest==0),0,'bx');
plot(dataTest(~isCorrect & predictionsTest==1),0,'rx');
xlabel('CVHRI')
xline(threshold15,'--')


%% Linear regression
% https://es.mathworks.com/matlabcentral/answers/183690-what-is-the-difference-between-lar-and-the-bisquare-remain-robust-in-regression-curve-fitting-tool

% [xData, yData] = prepareCurveData( cvhriTest, ahiDataset);
% 
% % Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Robust = 'Bisquare';
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'cvhriTest vs. ahiDataset', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'cvhriTest', 'Interpreter', 'none' );
% ylabel( 'ahiDataset', 'Interpreter', 'none' );
% grid on
% 
% ahiFromCVHRI = fitresult.p1*xData+fitresult.p2;
% ahiFromCVHRIHigherThan5 = ahiFromCVHRI>=5;
% ahiHigherThan5 = yData>=5;
% 
% c5 = confusionmat(ahiHigherThan5,ahiFromCVHRIHigherThan5)
% % plotconfusion(ahiHigherThan5',ahiFromCVHRIHigherThan5');
% 
% 
% ahiFromCVHRIHigherThan15 = ahiFromCVHRI>=15;
% ahiHigherThan15 = yData>=15;
% 
% c15 = confusionmat(ahiHigherThan15,ahiFromCVHRIHigherThan15)
% % plotconfusion(ahiHigherThan15',ahiFromCVHRIHigherThan15');


