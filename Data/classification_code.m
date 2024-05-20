%% Feature analysis for comparing control and treated h1RPE7 cells
% Updated 16th May 2024

% Code for classification of healthy and treated h1RPE7 cells using
% mitochondrial network features developed from AutoMitoNetwork. 

%This classification algorithm uses PCA and LDA for dimension reduction,
%followed by a supprot vector machine classifier for classification.

%10-fold cross validation is used. 

% The data was split into two groups: training data (80% of total data) to train the classifier and test data (20% of total data) to evaluate the classifier accuracy. 

% Users who use this code are required to change the excel file names, as
% per line 35 and 36. Users may also change the startfeature and
% stopfeature (lines 28 and 29), indicating the columns from which to read
% the data from. 

%%
close all
clear all
clc

%% for all and positive cells only - change excel file
% Files to evaluate
% File_Loc = 'C:\Users\z5160427\OneDrive - UNSW\Documents\h1RPE7 IOA 1uM treatment data';

startfeature = 2; %First column to read data from excel file
stopfeature = 42; %last column to read data from excel file

%Morphology = 1:18
%Intensity = 20:27
%Texture = 28:41

control = readtable('Control_combined_features_per_mito_network.xlsx'); %Change according to file name
treated = readtable('IOA_combined_features_per_mito_network.xlsx'); %Change according to file name

% Define group labels
myGroups = ["Control";"Treated"];

% First channel to evaluate (generally 1) 
myFirst = 1;

%% Data Evaluation
% Load data
myStatsArray1 = table2array(control(:,startfeature:stopfeature));
myStatsArray2 = table2array(treated(:,startfeature:stopfeature));

%Delete NaN values
myStatsArray1 = myStatsArray1(all(~isnan(myStatsArray1),2),:);
myStatsArray2 = myStatsArray2(all(~isnan(myStatsArray2),2),:);

%Delete infinity values
myStatsArray1 = myStatsArray1(all(~isinf(myStatsArray1),2),:);
myStatsArray2 = myStatsArray2(all(~isinf(myStatsArray2),2),:);


%% combine to make one big matrix for control and treated

evalMatrix = [myStatsArray1; myStatsArray2];


%% Split data for 80% training and 20% for testing
% Mixing data 
myRSelect1 = randperm(size(myStatsArray1,1)).'; num1_80 = ceil(.8*size(myStatsArray1,1)); num1_20 = size(myStatsArray1,1)-num1_80-1;
myRSelect2 = randperm(size(myStatsArray2,1)).'; num2_80 = ceil(.8*size(myStatsArray2,1)); num2_20 = size(myStatsArray2,1)-num2_80-1;

% Organising classes into correct order
myStats1r = myStatsArray1(myRSelect1,:);
myStats2r = myStatsArray2(myRSelect2,:);

% Splitting for test and training using mixing data
myStatsArray1 = myStats1r(1:num1_80,:); % 80%, group 1
myStats1Rest = myStats1r(num1_80+1:end,:); % 20%, group 1
myStatsArray2 = myStats2r(1:num2_80,:); % 80%, group 2
myStats2Rest = myStats2r(num2_80+1:end,:); % 20%, group 2

evalMatrix = [myStatsArray1; myStatsArray2];


% Giving classes same order
myClasses = [repmat(myGroups(1),size(myStatsArray1,1),1);repmat(myGroups(2),size(myStatsArray2,1),1)]; % 80%
myClassesC = categorical(cellstr(myClasses),cellstr(myGroups));

myClasses20 = [repmat(myGroups(1),size(myStats1Rest,1),1);repmat(myGroups(2),size(myStats2Rest,1),1)]; % 20%
myClassesC20 = categorical(cellstr(myClasses20),cellstr(myGroups));

myClassesN = grp2idx(myClassesC);
classIndices = find([myClassesN;0]-[0;myClassesN]~=0);
classIndices = classIndices(1:end-1);


myClassesTest = grp2idx(myClassesC20);
classIndicesTest = find([myClassesTest;0]-[0;myClassesTest]~=0);
classIndicesTest = classIndicesTest(1:end-1);

%% PCA 
myMean = mean(evalMatrix(:,:));    
myStd = std(evalMatrix(:,:));
meas_norm = zscore(evalMatrix(:,:)); 
[coeff,score,latent,tsquared,explained,mu] = pca(meas_norm);

cutOff = 0.95; % keep 95 % of information

if cutOff ~= 0
    myIndex = 0;
    mySum = 0;
    while mySum <= cutOff && myIndex+1 <= size(explained,1)
        myIndex = myIndex+1;
        mySum = mySum + explained(myIndex)/100;
    end
end
if myIndex <3
    myIndex = 3;
end   


score = score(:,1:myIndex);

%Plot PCA
figure(); hold on;
plot3(score(1:size(myStatsArray1,1),1),score(1:size(myStatsArray1,1),2),score(1:size(myStatsArray1,1),3),'.r');
plot3(score(size(myStatsArray1,1)+1:end,1),score(size(myStatsArray1,1)+1:end,2),score(size(myStatsArray1,1)+1:end,3),'.b');
xlabel("Cannonical variable 1"); ylabel("Cannonical variable 2"); zlabel("Cannonical variable 3");
legend(myGroups(1),myGroups(2));
title("PCA")

%% LDA
MdlLinear = fitcdiscr(score,myClassesC,'DiscrimType','pseudolinear');
% Project into subspace
mySigma = MdlLinear.Sigma;
if size(mySigma,1) ~= size(mySigma,2)
    mySigma = diag(mySigma);
end

[W, LAMBDA] = eig(MdlLinear.BetweenSigma, MdlLinear.Sigma); 
lambda = diag(LAMBDA);
[lambda, SortOrder] = sort(lambda, 'descend');
W = W(:, SortOrder);    
evalMatrix_LDA = score*W;


%% SVM classifier
MdlLinear_SVM = fitcsvm(evalMatrix_LDA,myClassesC,'Crossval','on','KFold',10);
 
% Validation prediction
[validationPredictions, validationScores] = kfoldPredict(MdlLinear_SVM);

rocObj = rocmetrics(myClassesC,validationScores,MdlLinear_SVM.ClassNames);

figure; plot(rocObj,ClassNames=MdlLinear_SVM.ClassNames(2))
xlabel('False positive rate') 
ylabel('True positive rate')
title("ROC for training data, AUC: " + rocObj.AUC(2))


%% Visualize
% LDA
figure(); 
hold on;     
plot3(evalMatrix_LDA(classIndices(1):classIndices(2)-1,1), evalMatrix_LDA(classIndices(1):classIndices(2)-1,2), evalMatrix_LDA(classIndices(1):classIndices(2)-1,3),'.','color',[0.5 0.7 1],'MarkerSize', 12);
plot3(evalMatrix_LDA(classIndices(2):end,1), evalMatrix_LDA(classIndices(2):end,2), evalMatrix_LDA(classIndices(2):end,3),'.','color',[0.6350 0.0780 0.1840],'MarkerSize', 12);
legend({myGroups(1), myGroups(2)});
xlabel('Linear discriminant 1')
ylabel('Linear discriminant 2')
zlabel('Linear discriminant 3')
title("Training");
set( legend , 'fontsize' , 16 , 'location' , 'east', 'box' , 'off','edgecolor','g' )
set(gca,'Fontsize',16)


%% Test model on 20% of unseen data
myStatsRest_z = [myStats1Rest;myStats2Rest];
myStatsRest_z = myStatsRest_z(:,:);

% Converting data into same space
myStats1Rest_z = (myStatsRest_z-repmat(myMean,size(myStatsRest_z,1),1))./repmat(myStd,size(myStatsRest_z,1),1);
myStatsRest_pca = myStats1Rest_z*coeff;
myStatsRest_pca = myStatsRest_pca(:,1:myIndex);
myStatsRest_LDA = myStatsRest_pca*W;


%LDA plot for test data
figure(); 
hold on;     

plot3(myStatsRest_LDA(classIndicesTest(1):classIndicesTest(2)-1,1), myStatsRest_LDA(classIndicesTest(1):classIndicesTest(2)-1,2), myStatsRest_LDA(classIndicesTest(1):classIndicesTest(2)-1,3),'.','color',[0.5 0.7 1],'MarkerSize', 12);
plot3(myStatsRest_LDA(classIndicesTest(2):end,1), myStatsRest_LDA(classIndicesTest(2):end,2), myStatsRest_LDA(classIndicesTest(2):end,3),'.','color',[0.6350 0.0780 0.1840],'MarkerSize', 12);
legend({myGroups(1), myGroups(2)});
xlabel('Linear discriminant 1')
ylabel('Linear discriminant 2')
zlabel('Linear discriminant 3')
title("Test");
set( legend , 'fontsize' , 12 , 'location' , 'east', 'box' , 'off','edgecolor','g' )
set(gca,'Fontsize',12)


% Predicting class
%Note can select which fold number for MdlLinear_SVM.Trained. Here have
%selected the first. 
[ypred_test,yci_test] = predict(MdlLinear_SVM.Trained{1}, myStatsRest_LDA);


% Receiver operating characteristics curve
[Xtest,Ytest,Tntest,auctest] = perfcurve(myClassesC20,yci_test(:,1),myGroups(1));


figure()
plot(Xtest,Ytest)
xlabel('False positive rate') 
ylabel('True positive rate')
title("ROC for Test data, AUC: " + auctest)


