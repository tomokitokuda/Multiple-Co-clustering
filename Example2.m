%% Example 2 (Gausian type + Categorical type)

%% Data generation
% Getting data matrix Xg and true cluster solutions
load('Xg.mat'); % Use the same Xg as in Example 1. 
load('objectcluster1.mat'); % Object clusters in Example 1
load('objectcluster2.mat'); % Object clusters in Example 1

% Generation of data matrix Xb
rng(1);
% View 1 (one feataure cluster and two object clusters)
kk1 = randsample(1:2, 150*5, true, [0.1 0.9]);
Xb1 = reshape(kk1, 150, 5);
kk2 = randsample(1:2, 150*5, true, [0.9 0.1]);
Xb2 = reshape(kk2, 150, 5);
objectcluster3 = [ones(150, 1); 2*ones(150, 1)];
per3 = randperm(300); 
Xb12 = [Xb1; Xb2];

Xb12 = Xb12(per3, :);
objectcluster3 = objectcluster3(per3);

% View 2 (one feature clsuter and three object clusters)
%        Note; the object cluster soluion is the same as that in Xg
kk3 = randsample(1:2, 100*3, true, [0.1 0.9]);
Xb3 = reshape(kk3, 100, 3);
kk4 = randsample(1:2, 100*3, true, [0.5 0.5]);
Xb4 = reshape(kk4, 100, 3);
kk5 = randsample(1:2, 100*3, true, [0.9 0.1]);
Xb5 = reshape(kk5, 100, 3);

Xb345 = [Xb3; Xb4; Xb5];
Xb = [Xb12 Xb345];

%% Implementation of the algorithm
options = setVBCCGauss('MaxRuns',30, 'M', 5);
Xp = []; % No data for Poisson
model = runVBCCGaussS(Xg, Xp, Xb, options); % main function
viewall = summaryModel(model); % Summarize results

%% Results: Examine the recovery of views and clusters.
% Estimated number of views
length(viewall)

% Feature clusters
viewall{1}.features.Gauss % Feature clusters
viewall{1}.features.Bernoulli

viewall{2}.features.Gauss % Feature clusters
viewall{2}.features.Bernoulli

viewall{3}.features.Gauss % Feature clusters
viewall{3}.features.Bernoulli

% Object clusters
[viewall{1}.objects objectcluster1] 
[viewall{2}.objects objectcluster2] 
[viewall{3}.objects objectcluster3] 
