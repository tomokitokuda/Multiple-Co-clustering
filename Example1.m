%% Example 1 (only Gaussian type)


%% Data generation
% view 1 (Two feature clusters and two object clusters)
rng(1);
Xg1 = [randn(150, 10)+3 randn(150, 5)-3];
Xg2 = [randn(150, 10)-3 randn(150, 5)+3];
Xg12 = [Xg1;Xg2];
objectcluster1 = [ones(150, 1); 2*ones(150, 1)]; % true object clusters

per2 = randperm(300); 
Xg12 = Xg12(per2, :); % random permutation of objects
objectcluster1  = objectcluster1(per2); 

% view 2 (One feature cluster and three object clusters)
Xg3 = randn(100, 10)-3;
Xg4 = randn(100, 10);
Xg5 = randn(100, 10)+3;
Xg345 = [Xg3; Xg4; Xg5];
objectcluster2 = [ones(100, 1); 2*ones(100, 1); 3*ones(100 , 1)];

% Combined view 1 and view 2. 
Xg12345 = [Xg12 Xg345];
Xg = zscore(Xg12345); % Standardize

%% Implementation of the algorithm
options = setVBCCGauss('MaxRuns',30, 'M', 5);
Xp = []; % No data for Poisson
Xb = []; % No data for Categorical
model = runVBCCGaussS(Xg, Xp, Xb, options); % main function
viewall = summaryModel(model); % Summarize results

%% Results: Examine the recovery of views and clusters.
% Estimated number of views
length(viewall)

% First view
viewall{1}.features.Gauss % Feature clusters
[viewall{1}.objects objectcluster1] % Object clusters. Compare with the true one

% Second view
viewall{2}.features.Gauss
[viewall{2}.objects objectcluster2] 