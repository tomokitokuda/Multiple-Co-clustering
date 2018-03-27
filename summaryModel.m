function view = summaryModel(model)
% This function summarizes results of clustering by the multiple co-clustering method.
% Input: esimtaed model yielded by runVBCCGaussS.m
%
% Output: view and clustering results 
% view{v} contains results of feature clustering and object clustering in
%         view v. The following structure is provided:
%      view{v}.objects: object cluster memberships          
%      view{v}.feature.Gauss: feature cluster memberships for Gaussian type
%      view{v}.feature.Poisson: feature cluster memberships for Poisson type
%      view{v}.feature.Bernoulli: feature cluster memberships for categorical type
% 
%  Note 1: Cluster memberships are given by positive integers.  
%       2: If a feature membership is zero in view{v}, the feature does not 
%          belong to view v.


% Number of views
numview = model.M2;

sprintf('The number of views : %d', numview)

view = cell(1, numview);
for i=1:numview
    view{i}.objects = model.MapZ{i};
    
    % Gauss
    if model.emptyInd(1)
       view{i}.features.Gauss = nan; % give nan if data matrix is empty
    else
      view{i}.features.Gauss = model.MapV{i};
    end
   
    % Poisson
    if model.emptyInd(2)
        view{i}.features.Poisson = nan;
    else
        view{i}.features.Poisson = model.MapVp{i};
    end
    
    % Bernoulli
    if model.emptyInd(3)
        view{i}.features.Bernoulli = nan;
    else
        view{i}.features.Bernoulli = model.MapVb{i};
    end
end


end