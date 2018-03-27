function Model=runVBCCGaussS(Xgauss, Xpois, Xber, options)

%% Introduction
% This is an algorithm to carry out multiple-clustering with co-clustering
% structure. This algorithm has the following outlook:
% - Optimally divide datasets into several views. In each view, both objects and 
%   and features are clustered. 
% - Three types of features can be simultaneously analyzed: 
%    Gaussian, Poisson and categorical
% - Number of views and clusters are automatically inferred.
%
% See more details in the following paper: 
%   Tokuda, T., Yoshimoto, J., Shimizu, Y., Okada, G., Takamura, M., 
%   Okamoto, Y., Yamawaki, S., & Doya, K. 
%   PLOS ONE, 12(10):e0186566 (2017)
%  
%% Input
% Xgauss: Gaussian type of data matrix N x Dg
% Xpois: Poisson type of data matrix N x Dp
% Xber: Bernouilli type (including categorical) of data matrix N x Db
%       Note: This data should inlclude only positive interests starting from
%             1. Do not use zero.
% Options: options for parameters (if it is not given, then default)
%          e.g. options = setVBCCGauss('MaxRuns',30, 'M', 5);
% 
% Note
%  1) The order of input should be: Xgauss, Xpois, and Xber (and options)
%  2) In case that some type of data are not used, the corresponding input should
%      be 'empty matrix'. 
%     e.g., if there is no data for Poisson type, then let Xpois=[]
%     Anyway, three datasets are always required as an input!
%  3) For Xber, the category should be expressed by 1:T where T is 
%        the maximum number of categories. (Do not start from zero)
%  4) For missing entries (nan), just leave as they are. 
%
%% Output
% Clustering results:
%   M2 : Estimate of the number of views
%        Note: The view is sorted in the descending order of the number 
%        of features (irrespective of type of features) in the view.
%   EZ : posterior probabilities of allocations of objects
%   EV : posterior probabilities of allocations of Gaussian features
%   EVp: posterior probabilities of allocations of Poisson features
%   EVb: posterior probabilities of allocations of Bernoulli features
%
%   MapZ : cluster solutions of objects based on EZ (allocation for a cluster with maximum probability)
%   MapV : cluster solutions of Gaussian features
%   MapVp: cluster solutions of Poisson features
%   MapVb: cluster solutions of Bernoulli features
%
% Hyperparameters of probabilistic distribution
%   (See the paper for details.)
%   Gaussian
%     Prior (see SSupportingupporiting Information S1 in the paper)
%      L0 : lambda_0
%      M0 : mu_0 
%      G0 : gamma_0
%      S0 : sigma_0^2
%     Posterior
%      LP: LP{v}{a, b} is the value of lambda in a normal-inverse gamma
%            for object cluster a, feature cluster, and view v.
%      MP: MP{v}{a, b} is the value of mu in a normal-inverse gamma
%            for object cluster a, feature cluster, and view v.
%      GP: GP{v}{a, b} is the value of gamma in a normal-inverse gamma
%            for object cluster a, feature cluster, and view v.
%      SP: SP{v}{a, b} is the value of sigma^2 in a normal-inverse gamma
%            for object cluster a, feature cluster, and view v.
%
%   Poisson
%     Prior 
%      gamma1: alpha_0 
%      gamma2: beta_0
%     Posterior
%      gamma1p: gammap1p{v}{a, b} is the value of alpha in a gamma distribution
%            for object cluster a, feature cluster, and view v.
%      gamma2p: gammap1p{v}{a, b} is the value of beta in a gamma distribution
%            for object cluster a, feature cluster, and view v.
%
%   Bernoulli (including categorical)
%     Prior
%      beta1: elements of rho. We suppose symmetric Dirichlet
%     Posterior
%      beta1p: beta1p{v}{t}(a, b) is the value of Dirichlet for 
%           category t in object cluster a, feature cluster b, and view v. 
%        
%
% Others
% LE : Maximum log likelihood
% emptyInd : Indicator of empty data set [Gauss Poisson Bernoulli]

%% Preprocessing
% If options is omitted, set them at default values
if(nargin<=3)
    options=setVBCCGauss;
end

% Get the size of data matrix. Empty data is replaced by nan vector

emptyInd = [isempty(Xgauss) isempty(Xpois) isempty(Xber)];

if ~emptyInd(1)
    [N,~]=size(Xgauss);
end
if ~emptyInd(2)
    [N,~]=size(Xpois);
end
if ~emptyInd(3)
    [N, ~] = size(Xber);
end
if emptyInd(1)
    Xgauss = nan*ones(N, 1);
end
if emptyInd(2)
    Xpois = nan*ones(N, 1);
end
if emptyInd(3)
    Xber = nan*ones(N, 1);
end

% Get maximum number of VB runs
MaxRuns=options.MaxRuns;

% Get maximum number of VB iteration loops for each run
MaxS=options.MaxVBIter;

% Get prior hyperparameters for Dirichlet Process 
M = options.M; % Number of views
AVW0 = options.AVW0; 
AZ0=options.AZ0; % Object clusters
AV0=options.AV0; % Feature clusters

% Get prior hyperparameters for Gaussian 
L0=options.L0; 
M0=options.M0;
G0=options.G0;
S0=options.S0;

% Poisson
gamma1 = options.gamma1;
gamma2 = options.gamma2;

% Bernouilli
beta1 = options.beta1;

% Option for restricting only one feature cluster in a view
onefcluster = options.onefcluster;

%% Multiple runs with different initial configurations
LE=-inf;
Model=[];

EZ=cell(1, M);
EVg=cell(1, M);
EVp=cell(1, M);
EVb=cell(1, M);

for k=1:MaxRuns
    disp(sprintf('Runs: %d out of %d', k, MaxRuns));
    if onefcluster % special case: One feature cluster in a view
        for m=1:M
            EZ{m}=randClassMatrix(size(Xgauss,1)); % Object cluster
            EVg{m}=ones(1, (size(Xgauss,2))); % Gauss Feature one cluster
            EVp{m}=ones(1, (size(Xpois,2))); % Poisson 
            EVb{m}=ones(1, (size(Xber,2))); % Bernouilli
        end
    else
        for m=1:M
            EZ{m}=randClassMatrix(size(Xgauss,1)); % Object cluster
            EVg{m}=randClassMatrix(size(Xgauss,2))'; % Gauss feature multiple clusters
            EVp{m}=randClassMatrix(size(Xpois,2))'; % Poisson
            EVb{m}=randClassMatrix(size(Xber,2))'; % Bernoulli
        end
    end

    % Run VB algorithm
    tmpModel=coreVBCCGaussS(Xgauss,Xpois,Xber, EZ,EVg,EVp,EVb, AZ0,AV0,AVW0, L0,M0,G0,S0,gamma1, gamma2, beta1, MaxS,M, onefcluster);
    
    % Update maximum LE
    if(tmpModel.LE>LE)
        LE=tmpModel.LE;
        Model=tmpModel;
    end
end

Model.emptyInd = emptyInd;

Model = rmfield(Model, 'MaxZ');
Model = rmfield(Model, 'MaxV');
Model = rmfield(Model, 'MaxVp');
Model = rmfield(Model, 'MaxVb');
Model = rmfield(Model, 'SZ');
Model = rmfield(Model, 'SV');


% Eliminate unnecessary clusters
sall = zeros(1, length(Model.MapZ));
for i=1:length(Model.MapZ)
    s1 = sum(Model.MapV{i})>0;%
    s2 = sum(Model.MapVp{i})>0;
    s3 = sum(Model.MapVb{i})>0;
    sall(i) = s1+s2+s3;
end
posiall = sall>0;
M2 = sum(posiall); % Estimated number of views

Model.M2 = M2;
fieldname = {'M2','EZ', 'EV', 'EVp', 'EVb', 'MapZ', 'MapV', 'MapVp', ...
    'MapVb', 'L0', 'M0', 'G0', 'S0', 'LP', 'MP', 'GP', 'SP', ...
    'gamma1', 'gamma2', 'gamma1p', 'gamma2p', 'beta1', ...
    'beta1p', 'LE', 'emptyInd'};

Model = orderfields(Model, fieldname);

% Extract only valid views
Model.EZ=Model.EZ(1:M2);
Model.EV=Model.EV(1:M2);
Model.EVp=Model.EVp(1:M2);
Model.EVb=Model.EVb(1:M2);
Model.MapZ=Model.MapZ(1:M2);
Model.MapV=Model.MapV(1:M2);
Model.MapVp=Model.MapVp(1:M2);
Model.MapVb=Model.MapVb(1:M2);
Model.LP=Model.LP(1:M2);
Model.MP=Model.MP(1:M2);
Model.GP=Model.GP(1:M2);
Model.SP=Model.SP(1:M2);
Model.gamma1p=Model.gamma1p(1:M2);
Model.gamma2p=Model.gamma2p(1:M2);
Model.beta1p=Model.beta1p(1:M2);




end

