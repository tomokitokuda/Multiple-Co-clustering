function options=setVBCCGauss(varargin)

%% Set default options
% Hyperparameters in Dirichlet Process
options.M = 5;   % Maximum number of views
options.AVW0 = 1;  % For views. alpha_1 in the paper
options.AV0 = 1; % For Gaussian features. alpha_2 in the paper
options.AZ0 = 1; % For objects. beta in the paper

% Gaussian
options.L0 = 1e-4; % lambda_0
options.M0 = 0; % mu_0 
options.G0 = 1; % gamma_0
options.S0 = 1e-4; % sigma_0^2

% Poisson
options.gamma1 = 1; % alpha_0 
options.gamma2 = 1; % beta_0

% Bernoulli 
options.beta1 = 1; % elements of rho. We uppose symmetric Dirichlet

% Number of iterations and runs
options.MaxRuns = 100; % Number of runs
options.MaxVBIter = 100; % Iteration for a single run.

% Option for allowing only single feature cluster in a view 
options.onefcluster = false; % true means only one f.cluster in a view.


%% Return the default options if input argments are empty
if(isempty(varargin))
    return;
end

%% Change the options according to the input arguments
for k=1:2:length(varargin)
    fname = varargin{k};
    if(isfield(options,fname))
        value=varargin{k+1};
        options.(fname)=value;
    end
end

end
