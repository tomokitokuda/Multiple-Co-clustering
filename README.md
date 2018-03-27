## Multiple co-clustering based on nonparametric mixture models


This algorithm carries out multiple-clustering with co-clustering structures as follows:

- Optimally divides datasets into several views. In each view, both objects and features are clustered.
- Simultaneously analyzes three types of features: Gaussian, Poisson and categorical features.
- Automatically infers the number of views and clusters.

See more details in <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0186566>


### Files
- *runVBCCGaussS.m*  
  Run the algorithm
- *summaryModel.m*  
  Summarize results, yielding view and cluster memberships. 
- *setVBCCGauss.m*  
      Set configurations of the algorithm 
- *coreVBCCGaussS.m*  
      Carry out the algorithm for a single run
- *randClassMatrix.m*  
      Randomize initial memberships
- *Example1.m*, *Example2.m*  
      Examples of the algorithm

### Main function: *runVBCCGaussS.m*
#### Input:
 - Xgauss: Gaussian type of data matrix N x Dg (sample size N, dimension of data Dg)
 - Xpois: Poisson type of data matrix N x Dp (sample size N, dimension of data is Dp)
 - Xber: Bernoulli type (including categorical) of data matrix N x Db (sample size N, dimension of data Db)
 - Options (can be omitted): options for parameters (if it is not given, then default)
              e.g. options = setVBCCGauss('MaxRuns', 30, 'M', 5);
              See more details in *setVBCCGauss.m*
#### Note
 - The order of input should be as follows: Xgauss, Xpois, Xber, and options (optional)
 - In case that some type of dataset are not available, the corresponding data matrix should be an 'empty matrix'.  
e.g., if there is no data for Poisson type, then let Xpois=[].  
The point is that *three* datasets are always required for input, whether they are avaialbe or not. 
 - For Xber, the category must be numerically represented by positive integers from 1 to T where T is the maximum number of categories (Do not start from zero).
 - For missing entries (nan), just leave as they are. The algorithm deals with missing entries in a Bayesian manner. 

#### Output:
- Relevant parameter estimates

See more details in *runVBCCGaussS.m*

### Summary function: *summaryModel.m*
#### Input:
- Object yielded by *runVBCCGaussS.m*
#### Output:
- View memberships
- Feature cluster memberships
- Object cluster memberships

See more details in *summaryModel.m*

### Examples
See *Example1.m* and *Example2.m*

### References
Tokuda, T., Yoshimoto, J., Shimizu, Y., Okada, G., Takamura, M., Okamoto, Y., Yamawaki, S., & Doya, K.  
  "Multiple co-clustering based on nonparametric mixture models with heterogeneous marginal distributions"  
  PLOS ONE, 12(10): e0186566 (2017)  
 
### Authors
**Tomoki Tokuda** <tomoki.tokuda@oist.jp> (1)  
**Junichiro Yoshimoto** <juniti-y@is.naist.jp> (2, 1)  
1. Okinawa Institute of Science and Techonology Graduate Univeristy, JAPAN  
2. Nara Institute of Science and Technology, JAPAN 

### Acknowledgments
This algorithm was developed with the support of the Strategic Research Program for Brain Sciences from Japan Agency for Medical Research and development, AMED. 
