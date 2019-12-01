# Stan files

There are two Stan files in this repository, the only difference being the covariance prior.
If you are indifferent to this prior choice, the LKJ prior results in significantly lower runtime.
Only the emphasized inputs are subject to user specification, for which the rest of the "inputs" are derived.  

### Shared user inputs for both files
* T - The number of observations.
* dim - The dimension of observations.
* **y** - A T x dim observation matrix.
* **ARdim** - The number of autoregressive terms in the mean structure.
* **nreg** - The number of regimes.
* **mean_reg** - A boolean value indicating if there is a unique mean structure for each regime.
* **sigma_reg** - A boolean value indicating if there is a unique covariance structure for each regime.
* **Q_alpha** - Hyperparameter for the transition matrix, a vector of lenght nreg. The prior distribution for each row of the transition matrix is dirichlet(Q_alpha).
* **mu_mean** - Mean value of the normal prior for the constant means
* **mu_sd** - Standard deviation of the normal prior for the constant means
* **phi_mean** - Mean value of the normal prior for the autoregressive coefficient
* **phi_sd** - Standard deviation of the normal prior for the autoregressive coefficient

### User inputs spesific for the Wishart covariance prior (implemented using Bartlett decomposition)
* eye - Simply the dim x dim identity matrix.
* L - Cholesky factor of Wishart scale matrix.
* **nu** - Wishart degrees of freedom.

### User inputs spesific for the LKJ correlation prior, with inverse-gamma prior for the variances
* **eta** - LKJ hyperparameter
* **gamma_alpha** - Inverse gamma hyperparameter (shape)
* **gamma_beta** - Inverse gamma hyperparameter (scale)


# Example files

The repository includes two R-files that shows how the stan-files are used in practise:
* **MS_VAR_example_sim**, a two-regime MS-VAR with 1 autoregressive lag, fitted to a simulated two-dimensional data. The regimes are applied to both the mean structure and the covariance structure.
* **MS_VAR_example_USmacro**, a two-regime MS-VAR with 4 autoregressive lags, fitted to a three-dimensional quarterly US macro data set, consisting of inflation, unemployment and an interest rate. To avoid excessive parameterization, the regimes are only applied to the covariance structure. 

By default, both files use the Wishart prior, but the LKJ prior can be used by uncommenting the relevant code chunk.  
Number of chains, burnin and samples are specified by the user. 

# Label-switching
The prior distribution for the mean/covariance structures is equal for each regime, making the labelling arbitrary.
When using more than one chain, the arbitrary labelling for each chain may not coincide.
In the example files, unique labelling is introduces by ordering on the first element of the covariance matrices.
This is achieved using the **labsw_2reg** file, which is only valid for nreg=2. 
