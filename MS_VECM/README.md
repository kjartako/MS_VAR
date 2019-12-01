# Stan files

There are two Stan files in this folder, the first difference being the covariance prior.
The second difference is the error correction term, which is either alpha*(y[t,1]-y[t,2]) or alpha*(y[t,1]-beta*y[t,2]).

See the readme file in the MS-VAR folder for more details.

# Example files

The folder includes two R-files that shows how the stan-files are used in practise:
* **MS_VECM_example_oilgas**, a two-dimensional,two-regime MS-VECM with 1 autoregressive lag, estimating the joint dynamics of UK natural gas and Brent oil prices. The regimes are applied to both the mean structure and the covariance structure.
* **MS_VECM_example_oilgas_beta**, same example, but this time with the LKJ covariance prior and the error correction term which includes the additional beta parameter. 


# Label-switching
The prior distribution for the mean/covariance structures is equal for each regime, making the labelling arbitrary.
When using more than one chain, the arbitrary labelling for each chain may not coincide.
In the example files, unique labelling is introduces by ordering on the first element of the covariance matrices.
This is achieved using the **labsw** file. 
