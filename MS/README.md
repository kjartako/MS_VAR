# Stan files

There are two Stan files in this folder, the difference being the covariance prior.

See the readme file in the MS-VAR folder for more details.

# Example files

The folder includes one R-file that shows how the stan-files are used in practise:
* **MS_example_sim_2reg**, a standard two-dimensional,two-regime MS model (with no autoregressive lags).

# Label-switching
The prior distribution for the mean/covariance structures is equal for each regime, making the labelling arbitrary.
When using more than one chain, the arbitrary labelling for each chain may not coincide.
In the example files, unique labelling is introduces by ordering on the second element of the mean vectors.
This is achieved using the **labsw_2reg_mu** file. 

**Note that Stan may return warning messages indicating chains have not mixed.
This is due to the label-switching.** 
After unique labelling is applied to the Stan output, the chains should mix well.
