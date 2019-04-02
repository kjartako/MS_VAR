# MS_VAR

The repository contains code used in the following paper:

["MCMC for Markov-switching models — Gibbs sampling vs. marginalized likelihood"](https://www.tandfonline.com/doi/full/10.1080/03610918.2019.1565580),

by Kjartan Kloster Osmundsen, Tore Selland Kleppe & Atle Oglend.

The code fits a Markov-switching vector autoregressive model to data input. The user can specify the data, number of regimes and number of autoregressive terms. The regimes can be applied to the mean structure and/or the covariance structure.

See the wiki for detailed instructions. See also this [blog post](https://www.kjartako.no/post/stan-code-for-markov-switching-vector-autoregressive-models/).

The code assumes that the R-packages rstan and coda are installed (and rstudioapi if you are using Rstudio).
