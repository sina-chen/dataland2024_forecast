# Dataland 2024 forecast
This repository houses the code and results pertaining to the 2024 election forecast in Dataland, covering scenarios A through F. This forecast was conducted using R version 4.3.1 and rStan version 2.26.22.

## Model

$$p_{ki} {\sim N(\pi_{ki}^{'}, \sigma_{ki}^{2})},$$

$$\pi_{ki}^{'} = {\pi_{ki}/\sum_{k=1}^{K} \pi_{ki}},$$

$$log(\pi_{ki}) = {log(\theta_{r[i],t[i],k}^{'}) + \alpha_{1r[i]k}},$$

$$\theta_{r[i],t[i],k}^{'} = {\theta_{r[i],t[i],k}/\sum_{k=1}^{K} \theta_{r[i],t[i],k}},$$

$$log(\theta_{r[i],t[i]}) {\sim N_{K}(log(\theta_{r[i],t[i]-1}),\Sigma_{\theta})}.$$

Of particular interest here is $\theta_{r[i],t[i],k}^{'}$, representing the mean expected vote share for party $k$ on day $t$ during election $r," with $i$ denoting the specific poll. A softmax transformation is applied to ensure that the mean expected vote shares sum to 1 across all $K$ parties. The original, untransformed mean expected vote shares follow a multivariate normal distribution and follow a random walk. The initial values are  a K-length vector denoted as $\theta_{r[i], 0}$, which corresponds to the logarithm of the previous election results.

## Notes

- The model calculates party support for the four primary parties: CC, DGM, PDAl, and SSP (where applicable in provinces). Undecided voters are excluded from the analysis, although they can be readily incorporated as a fifth party.
- The national mean vote shares and win probabilities for scenarios with limited polling data, covering only a subset of provinces, are derived from that specific subset.
- Owing to computational constraints, the model was reduced in the following ways:
    - Three chains were executed for 500 iterations (comprising 250 warmup iterations). Full convergence was not achieved, and increasing the number of iterations may enhance convergence.
    - In the mean equation of the model, only party-election specific bias was integrated as a covariate. It can be readily expanded to incorporate covariates like the mode, pollster, and others, akin to the party-election specific bias denoted as $\alpha1_{r}.
 
## Folder structure

- The [code](https://github.com/sina-chen/dataland2024_forecast/edit/main/code/) folder contains all code
    - the [prepare_data_dataland.R](https://github.com/sina-chen/dataland2024_forecast/edit/main/code/prepare_data_dataland.R) file prepares the data input for the analysis.
    - the [fit_dataland1985_2024_scenario.R](https://github.com/sina-chen/dataland2024_forecast/edit/main/code/fit_dataland1985_2024_scenario.R) file fits the Stan model. 
    - the [res_dataland.R](https://github.com/sina-chen/dataland2024_forecast/edit/main/code/res_dataland.R) computes and saves the  mean expected vote shares and win probybilities for all province- and national-level forecasts. 
    - the [forecast_model_sandbox.stan](https://github.com/sina-chen/dataland2024_forecast/edit/main/code/forecast_model_sandbox.stan) contains the hierarchical Bayesian stan model. 
  


