# Dataland 2024 forecast
This Repo contains code and results for forecasting the 2024 election in Dataland for scenario A to F.

## Model

$$p_{ki} {\sim N(\pi_{ki}^{'}, \sigma_{ki}^{2})},$$

$$\pi_{ki}^{'} = {\pi_{ki}/\sum_{k=1}^{K} \pi_{ki}},$$

$$p_{ki} = {log(\theta_{r[i],t[i],k}^{'}) + \alpha_{1r[i]k}},$$

$$\theta_{r[i],t[i],k}^{'} = {\theta_{r[i],t[i],k}/\sum_{k=1}^{K} \theta_{r[i],t[i],k}},$$

$$log(\theta_{r[i],t[i]}) {\sim N_{K}(log(\theta_{r[i],t[i]-1}),\Sigma_{\theta})}.$$

Of iterest here is $\tehta_{r[i],t[i],k}^{'}$, which is the the mean expected vote share for party $k$, on day $t$, in election $r$, with $i$ identifying the poll. A softmax transfornation is applied, to ensure that the mean expected vote shares sum to 1 across all $K$ parties. The untransformed mean expected vote shares follow
a random walk, with the K length vector of starting values $\theta_{r[i], 0}$ corresponding to the log of the previous election results. 

## Notes

- The model estimates party support for the four main parties CC, DGM, PDAl and SSP (in provinces where applicable), undecided are excluded from the analysis but can easily be added as fivth party.
- Due to computational limitations, the model was reduced:
    - Three chains run for 500 iterations (of which 250 were warmup)
    - Only party-eletion specific bias was included as a covariate in the mean equation of thhe model. It can easily be extended to additionally account for mode, pollster, etc. analogues to the party-election specific bias $\alpha1_{r}$.
 
## Folder structure


