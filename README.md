# Dataland 2024 forecast
This Repo contains code and results for forecasting the 2024 election in Dataland for scenario A to F.

## Model specifies

$$p_{ki} {\sim N(\pi_{ki}^{'}, \sigma_{ki}^{2})},$$

$$\pi_{ki}^' = {\frac{\pi_{ki}}{\sum_{k=1}^K \pi_{ki}}}$$,

$$p_{ki} {= log(\theta_{r[i],t[i],k}^') + \alpha_{1r[i]k}},$$
$$\theta_{r[i],t[i],k[i]}^' {= \frac{\theta_{r[i],t[i],k[i]}}{\sum_{k=1}^K \theta_{r[i],t[i],k[i]}}},$$
$$log(\theta_{r[i],t[i]}) {\sim N_K(log(\theta_{r[i],t[i]-1,\Sigma_{\theta})}.$$

## Notes

- The model estimates party support for the four main parties CC, DGM, PDAl and SSP (in provinces where applicable), undecided are excluded from the analysis but can easily be added as fivth party.
- Due to computational limitations, the model was reduced:
    - Three chains run for 500 iterations (of which 250 were warmup)
    - Only party-eletion specific bias was included as a covariate in the mean equation of thhe model. It can easily be extended to additionally account for mode, pollster, etc. analogues to the party-election specific bias $\alpha1_{r}$.
 
## Folder structure


