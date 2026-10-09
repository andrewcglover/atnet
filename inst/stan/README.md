# Stan models

Four models, all fitted to the laboratory assays. Two feed the transmission
simulations; the other two are fitted only for comparison in the
supplementary information.

| File | Fits | Data | Used by |
|---|---|---|---|
| `eip_fit_hill.stan` | lengthening of the extrinsic incubation period, Hill recovery | sporozoite positivity, one row per experiment | the simulations, and the SI |
| `int_fit_nb_rate_hill_bidir_splitnH.stan` | transmission-reducing activity (TRA): the reduction in mean oocyst count | oocyst count of each mosquito | the simulations, and the SI |
| `si_comparison/eip_fit.stan` | the same as `eip_fit_hill.stan`, with an exponential recovery | as above | the SI only |
| `si_comparison/int_fit_linear_hill_bidir_splitnH.stan` | laboratory transmission-blocking activity (TBA): the reduction in the proportion of mosquitoes with any oocyst | as above, reduced to whether each mosquito is positive | the SI only |

The two blocking models share the Hill blocking function and its priors, and
differ only in what they are fitted to. Their names describe the likelihood:
in `nb_rate` the treated arm's mean count is the control mean thinned by the
blocking function, and in `linear` the treated arm's probability of infection
is the control probability multiplied by it.

The field transmission-blocking activity used in the simulations is not fitted.
It is obtained from the TRA posterior inside malariasimulation.

## Origin

Copied unchanged on 9 Oct 2026 from the repositories the current posteriors
were fitted in, so the first commit of each file here is that version:

| File | Copied from |
|---|---|
| `eip_fit_hill.stan` | `malariasimple_ATNs/dev/` |
| `int_fit_nb_rate_hill_bidir_splitnH.stan` | `malariasimple_ATNs/dev/` |
| `si_comparison/eip_fit.stan` | `malariasimple_ATNs_190526backup/dev/` |
| `si_comparison/int_fit_linear_hill_bidir_splitnH.stan` | `malariasimple_ATNs/dev/` |

The header comment in each file predates the move.

## Priors

Every prior value comes from `atn_priors()` in `R/priors.R`, passed in as data,
and a test fails if a model fixes one itself. Two were fixed inside the EIP
models when they were copied over, and have since been moved out at the same
values: the half-normal(0, 1) prior on the random-effect standard deviation in
both, and the truncation of the Hill exponent at exp(0.5), about 1.65, in
`eip_fit_hill.stan`.
