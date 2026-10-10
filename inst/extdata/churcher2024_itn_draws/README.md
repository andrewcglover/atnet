# ITN parameter draws of Churcher et al. (2024)

Unmodified copies of three files from the `parameters/` folder of
<https://github.com/EllieSherrardSmith/Mosq-Net-Efficacy>, at commit
`057ac811879b2a58722c542d535cd829328eafb3`, kept here so that the projections
do not depend on that repository remaining available. They are redistributed
under its MIT licence, copied alongside as `LICENSE`.

| File | Net | MD5 |
|---|---|---|
| `pyrethroid_uncertainty_LancetGH2024` | pyrethroid-only | `912b291e27ad07c35b325912e5059971` |
| `pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024` | pyrethroid-PBO | `f8b28649337ed5783bfba7ce1e7af36e` |
| `pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024` | pyrethroid-pyrrole (chlorfenapyr) | `baf1ab0fbad49f12e53b67c5278d115c` |

Each is an R data frame saved with `saveRDS()` (read it with `readRDS()`),
holding 1000 posterior draws of the net parameters at each of 101 levels of
pyrethroid resistance, 0 to 1 in steps of 0.01. They are read and checked by
`read_itn_draws()`, which extracts `dn0`, `rn0` and `gamman` as the
repository's script "R code/Part 6 Functions to output parameters with
uncertainty.R" does.
