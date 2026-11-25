# =====================================================================
# 07g_build_generic_site.R
#   Builds the synthetic "site file" for the indicative 3x3 figure: nine
#   hypothetical settings, NOT Mali, each one cell of
#       transmission (EIR 5 / 50 / 500)  x  pyrethroid resistance (50 / 80 / 95%).
#
#   Each cell is a "region" in the site-file sense, so 07g_run_generic.R drives
#   them exactly as 07 drives Mali's nine admin-1 regions: one shared history
#   per cell, then per-realisation segments, then the arms.
#
#   The object is built by CLONING the Mali site file's structure and
#   overwriting its contents, so site::subset_site() and site::site_parameters()
#   keep working on a shape they already accept. Nothing of Mali's data survives
#   except the demographic life table and the column layout.
#
#   Timeline (365-day years from 1 Jan 2000, as the model counts them):
#       2000-2002  burn-in, no nets, case management only
#       2003/06/09 historical Pyr-CFP mass campaigns (triennial)
#       2009       branch year: realisations get their own ITN draw from here
#       2012-2018  the arms, campaigns 2012 / 2015 / 2018
#       2015-2017  the reporting window, = campaign_cycle(cycle = 2) as in Mali
#
#   Seasonality is MINT's "seasonal" Fourier profile for every cell (single
#   peak, mid-May), so only EIR and resistance vary across the grid.
#
#   Run from the atnet package root:  Rscript analysis/07g_build_generic_site.R
#   Writes data_private/site_files/GEN.rds and prints a validation summary.
# =====================================================================

pkgload::load_all(normalizePath("."), quiet = TRUE)   # atnet: read_itn_draws()

TEMPLATE <- file.path("data_private", "site_files", "MLI.rds")
OUT      <- file.path("data_private", "site_files", "GEN.rds")

# --- the grid -------------------------------------------------------
EIRS <- c(5, 50, 500)
RES  <- c(0.50, 0.80, 0.95)
grid <- expand.grid(eir = EIRS, res = RES, KEEP.OUT.ATTRS = FALSE)
grid$name_1 <- sprintf("E%03dR%02d", grid$eir, round(grid$res * 100))
stopifnot(!anyDuplicated(grid$name_1))

# --- timeline -------------------------------------------------------
START_YEAR    <- 2000L
HIST_LAST     <- 2011L                  # future_yr0 = 2012, n_future_years = 7
CAMPAIGN_YRS  <- c(2003L, 2006L, 2009L) # historical Pyr-CFP, triennial
YEARS         <- START_YEAR:HIST_LAST

# MINT's campaign month is March (peak day 138 => May; 2 months before).
# cd_interval = 365/12, and 07's grid puts month m at (m-1)*cd_interval.
CAMPAIGN_DAY  <- round(2 * 365 / 12)    # 61, early March
CAMPAIGN_COV  <- 0.95
TX_COV        <- 0.46
PROP_ACT      <- 0.31
NET_TYPE      <- "pyrethroid_pyrrole"   # = Pyr-CFP in net_of_site_type

# --- MINT's "seasonal" Fourier profile ------------------------------
# From the MINT workflow's set_seasonality(); g0 = [1], g = [2:4], h = [5:7].
# Scale is irrelevant: carrying_capacity() is K0 * rainfall(t) / R_bar.
MINT_SEASONAL <- c(0.285505, -0.1328150, 0.0109352, 0.0139190,
                   0.3253520, -0.1046750, 0.0779865)

# --- MINT's species mix ---------------------------------------------
# Churcher et al. 2024 supplementary appendix 6, sec 1.5.1: the three modelled
# species at a fixed relative abundance, not user-changeable.
MINT_PROP <- c(gambiae = 0.50, arabiensis = 0.25, funestus = 0.25)

# =====================================================================
site <- readRDS(TEMPLATE)
cells <- grid$name_1
n     <- length(cells)

# --- sites ----------------------------------------------------------
site$country <- "Generic"
site$sites <- data.frame(country = "Generic", iso3c = "GEN", name_1 = cells,
                         stringsAsFactors = FALSE)

# --- interventions --------------------------------------------------
# One region's rows as the column template, then rebuilt for every cell x year.
tmpl <- site$interventions[site$interventions$name_1 == site$interventions$name_1[1], ][1, ]
iv <- tmpl[rep(1L, n * length(YEARS)), ]
iv$name_1 <- rep(cells, each = length(YEARS))
iv$year   <- rep(YEARS, times = n)
iv$country <- "Generic"; iv$iso3c <- "GEN"
rownames(iv) <- NULL

# Case management on throughout; every other intervention off.
iv$tx_cov   <- TX_COV
iv$prop_act <- PROP_ACT
for (cc in c("irs_cov", "smc_cov", "rtss_cov", "r21_cov", "lsm_cov", "pmc_cov")) {
  iv[[cc]] <- 0
}
iv$irs_spray_rounds <- 0
iv$smc_n_rounds     <- 0

# Resistance: fixed for the whole run, by cell (MINT holds it constant too).
iv$pyrethroid_resistance <- grid$res[match(iv$name_1, grid$name_1)]

# Historical Pyr-CFP campaigns; nothing in any other year.
is_camp <- iv$year %in% CAMPAIGN_YRS
iv$itn_input_dist <- ifelse(is_camp, CAMPAIGN_COV, 0)
iv$itn_use        <- 0
iv$predicted_use  <- 0
iv$net_type       <- NET_TYPE
iv$itn_distribution_day <- CAMPAIGN_DAY
iv$bioassay_mortality   <- 1 - iv$pyrethroid_resistance

# Pyr-CFP efficacy for the shared-history campaigns (2003, 2006): the MEDIAN
# across ITN draws at each cell's resistance. Campaigns from branch_year take
# the realisation's own draw instead, so these only set the pre-branch nets.
itn <- read_itn_draws()
cfp <- itn[["cfp"]]
med_at <- function(res) {
  j <- which.min(abs(sort(unique(cfp$resistance)) - res))
  r <- sort(unique(cfp$resistance))[j]
  x <- cfp[abs(cfp$resistance - r) < 1e-9, ]
  c(dn0 = median(x$dn0), rn0 = median(x$rn0), gamman = median(x$gamman) * 365)
}
eff <- t(vapply(grid$res, med_at, numeric(3)))
rownames(eff) <- grid$name_1
iv$dn0    <- eff[iv$name_1, "dn0"]
iv$rn0    <- eff[iv$name_1, "rn0"]
iv$gamman <- eff[iv$name_1, "gamman"]
iv$rnm    <- 0.24   # Mali uses this constant for every distribution; it is below the
                    # minimum drawn Pyr-CFP rn0 (0.4164) at all three resistances, so it
                    # stays < rn when build_past_schedule() overrides rn0 from a draw.
iv$mean_retention <- tmpl$mean_retention # unused: 07 sets retention itself

site$interventions <- iv

# --- seasonality: MINT's seasonal profile for every cell -------------
site$seasonality$seasonality_parameters <- data.frame(
  country = "Generic", iso3c = "GEN", name_1 = cells,
  g0 = MINT_SEASONAL[1],
  g1 = MINT_SEASONAL[2], g2 = MINT_SEASONAL[3], g3 = MINT_SEASONAL[4],
  h1 = MINT_SEASONAL[5], h2 = MINT_SEASONAL[6], h3 = MINT_SEASONAL[7],
  stringsAsFactors = FALSE
)

# --- EIR ------------------------------------------------------------
site$eir <- data.frame(country = "Generic", iso3c = "GEN", name_1 = cells,
                       sp = "pf", eir = grid$eir, stringsAsFactors = FALSE)

# --- vectors --------------------------------------------------------
vs <- site$vectors$vector_species
spp <- names(MINT_PROP)
one <- do.call(rbind, lapply(spp, function(s) vs[vs$species == s, ][1, ]))
vsn <- one[rep(seq_len(nrow(one)), times = n), ]
vsn$name_1  <- rep(cells, each = length(spp))
vsn$country <- "Generic"; vsn$iso3c <- "GEN"
vsn$prop    <- rep(MINT_PROP[spp], times = n)
rownames(vsn) <- NULL
site$vectors$vector_species <- vsn

# Resistance must span the PROJECTION years too: build_future_schedule() matches
# grid_year (through (YEARS) + n_future_years) against this table, and an
# unmatched year gives NA, which makes net_draw()'s which.min() return integer(0).
# Mali's own table runs 2000-2050, so do the same.
RES_YEARS <- START_YEAR:2050L
pr <- expand.grid(name_1 = cells, year = RES_YEARS, stringsAsFactors = FALSE)
pr$iso3c <- "GEN"; pr$country <- "Generic"
pr$pyrethroid_resistance <- grid$res[match(pr$name_1, grid$name_1)]
site$vectors$pyrethroid_resistance <-
  pr[, c("iso3c", "country", "name_1", "year", "pyrethroid_resistance")]

# --- population: equal across cells ---------------------------------
pt <- site$population$population_total
pt1 <- pt[pt$name_1 == pt$name_1[1], ][1, ]
ptn <- pt1[rep(1L, n * length(YEARS)), ]
ptn$name_1 <- rep(cells, each = length(YEARS))
ptn$year   <- rep(YEARS, times = n)
ptn$country <- "Generic"; ptn$iso3c <- "GEN"
if ("pop" %in% names(ptn)) ptn$pop <- 1e6
if ("population" %in% names(ptn)) ptn$population <- 1e6
rownames(ptn) <- NULL
site$population$population_total <- ptn

# demography: the template's life table, restricted to the run's years.
site$demography <- site$demography[site$demography$year %in% YEARS, ]
site$demography$iso3c <- "GEN"

# Anything else keyed by region would be stale; drop what 07 never reads.
for (nm in c("boundary", "shape", "cases_deaths", "prevalence", "accessibility",
             "blood_disorders", "bias_correction")) {
  site[[nm]] <- NULL
}
site$metadata <- list(built_by = "analysis/07g_build_generic_site.R",
                      built_on = Sys.time(),
                      template_md5 = unname(tools::md5sum(TEMPLATE)),
                      note = "Synthetic indicative settings; NOT Mali.")

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(site, OUT)

# =====================================================================
# Validation summary
# =====================================================================
cat("\nWrote", OUT, " md5", unname(tools::md5sum(OUT)), "\n\n")
cat("Cells (", n, "):\n", sep = "")
print(data.frame(cell = grid$name_1, eir = grid$eir, resistance = grid$res,
                 dn0 = round(eff[, "dn0"], 3), rn0 = round(eff[, "rn0"], 3),
                 gamman_d = round(eff[, "gamman"], 0)), row.names = FALSE)

cat("\nYears", min(YEARS), "-", max(YEARS), "(history); campaigns",
    toString(CAMPAIGN_YRS), "on day", CAMPAIGN_DAY, "\n")

pk <- malariasimulation::peak_season_offset(list(
  model_seasonality = TRUE, g0 = MINT_SEASONAL[1], g = MINT_SEASONAL[2:4],
  h = MINT_SEASONAL[5:7], rainfall_floor = 0.001))
month_starts <- cumsum(c(1L, 31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L))
cat("Seasonality: MINT 'seasonal', peak day", pk,
    paste0("(", month.name[findInterval(pk, month_starts)], ")"),
    "-> campaigns in", month.name[(findInterval(pk, month_starts) - 1 - 2) %% 12 + 1], "\n")
cat("Species:", paste(sprintf("%s %.0f%%", names(MINT_PROP), 100 * MINT_PROP),
                      collapse = ", "), "\n")
cat("tx_cov", TX_COV, "| prop_act", PROP_ACT, "| SMC, IRS, RTS,S, R21, PMC, LSM all 0\n")

stopifnot(
  nrow(site$interventions) == n * length(YEARS),
  all(site$interventions$smc_cov == 0), all(site$interventions$irs_cov == 0),
  all(site$interventions$rnm < site$interventions$rn0),
  nrow(site$eir) == n,
  nrow(site$seasonality$seasonality_parameters) == n,
  nrow(site$vectors$vector_species) == n * length(spp),
  all(abs(tapply(vsn$prop, vsn$name_1, sum) - 1) < 1e-9)
)
cat("\nAll structural checks passed.\n")
