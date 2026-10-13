# Second provisioning step for the DIDE cluster, run after pkgdepends.txt:
#   hipercow::hipercow_provision(method = "script", script = "provision_site.R")
#
# Installs the site package WITHOUT its dependencies. Its Remotes field asks for
# upstream mrc-ide/malariasimulation, which would conflict with (or replace) the
# ATN fork; its imports are already installed by pkgdepends.txt. repos = NULL
# makes install.packages() install this one package and nothing else.
#
# Pinned to site 1.0.10, the frozen snapshot tagged malariaverse_01_2025: the
# version every earlier Mali run used. Later versions (from June 2026) refuse the
# Mali site file (malariaverse_no_urban_rural_06_2025) as an older site-file
# structure.
site_sha <- "b32be32aceb86ab53734a776889ec3be33dbddf6"   # site 1.0.10
install.packages(sprintf("https://github.com/mrc-ide/site/archive/%s.tar.gz", site_sha),
                 repos = NULL, type = "source")

# Check: site loads, and the malariasimulation it sees is still the fork.
stopifnot(packageVersion("site") == "1.0.10",
          "omega_atn" %in% names(malariasimulation::get_parameters()))
