# Second provisioning step for the DIDE cluster, run after pkgdepends.txt:
#   hipercow::hipercow_provision(method = "script", script = "provision_site.R")
#
# Installs the site package WITHOUT its dependencies. Its Remotes field asks for
# upstream mrc-ide/malariasimulation, which would conflict with (or replace) the
# ATN fork; its imports are already installed by pkgdepends.txt. repos = NULL
# makes install.packages() install this one package and nothing else.
site_sha <- "b1048b4884b0f19cd9d50c496d283f675876ec6c"   # site 1.0.13
install.packages(sprintf("https://github.com/mrc-ide/site/archive/%s.tar.gz", site_sha),
                 repos = NULL, type = "source")

# Check: site loads, and the malariasimulation it sees is still the fork.
stopifnot(packageVersion("site") == "1.0.13",
          "omega_atn" %in% names(malariasimulation::get_parameters()))
