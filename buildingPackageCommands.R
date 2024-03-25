sessionInfo()
setwd(here::here())
# detach("package:cutoff")
detach("package:cutoff2")
# unloadNamespace("cutoff")
unloadNamespace("cutoff2")
remove.packages("cutoff2", lib="/home/pleydell/R/x86_64-pc-linux-gnu-library/4.3")
#
setwd(here::here('cutoff2'))
# Remove Rd files (which cause bugs with deeted cutoff2.rdb file for some reason)
system(paste0("rm ", here::here("cutoff2/man"), "/*"))
#
roxygen2::roxygenise()
#
setwd(here::here())
devtools::build("cutoff2", vignettes=TRUE) # FALSE

# devtools::check("cutoff2")

devtools::install("cutoff2", build_vignettes = TRUE) # FALSE



library("cutoff2")
help.start()

setwd("cutoff2")
devtools::build_vignettes()
