# Base-R checks for the published illustrative PRM methods.
source("selected_prm_methods.R")
stopifnot(is.finite(oblicz_mz_PRM("ACDK", 2L)))
stopifnot(oblicz_mase_peptydu_PRM("C") > oblicz_mase_peptydu_PRM("A"))
stopifnot(is.na(oblicz_mz_PRM("ACDK", 0)))
stopifnot(identical(trawienie_bialka("AKRPR", "trypsyna", 0L), c("AK", "RPR")))
stopifnot(isTRUE(all.equal(RT60_na_RT30(10, a60=2, b60=1, a30=4, b30=1), 5)))
cat("Selected PRM helper checks passed\n")
