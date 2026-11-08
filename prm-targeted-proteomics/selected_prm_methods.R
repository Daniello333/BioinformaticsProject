# Selected functions adapted from the five-stage MSc targeted-proteomics PRM pipeline.
# The original pipeline reads external FASTA, target-panel and DIA-NN datasets.
# This small self-contained excerpt demonstrates sequence digestion, peptide mass,
# precursor m/z calculation with fixed carbamidomethyl C, and iRT time transfer.
# NOT a replacement for complete peptide ranking / DIA filtering / scheduling.

masa_reszt_PRM <- c(
  A = 71.037113805,
  R = 156.101111050,
  N = 114.042927470,
  D = 115.026943065,
  C = 103.009184505,
  E = 129.042593135,
  Q = 128.058577540,
  G = 57.021463735,
  H = 137.058911875,
  I = 113.084064015,
  L = 113.084064015,
  K = 128.094963050,
  M = 131.040484645,
  F = 147.068413945,
  P = 97.052763875,
  S = 87.032028435,
  T = 101.047678505,
  W = 186.079312980,
  Y = 163.063328575,
  V = 99.068413945
)
masa_H2O_PRM <- 18.010564684
masa_protonu_PRM <- 1.007276466812
delta_CAM_C <- 57.021463735

# Original method from stage 02 of the submitted thesis source.
oblicz_mase_peptydu_PRM <- function(peptyd) {
  aa <- strsplit(as.character(peptyd), "", fixed = TRUE)[[1]]
  if (length(aa) == 0L || any(!aa %in% names(masa_reszt_PRM))) return(NA_real_)

  masa_bez_modyfikacji <- sum(masa_reszt_PRM[aa]) + masa_H2O_PRM
  liczba_C <- sum(aa == "C")
  masa_bez_modyfikacji + liczba_C * delta_CAM_C
}

# Original method from stage 02.
oblicz_mz_PRM <- function(peptyd, z) {
  z <- suppressWarnings(as.integer(z))
  if (length(z) != 1L || is.na(z) || z <= 0L) return(NA_real_)
  masa <- oblicz_mase_peptydu_PRM(peptyd)
  if (!is.finite(masa)) return(NA_real_)
  (masa + z * masa_protonu_PRM) / z
}

# Original digestion method from stage 01: trypsin cleaves after K/R except
# before P; Lys-C cleaves after K except before P (as coded in this thesis).
# A missed cleavage extends a peptide across adjacent digestion fragments.
trawienie_bialka <- function(sekwencja, enzym = "trypsyna", max_dpmc = 0L) {
  aa <- strsplit(sekwencja, "", fixed = TRUE)[[1]]
  n <- length(aa)
  if (n < 1) return(character(0))
  reszty_ciecia <- if (enzym == "trypsyna") c("K", "R") else "K"
  miejsca_ciecia <- integer(0)
  if (n >= 2) {
    miejsca_ciecia <- which(
      aa[-n] %in% reszty_ciecia & aa[-1] != "P"
    )
  }
  granice <- c(0L, miejsca_ciecia, n)
  liczba_fragmentow <- length(granice) - 1L
  peptydy <- character(0)
  for (dpmc in 0:max_dpmc) {
    if ((liczba_fragmentow - dpmc) <= 0) next
    for (i in seq_len(liczba_fragmentow - dpmc)) {
      poczatek <- granice[i] + 1L
      koniec <- granice[i + 1L + dpmc]
      peptydy <- c(peptydy, paste0(aa[poczatek:koniec], collapse = ""))
    }
  }
  unique(peptydy)
}

# Equivalent to stage 05: iRT = a * RT + b in each chromatographic gradient.
RT60_na_RT30 <- function(rt60, a60, b60, a30, b30) {
  if (length(a30) != 1 || !is.finite(a30) || a30 == 0) {
    stop("Invalid 30-minute calibration slope")
  }
  (a60 * rt60 + b60 - b30) / a30
}
