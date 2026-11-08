# Targeted proteomics PRM planning — selected code

**Status: representative research functions, not the full five-stage thesis pipeline.**

This directory contains selected methods from the submitted five-script R project:

- [`selected_prm_methods.R`](selected_prm_methods.R) — peptide mass with fixed cysteine carbamidomethylation, precursor charge-state m/z, trypsin/Lys-C digestion (as implemented in the thesis), and gradient-to-gradient iRT transfer.
- [`test_selected_prm_methods.R`](test_selected_prm_methods.R) — basic checks for the stand-alone excerpt.

```bash
cd prm-targeted-proteomics
Rscript test_selected_prm_methods.R
```

**Full experimental workflow is not yet hosted here.** The source archive contains five substantive R analysis stages (FASTA processing; DIA-NN evidence; precursor list assembly; retention crowding; iRT-based extensions), and result summaries. Its original data and DIA-NN reports are external and are required for a full rerun.

The thesis candidate lists covered **596 / 741** panel proteins at baseline and **609 / 741** after adding eligible 60-minute evidence. The reported 30-minute RT transfer MAE is approximately **0.333 minutes**. These are exported research results, not rederived by the stand-alone excerpt.

**Further reading:** [PRM research case study](https://github.com/Daniello333/Daniello333/blob/main/projects/prm-proteomics.md).

**Author:** Daniel Drabik · Jagiellonian University
