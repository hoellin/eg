- `preliminary_analysis_v9.Rmd` contains the Analysis of E-G Networks starting from genes involved in hemochromatosis / iron metabolism, based on ABC data obtained from Nasser 2021 paper, in liver and intestine. It computes genes regulated by at least one enhancer regulating one of the 12 following initial genes:

  ```R
  c("HFE", "TFR2", "HFE2", "HAMP", "SLC40A1", "BMP6", "TMPRSS6", "TFRC", "SLC11A2", "CYBRD1", "NEO1", "CIAPIN1")
  ```

  Then it plots a few graphs, and saves a few output tables.

  Note that it is designed in a way that makes the knitted html document easy-to-read, with only most important code and outputs.

- `preliminary_analysis_v10.Rmd` is the very same R markdown, except it produces a more comprehensive html output when knitted, with almost all codes.

- `GO_enrichment_FDR.R` is a script to compute GO enrichment of any list of genes. See `custom_go.md` for explanations + example ; or see the content of `go_enrichment.md` in the current repository to see its applications to the list of genes obtained with `preliminary_analysis_v9.Rmd`

- `results` contains files computed at different levels of our analysis, the details of which are given when they are introduced in the dedicated markdowns.

