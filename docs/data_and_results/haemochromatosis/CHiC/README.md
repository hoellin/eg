- [`preliminary_analysis_chic_v2.Rmd`](preliminary_analysis_chic_v2.Rmd) contains the Analysis of E-G Networks starting from genes involved in hemochromatosis / iron metabolism, based on CHiC data. It computes genes regulated by at least one enhancer regulating one of the 13 following initial genes:

  ```R
  c("HFE", "TFR2", "HFE2", "HAMP", "SLC40A1", "BMP6", "TMPRSS6", "TFRC", "SLC11A2", "CYBRD1", "NEO1", "CIAPIN1", "SLC39A14")
  ```

(Actually there are no data for the gene HAMP in our datasets).

Then it plots a few graphs, and saves a few output tables.

Note that it is designed in a way that makes the knitted html document easy-to-read, with only most important code and outputs.

Details on CHiC data pre-processing can be found in [`preprocess_data.md`](preprocess_data.md) or here in [preprocess_data.html](preprocess_data.html).

- [`GO_enrichment_FDR.R`](GO_enrichment_FDR.R) is a script to compute GO enrichment of any list of genes. See http://genoweb.toulouse.inra.fr/~thoellinger/fall_2021/custom_go.html for explanations + example ; or see the content of [`go_enrichment.md`](go_enrichment.md) in the current repository to see its applications to the list of genes obtained with [`preliminary_analysis_chic_v2.Rmd`](preliminary_analysis_chic_v2.Rmd).

- [`results/`](results/) contains files computed at different levels of our analysis, the details of which are given when they are introduced in the dedicated markdowns.

