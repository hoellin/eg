# GO enrichment of all genes obtained

We use `GO_enrichment_FDR.R` which allows us to choose between multiple procedures for multiple testing correction. Here we are interested in controlling the FDR ; we can use either `pAdjustMethod = "BH"` (Benjamini-Hochberg, default parameter) or `pAdjustMethod = "fdr"`.

All results (list of genes and detailed GO enrichment) can be found here:

```
/work2/project/regenet/workspace/thoellinger/shared/2022/networks_hemochromatosis
```

## Dependencies

We use R 4.0. We need the `clusterProfiler` package. See `custom_go.md` for installation details.

## Main remarks

- From time to time, one may get the following error while using `GO_enrichment_FDR.R`:

  ```R
  Error in `[.data.frame`(d, , 2) : undefined columns selected
  ```

  This is expected, and should be considered as a warning only (I did not take the time to handle the case with an exception). This happens when no significant GO terms are found ; in such case, no graphical outputs are computed.

## Collecting the list of genes of interest

```bash
cd /work2/project/regenet/workspace/thoellinger/shared/2022/networks_hemochromatosis
```

We obtained this list using our R markdown doing network analysis of E-G pairs starting from genes directly involved in hemochromatosis or involved in the regulation of iron metabolism.

We saved this list of 457 (13 known + 444 inferred) genes as `new_genes_v7.list`. For each one of the 13 initial genes, we also saved separately the corresponding inferred genes:

> ```bash
> └── results
>  ...
>  ├── new_genes_v7.list # 457
>  └── separate
>      ├── BMP6.list # 17
>      ├── CIAPIN1.list # 56
>      ├── CYBRD1.list # 35
>      ├── HAMP.list # 35
>      ├── HFE2.list # 27
>      ├── HFE.list # 61
>      ├── NEO1.list # 39
>      ├── SLC11A2.list # 68
>      ├── SLC39A14.list # 35
>      ├── SLC40A1.list # 26
>      ├── TFR2.list # 34
>      ├── TFRC.list # 10
>      └── TMPRSS6.list # 14
> ```

## GO enrichment

```bash
module load system/R-4.0.4_gcc-9.3.0
```

### All 457 genes (original + inferred)

```bash
mkdir -p GO_FDR/all_genes/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/new_genes_v7.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/all_genes/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "457 provided genes; 336 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 3457 distinct GO terms"
> [1] "Of those 3457 GO terms, 29 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

![](images_go_enrichment_md/GO_FDR/all_genes/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/all_genes/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/all_genes/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/all_genes/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/all_genes/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/all_genes/symbol/default_universe.BP.0.05BH.pvalues.png)



### All 444 genes (inferred only ; w/o original genes)

```bash
mkdir -p GO_FDR/all_genes_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/new_genes_v7_without_original.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/all_genes_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "444 provided genes; 324 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 3297 distinct GO terms"
> [1] "Of those 3297 GO terms, 23 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

![](images_go_enrichment_md/GO_FDR/all_genes_wo_orig/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/all_genes_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/all_genes_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/all_genes_wo_orig/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/all_genes_wo_orig/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/all_genes_wo_orig/symbol/default_universe.BP.0.05BH.pvalues.png)

### All genes co-expressed in liver

```bash
mkdir -p GO_FDR/genes_inferred_from_that_expressed_in_liver/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/genes_inferred_from_that_expressed_in_liver.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/genes_inferred_from_that_expressed_in_liver/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "167 provided genes; 110 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 1595 distinct GO terms"
> [1] "Of those 1595 GO terms, 0 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

### All genes co-expressed in intestine

```bash
mkdir -p GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/genes_inferred_from_that_expressed_in_intestine.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "126 provided genes; 97 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 1522 distinct GO terms"
> [1] "Of those 1522 GO terms, 6 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

![](images_go_enrichment_md/GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/genes_inferred_from_that_expressed_in_intestine/symbol/default_universe.BP.0.05BH.pvalues.png)

### Each original gene separately

#### Inferred from BMP6

```bash
mkdir -p GO_FDR/separate/BMP6/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/BMP6.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/BMP6/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "17 provided genes; 12 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 420 distinct GO terms"
> [1] "Of those 420 GO terms, 151 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

Note that we tried at a more stringent fdr (`-f 0.05`), and it yielded no significant results ; so we kept the default 0.1 threshold for this one, and a few other ones later.

![](images_go_enrichment_md/GO_FDR/separate/BMP6/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from BMP6 w/o BMP6

```bash
mkdir -p GO_FDR/separate/BMP6_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/BMP6_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/BMP6_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "16 provided genes; 11 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 242 distinct GO terms"
> [1] "Of those 242 GO terms, 108 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

![](images_go_enrichment_md/GO_FDR/separate/BMP6_wo_orig/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6_wo_orig/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6_wo_orig/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/BMP6_wo_orig/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from CIAPIN1

```bash
mkdir -p GO_FDR/separate/CIAPIN1/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/CIAPIN1.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/CIAPIN1/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "56 provided genes; 48 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 1141 distinct GO terms"
> [1] "Of those 1141 GO terms, 23 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1/symbol/default_universe.BP.0.05BH.pvalues.png)

#### Inferred from CIAPIN1 w/o CIAPIN1

```bash
mkdir -p GO_FDR/separate/CIAPIN1_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/CIAPIN1_wo_original.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/CIAPIN1_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "55 provided genes; 47 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 1138 distinct GO terms"
> [1] "Of those 1138 GO terms, 23 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1_wo_orig/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1_wo_orig/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1_wo_orig/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/CIAPIN1_wo_orig/symbol/default_universe.BP.0.05BH.pvalues.png)

#### Inferred from CYBRD1

```bash
mkdir -p GO_FDR/separate/CYBRD1/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/CYBRD1.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/CYBRD1/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "35 provided genes; 26 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 471 distinct GO terms"
> [1] "Of those 471 GO terms, 16 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/CYBRD1/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1/symbol/default_universe.BP.0.05BH.pvalues.png)

#### Inferred from CYBRD1 w/o CYBRD1

```bash
mkdir -p GO_FDR/separate/CYBRD1_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/CYBRD1_wo_original.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/CYBRD1_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "34 provided genes; 25 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 466 distinct GO terms"
> [1] "Of those 466 GO terms, 17 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/CYBRD1_wo_orig/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1_wo_orig/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1_wo_orig/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/CYBRD1_wo_orig/symbol/default_universe.BP.0.05BH.pvalues.png)

#### Inferred from HAMP

```bash
mkdir -p GO_FDR/separate/HAMP/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/HAMP.list -f 0.2 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/HAMP/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "35 provided genes; 28 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 693 distinct GO terms"
> [1] "Of those 693 GO terms, 0 have a BH-adjusted p-val < 0.2"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

#### Inferred from HAMP w/o HAMP

```bash
mkdir -p GO_FDR/separate/HAMP_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/HAMP_wo_original.list -f 0.2 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/HAMP_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "34 provided genes; 27 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 591 distinct GO terms"
> [1] "Of those 591 GO terms, 0 have a BH-adjusted p-val < 0.2"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

#### Inferred from HFE2

```bash
mkdir -p GO_FDR/separate/HFE2/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/HFE2.list -f 0.2 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/HFE2/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "27 provided genes; 19 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 313 distinct GO terms"
> [1] "Of those 313 GO terms, 2 have a BH-adjusted p-val < 0.2"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

Significant enrichment found for only 2 GO terms. Not sure why it was not enough to compute graphical outputs...? Anyways, output tables we successfully computed.

#### Inferred from HFE2 w/o HFE2

```bash
mkdir -p GO_FDR/separate/HFE2_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/HFE2_wo_original.list -f 0.2 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/HFE2_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "26 provided genes; 19 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 313 distinct GO terms"
> [1] "Of those 313 GO terms, 2 have a BH-adjusted p-val < 0.2"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

Significant enrichment found for only 2 GO terms. Not sure why it was not enough to compute graphical outputs...? Anyways, output tables we successfully computed.

#### Inferred from NEO1

```bash
mkdir -p GO_FDR/separate/NEO1/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/NEO1.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/NEO1/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "39 provided genes; 31 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 942 distinct GO terms"
> [1] "Of those 942 GO terms, 0 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

#### Inferred from NEO1 w/o NEO1

```bash
mkdir -p GO_FDR/separate/NEO1_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/NEO1_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/NEO1_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "38 provided genes; 30 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 933 distinct GO terms"
> [1] "Of those 933 GO terms, 0 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

#### Inferred from SLC11A2

```bash
mkdir -p GO_FDR/separate/SLC11A2/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/SLC11A2.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/SLC11A2/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "68 provided genes; 55 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 1022 distinct GO terms"
> [1] "Of those 1022 GO terms, 7 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/SLC11A2/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from SLC11A2 w/o SLC11A2

```bash
mkdir -p GO_FDR/separate/SLC11A2_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/SLC11A2_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/SLC11A2_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "67 provided genes; 54 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 998 distinct GO terms"
> [1] "Of those 998 GO terms, 7 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/SLC11A2_wo_orig/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2_wo_orig/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2_wo_orig/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC11A2_wo_orig/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from SLC39A14

```bash
mkdir -p GO_FDR/separate/SLC39A14/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/SLC39A14.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/SLC11A2/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "35 provided genes; 26 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 584 distinct GO terms"
> [1] "Of those 584 GO terms, 3 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

#### Inferred from SLC39A14 w/o SLC39A14

```bash
mkdir -p GO_FDR/separate/SLC39A14_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/SLC39A14_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/SLC11A2_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "34 provided genes; 25 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 541 distinct GO terms"
> [1] "Of those 541 GO terms, 2 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```

#### Inferred from SLC40A1

```bash
mkdir -p GO_FDR/separate/SLC40A1/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/SLC40A1.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/SLC40A1/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "26 provided genes; 19 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 401 distinct GO terms"
> [1] "Of those 401 GO terms, 2 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/SLC40A1/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from SLC40A1 w/o SLC40A1

```bash
mkdir -p GO_FDR/separate/SLC40A1_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/SLC40A1_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/SLC40A1_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "25 provided genes; 18 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 381 distinct GO terms"
> [1] "Of those 381 GO terms, 2 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/SLC40A1_wo_orig/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1_wo_orig/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1_wo_orig/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/SLC40A1_wo_orig/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from TFR2

```bash
mkdir -p GO_FDR/separate/TFR2/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/TFR2.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/TFR2/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "34 provided genes; 27 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 386 distinct GO terms"
> [1] "Of those 386 GO terms, 0 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

#### Inferred from TFR2 w/o TFR2

```bash
mkdir -p GO_FDR/separate/TFR2_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/TFR2_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/TFR2_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "33 provided genes; 26 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 360 distinct GO terms"
> [1] "Of those 360 GO terms, 0 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> Error in `[.data.frame`(d, , 2) : undefined columns selected
> Calls: upsetplot ... eval_tidy -> split -> split.default -> [ -> [.data.frame
> Execution halted
> ```

No significant enrichment found in any GO term.

#### Inferred from TFRC

```bash
mkdir -p GO_FDR/separate/TFRC/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/TFRC.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/TFRC/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "10 provided genes; 5 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 172 distinct GO terms"
> [1] "Of those 172 GO terms, 141 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/TFRC/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC/symbol/default_universe.BP.0.1BH.pvalues.png)

#### Inferred from TFRC w/o TFRC

```bash
mkdir -p GO_FDR/separate/TFRC_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/TFRC_wo_original.list -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/TFRC_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "9 provided genes; 4 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 66 distinct GO terms"
> [1] "Of those 66 GO terms, 66 have a BH-adjusted p-val < 0.1"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/TFRC_wo_orig/symbol/default_universe.BP.0.1BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC_wo_orig/symbol/default_universe.BP.0.1BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC_wo_orig/symbol/default_universe.BP.0.1BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC_wo_orig/symbol/default_universe.BP.0.1BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/TFRC_wo_orig/symbol/default_universe.BP.0.1BH.pvalues.png)


#### Inferred from TMPRSS6

```bash
mkdir -p GO_FDR/separate/TMPRSS6/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/TMPRSS6.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/TMPRSS6/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "14 provided genes; 13 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 257 distinct GO terms"
> [1] "Of those 257 GO terms, 8 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6/symbol/default_universe.BP.0.05BH.pvalues.png)


#### Inferred from TMPRSS6 w/o TMPRSS6

```bash
mkdir -p GO_FDR/separate/TMPRSS6_wo_orig/symbol
```

```bash
Rscript GO_enrichment_FDR.R -k "SYMBOL" -G results/separate/TMPRSS6_wo_original.list -f 0.05 -c "BP" -a "BH" -o "default_universe" -d "GO_FDR/separate/TMPRSS6_wo_orig/symbol"
```

> ```R
> ...
> [1] "Loading input data..."
> [1] "Warning: using defaut universe automatically provided by the clusterProfiler package"
> [1] "Done."
> [1] "Computing GO enrichment..."
> `universe` is not in character and will be ignored...
> [1] "Done."
> [1] "18866 (default) background genes"
> [1] "13 provided genes; 12 found by `enrichGO`"
> [1] "Computed GO enrichment (whether significant or not) for 237 distinct GO terms"
> [1] "Of those 237 GO terms, 10 have a BH-adjusted p-val < 0.05"
> [1] "Writing outputs tables..."
> [1] "Done. Writing output images..."
> ...
> [1] "Done."
> ```


![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6_wo_orig/symbol/default_universe.BP.0.05BH.dotplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6_wo_orig/symbol/default_universe.BP.0.05BH.gene-concept.circular.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6_wo_orig/symbol/default_universe.BP.0.05BH.heatplot.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6_wo_orig/symbol/default_universe.BP.0.05BH.upset.png)

![](images_go_enrichment_md/GO_FDR/separate/TMPRSS6_wo_orig/symbol/default_universe.BP.0.05BH.pvalues.png)