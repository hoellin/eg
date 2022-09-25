---
title: "Analysis of E-G Networks starting from genes involved in hemochromatosis / iron metabolism ; based on element-gene relation obtained with CHiC data"
date : "January 2022"
output:
  html_document:
    keep_md: true
    toc: true
    toc_float:
      collapsed: true
    toc_depth : 4
    number_sections : true
    theme: united
    hilight: tango
---


<style type="text/css">
.badCode {
background-color: #C9DDE4;
}
</style>







# Libraries & Version


```{.r .badCode}
library(tidyverse)
library(visNetwork)
```



```
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] visNetwork_2.1.0 forcats_0.5.1    stringr_1.4.0    dplyr_1.0.7     
 [5] purrr_0.3.4      readr_2.0.2      tidyr_1.1.4      tibble_3.1.6    
 [9] ggplot2_3.3.5    tidyverse_1.3.1  knitr_1.36      

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.1  xfun_0.28         bslib_0.3.1       haven_2.4.3      
 [5] colorspace_2.0-2  vctrs_0.3.8       generics_0.1.1    htmltools_0.5.2  
 [9] yaml_2.2.1        utf8_1.2.2        rlang_0.4.12      jquerylib_0.1.4  
[13] pillar_1.6.4      withr_2.4.2       glue_1.5.0        DBI_1.1.1        
[17] dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.1  
[21] munsell_0.5.0     gtable_0.3.0      cellranger_1.1.0  rvest_1.0.2      
[25] htmlwidgets_1.5.4 evaluate_0.14     tzdb_0.2.0        fastmap_1.1.0    
[29] fansi_0.5.0       broom_0.7.10      Rcpp_1.0.7        backports_1.3.0  
[33] scales_1.1.1      formatR_1.11      jsonlite_1.7.2    fs_1.5.0         
[37] hms_1.1.1         digest_0.6.28     stringi_1.7.5     grid_4.0.3       
[41] cli_3.1.0         tools_4.0.3       magrittr_2.0.1    sass_0.4.0       
[45] crayon_1.4.2      pkgconfig_2.0.3   ellipsis_0.3.2    xml2_1.3.2       
[49] reprex_2.0.1      lubridate_1.8.0   rstudioapi_0.13   assertthat_0.2.1 
[53] rmarkdown_2.11    httr_1.4.2        R6_2.5.1          compiler_4.0.3   
```

# Preliminary work

Have a look https://datastorm-open.github.io/visNetwork/nodes.html for documentation on `visNetwork`.

## Data importation

All the files imported here can be found on Genotoul, in `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/` and `/work2/project/regenet/results/phic/homo_sapiens/hg19/jung.ren.2019/LI11/nofilter/pp`.

The full code itself is available here: `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/` (as `preliminary_analysis_v1.Rmd`). Please do not modify directly this repository as it is backed up between multiple computers on a regular basis.


```{.r .badCode}
rm(list = ls())

# wd = 'data/'
wd = "/home/thoellinger/Documents/shared/2022/promoter_capture_hic/data/"
# wd =
# '/home/hoellinger/Documents/INSERM/shared/2022/promoter_capture_hic/data/'
egfile = paste(wd, "liver.all_putative_enhancers.overlapping_ccRE-ELS.sorted.bedpe",
    sep = "")
merged_efile = paste(wd, "list_all_enhancers.overlapping_ccRE-ELS.bed", sep = "")

############# Enhancers #
me = as.data.frame(read.table(merged_efile, sep = "\t"))

############ E-G list #

eg = as.data.frame(read.table(egfile, header = F, col.names = c("chrom1", "start1",
    "end1", "chrom2", "start2", "end2", "name", "score.contact", "strand1", "strand2",
    "tissue", "gene.symbol", "original.distance", "gene.id"), sep = "\t"))

############### known genes #
gene_list = c("HFE", "TFR2", "HFE2", "HAMP", "SLC40A1", "BMP6", "TMPRSS6", "TFRC",
    "SLC11A2", "CYBRD1", "NEO1", "CIAPIN1", "SLC39A14")

# genes known to be causally involved in hemochromatosis: 'HFE', 'TFR2',
# 'HFE2', 'HAMP', 'SLC40A1', 'BMP6' (the other genes are involved in iron
# metabolism regulation) genes co-expressed in the liver: HFE + TFR2 + HJV +
# HAMP + TMPRSS6 genes co-expressed in intestine: DCYTB + DMT1 + SLC40A1
```

## Data preprocessing

### Conversion to factors


```{.r .badCode}
# to_factor_cols = c('chrom1', 'chrom2', 'name', 'strand1', 'strand2',
# 'tissue', 'gene.symbol', 'gene.id')
to_factor_cols = c("chrom1", "chrom2", "name", "strand1", "strand2", "tissue")
eg[to_factor_cols] = lapply(eg[to_factor_cols], factor)
```

## Exploration

### Summary statistics on enhancer lists


```{.r .badCode}
length(me[, 1])
```

```
[1] 25600
```

Warning: we shall pay attention to the fact that those 25,600 putative enhancers correspond to those, among the 31,749 initial putative regulatory elements in the CHiC data, that intersect at least one ccRE-ELS. We filtered out the other elements / the E-G pairs involving such elements. We did so because they are pretty large (see next cell), so it is very unlikely for such large element to be an enhancer but not to intersect any ccRE-ELS.



```{.r .badCode}
me$length = abs(me$V3 - me$V2)

summary(me$length)
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1571    4386    5843    7324    8198  122914 
```

### E-G pairs

We extract the subsample of the E-G bedpe input list, where genes are contained in our list of genes involved either directly in hemochromatosis or in iron metabolism. In the variable name, "dist0" stands for "distance is 0 between the genes in `eg_dist0` and the list of initial genes".


```{.r .badCode}
eg_dist0 = eg[eg$gene.symbol %in% gene_list, ]
```



```{.r .badCode}
length(eg_dist0$name)
```

```
[1] 49
```

```{.r .badCode}
length(unique(eg_dist0$gene.symbol))
```

```
[1] 12
```


```{.r .badCode}
print(as.character(unique(eg_dist0$gene.symbol)))
```

```
 [1] "HFE2"     "CYBRD1"   "SLC40A1"  "TFRC"     "BMP6"     "HFE"     
 [7] "TFR2"     "SLC39A14" "SLC11A2"  "NEO1"     "CIAPIN1"  "TMPRSS6" 
```

```{.r .badCode}
print(gene_list)
```

```
 [1] "HFE"      "TFR2"     "HFE2"     "HAMP"     "SLC40A1"  "BMP6"    
 [7] "TMPRSS6"  "TFRC"     "SLC11A2"  "CYBRD1"   "NEO1"     "CIAPIN1" 
[13] "SLC39A14"
```

We have element-gene pairs in our data for 12 out of 13 initial genes.

### Genes

Chromosomes where the genes are located:


```
     HFE     TFR2     HFE2     HAMP  SLC40A1     BMP6  TMPRSS6     TFRC 
  "chr6"   "chr7"   "chr1"   "HAMP"   "chr2"   "chr6"  "chr22"   "chr3" 
 SLC11A2   CYBRD1     NEO1  CIAPIN1 SLC39A14 
 "chr12"   "chr2"  "chr15"  "chr16"   "chr8" 
```

# Networks

Note: in all subsequent graphs, size of nodes of type "gene" (and not "known_gene", for which the size is fixed) is proportional to the number of distinct enhancers regulating them.

## Find all genes regulated by the initial enhancers

The "initial enhancers" are the enhancers involved in eg_dist0, ie all the enhancers regulating the initial genes in `gene_list`.

`genes` is the subset of `gene_list` for which we have data (12 out of 13 genes here), and `enhancers` is the list of enhancers regulating those initial known genes.


```{.r .badCode}
enhancers = unique(paste(eg_dist0$chrom1, ":", eg_dist0$start1, "-", eg_dist0$end1,
    sep = ""))
genes = unique(eg_dist0$gene.symbol)
```


Now we compute the list of all genes regulated by enhancers in `eg_dist0`. Specifically, we extract, from the full E-G list `eg`, the list `eg_dist1` containing only the enhancers-genes pairs for which the gene G is regulated by any of the enhancers regulating a gene in `gene_list` ("dist1" stands for "distance is at most 1 between the genes in `eg_dist1` and the list of initial genes").



```{.r .badCode}
eg_enhancers_id = data.frame(source = paste(eg$chrom1, ":", eg$start1, "-", eg$end1,
    sep = ""), eg[, -c(1, 2, 3)])  # same as eg but columns 1-3 have been concatenated to make unique enhancers id
eg_dist1 = eg_enhancers_id[eg_enhancers_id$source %in% paste(eg_dist0$chrom1, ":",
    eg_dist0$start1, "-", eg_dist0$end1, sep = ""), ]
eg_dist1$from = lapply(eg_dist1$source, function(x) unique(as.character(eg_dist0[paste(eg_dist0$chrom1,
    ":", eg_dist0$start1, "-", eg_dist0$end1, sep = "") == x, "gene.symbol"])))
eg_dist1$from.id = lapply(eg_dist1$source, function(x) unique(as.character(eg_dist0[paste(eg_dist0$chrom1,
    ":", eg_dist0$start1, "-", eg_dist0$end1, sep = "") == x, "gene.id"])))
eg_dist1$score.contact.IE = left_join(data.frame(name = paste(eg_dist1$source, eg_dist1$from.id,
    eg_dist1$from, sep = "::")), eg_enhancers_id, by = "name")$score.contact  # max ABC score of the I-E pair (initialGene-Enhancer) corresponding to the E-G pair
genes_dist1 = unique(eg_dist1$gene.symbol)
```

`genes_dist1` is the list of genes regulated by `enhancers`.


Compute "contact product", ie the product of the contact scores of:
- the initial gene - enhancer pair (I-E)
- the enhancer - gene pair (E-G)


```{.r .badCode}
eg_dist1$contact.product = eg_dist1$score.contact * eg_dist1$score.contact.IE
eg_dist1$contact.product = eg_dist1$contact.product/max(eg_dist1$contact.product)
```



```{.r .badCode}
print(min(eg_dist1$contact.product))
```

```
[1] 0.1613721
```

```{.r .badCode}
print(mean(eg_dist1$contact.product))
```

```
[1] 0.3154293
```

```{.r .badCode}
print(max(eg_dist1$contact.product))
```

```
[1] 1
```

```{.r .badCode}
print(quantile(eg_dist1$contact.product, c(0.1, 0.4, 0.5, 0.6, 0.8, 0.9)))
```

```
      10%       40%       50%       60%       80%       90% 
0.1853717 0.2340752 0.2457107 0.2883543 0.4243637 0.5197058 
```



```{.r .badCode}
eg_dist1$contact.product.label = 1
eg_dist1[eg_dist1$contact.product >= median(eg_dist1$contact.product), ]$contact.product.label = 2
eg_dist1[eg_dist1$contact.product >= quantile(eg_dist1$contact.product, 0.9)[[1]],
    ]$contact.product.label = 3
```



```{.r .badCode}
table(eg_dist1$contact.product.label)
```

```

 1  2  3 
50 39 12 
```


In the following cell we re-organize `genes_dist1` into `genes_dist1.more` which is in a well-suited format for further "concatenation" with inferences made with other type of data (ABC or QTL -based).



```{.r .badCode}
genes_dist1.more = eg_dist1 %>%
    group_by(gene.symbol) %>%
    mutate(CHiC.sources = paste0(source, collapse = ",")) %>%
    mutate(CHiC.count = length(str_split(CHiC.sources, ",")[[1]])) %>%
    slice(which.max(contact.product.label)) %>%
    ungroup() %>%
    select(gene.symbol, gene.id, CHiC.sources, contact.product.label, CHiC.count,
        from)

genes_dist1.more = subset(genes_dist1.more, !(genes_dist1.more$gene.symbol %in% genes))
genes_dist1.more$from = as.character(genes_dist1.more$from)
genes_dist1.more
```

```
# A tibble: 42 × 6
   gene.symbol gene.id            CHiC.sources contact.product… CHiC.count from 
   <chr>       <chr>              <chr>                   <dbl>      <int> <chr>
 1 ACTL6B      ENSG00000077080.5  chr7:100151…                1          2 TFR2 
 2 ADPGK       ENSG00000159322.13 chr15:73089…                3          2 NEO1 
 3 AGFG2       ENSG00000106351.8  chr7:100151…                1          1 TFR2 
 4 ANKRD34A    ENSG00000181039.7  chr1:145508…                2          1 HFE2 
 5 ARL2BP      ENSG00000102931.3  chr16:57256…                2          1 CIAP…
 6 ATF1        ENSG00000123268.4  chr12:51316…                2          1 SLC1…
 7 C1QTNF6     ENSG00000133466.9  chr22:37461…                2          2 TMPR…
 8 C7orf61     ENSG00000185955.4  chr7:100151…                1          1 TFR2 
 9 CCDC135     ENSG00000159625.10 chr16:57491…                3          2 CIAP…
10 CCL22       ENSG00000102962.4  chr16:57428…                2          1 CIAP…
# … with 32 more rows
```

We re-arrange `eg_dist1` as an edge list `edges_list_dist1`, which will be more suitable to later construct the edges list for visualization as a graph.


```{.r .badCode}
edges_list_dist1 = data.frame(source = eg_dist1$source, target = eg_dist1$gene.symbol,
    distance_kb = floor(eg_dist1$original.distance/1000), inv_dist = 1/(eg_dist1$original.distance +
        1), rescaled_log_inv_dist = 1 - min(log(1/(eg_dist1$original.distance + 1))) +
        log(1/(eg_dist1$original.distance + 1)))
edges_list_dist1
```


Now we can compute the list `nodes_dist1` of (colored) nodes required for our graphs. There are 3 types of nodes: `enhancer`, `known_gene` and (unknown) `gene`.


```{.r .badCode}
nodes_dist1 = full_join(data.frame(label = unique(eg_dist1$source), group = "enhancer"),
    data.frame(label = unique(eg_dist1$gene.symbol), group = "gene")) %>%
    rowid_to_column("id")
nodes_dist1[nodes_dist1$label %in% genes, ]$group = "known_gene"
nodes_dist1 = unique(nodes_dist1)
```




```{.r .badCode}
table(nodes_dist1$group)
```

```

  enhancer       gene known_gene 
        49         42         12 
```

In the list of edges of the graph, `edges_dist1`, the `sample` column indicates in which family of tissues (liver, intestine or both) each E-G pair has been found.




We add to `nodes_dist1` a column `d_in` for plotting purpose: it contains 1 for each node of type `enhancer`, the number of incoming enhancers for each node of type `gene`, and the max of the latter for each node of type `known_gene`.










## Edge weight based on distance

Width proportional to distance. For each enhancer-gene pair $E$-$G$, the distance indicated in the `eg` dataframe as `original_distance.mean` is given in in base pairs. Here, `distance_kb` is the very same quantity but expressed in `kb`.





```{=html}
<div id="htmlwidget-acb42da9b1be409ff14f" style="width:100%;height:400px;" class="visNetwork html-widget"></div>
<script type="application/json" data-for="htmlwidget-acb42da9b1be409ff14f">{"x":{"nodes":{"label":["ACTL6B","ADPGK","AGFG2","ANKRD34A","ARL2BP","ATF1","BMP6","C1QTNF6","C7orf61","CCDC135","CCL22","CELA1","CERS5","chr1:144995955-144999800","chr1:145382770-145387982","chr1:145445045-145450003","chr1:145471844-145478745","chr1:145508505-145515028","chr1:145519579-145527792","chr12:50495719-50503727","chr12:51184932-51188867","chr12:51316042-51321932","chr12:51540821-51574422","chr12:51637790-51642516","chr12:51652001-51656616","chr12:51670781-51676690","chr15:73089434-73099467","chr15:73177101-73180543","chr15:73372310-73380150","chr16:57256998-57262124","chr16:57428584-57434799","chr16:57491859-57510132","chr16:57662892-57669586","chr2:171868841-171871980","chr2:172052913-172056719","chr2:172244003-172248890","chr2:172359194-172365051","chr2:172396460-172401984","chr2:190325771-190329372","chr2:190329373-190334233","chr2:190360995-190365156","chr2:190390528-190393955","chr2:190410593-190416234","chr2:190416235-190421209","chr2:190421210-190425029","chr2:190475143-190478818","chr2:190496245-190499388","chr22:36318565-36323139","chr22:37451969-37461031","chr22:37461032-37471647","chr22:37533344-37544456","chr22:37550713-37556000","chr3:195829327-195839503","chr3:195904445-195910895","chr6:26065333-26073782","chr6:7757920-7762092","chr6:8000333-8004738","chr7:100151130-100158146","chr7:100479857-100484394","chr8:22136614-22143599","chr8:22143600-22148473","chr8:22148474-22155405","CIAPIN1","CLDN15","COQ9","CSRNP2","CX3CL1","CYBRD1","DAZAP2","DOK4","EPO","FBXO24","FIS1","HFE","HFE2","IL2RB","KATNB1","KCTD17","LRCH4","METTL7A","MPST","NEO1","NR4A1","NUDT17","NYAP1","PCOLCE","PDE4DIP","PHYHIP","POLR3D","POLR3GL","RAC2","SAP25","SLC11A2","SLC39A14","SLC40A1","TFR2","TFRC","TMPRSS6","TRIM56","TSC22D4","TST","XPO7","ZNHIT1"],"id":[62,88,63,54,89,84,58,99,60,95,92,86,81,1,2,3,4,5,6,31,32,33,34,35,36,37,38,39,40,41,42,43,44,7,8,9,10,11,12,13,14,15,16,17,18,19,20,45,46,47,48,49,21,22,25,23,24,26,27,28,29,30,90,68,91,83,93,55,80,96,72,67,71,59,50,101,94,98,73,82,102,87,85,52,65,66,51,77,76,53,100,69,79,75,56,64,57,97,70,61,103,78,74],"group":["gene","gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","gene","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","known_gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","known_gene","known_gene","gene","gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","known_gene","known_gene","known_gene","known_gene","known_gene","known_gene","gene","gene","gene","gene","gene"],"d_in":["2","2","1","1","1","1","4","2","1","2","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","4","1","4","1","1","4","2","1","1","1","1","4","4","1","2","2","1","1","1","4","1","1","1","1","1","1","1","1","1","1","4","4","4","4","4","4","1","1","1","1","1"],"value":["2","2","1","1","1","1","4","2","1","2","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","4","1","4","1","1","4","2","1","1","1","1","4","4","1","2","2","1","1","1","4","1","1","1","1","1","1","1","1","1","1","4","4","4","4","4","4","1","1","1","1","1"],"x":[-0.539818751300899,-0.931913348437006,-0.802483689532863,0.042341379387264,0.131248551851264,-0.335719889281041,0.795838934090759,0.777394074358613,-0.741992069607261,0.612493147733946,0.535567835398878,-0.929201045039506,-0.42945760512469,0.250534949073871,0.496612543671098,0.43816573451323,0.31035503702018,0.170649906065797,0.388854514980547,-0.517947641690017,-0.674176452268129,-0.475134862481194,-0.660357120330651,-0.826765922304183,-0.507428915035857,-0.73212114324448,-1,-0.96428528441766,-0.857507540943602,0.24607872993314,0.420488778078205,0.541758050026986,0.453013530154901,0.261230038367272,0.504992470097493,0.357927667157214,0.449145444267993,0.235698453949501,-0.319774391792109,-0.0488209932200246,-0.0873347504653815,-0.147378088033297,-0.390949973292158,-0.255911094461514,-0.200828050812204,-0.0445512521836724,-0.352415607847606,0.687127519864245,0.580916252582676,0.871933029789944,0.828375049429453,0.649412498456972,0.950751305445568,0.983336487908838,0.140813713337155,0.862425816135643,0.730662655775512,-0.644113612725189,-0.344529171020632,0.0563763509191202,-0.12570176956381,-0.0796009229594368,0.399135737651386,-0.338669375801438,0.378819062252798,-0.414132657492003,0.350559219163499,0.361235414295681,-0.578465264438049,0.690364926426057,-0.211080009602393,-0.321083299405402,-0.146438099233644,0.224577497766095,0.333718394723932,1,0.520464098703783,0.912827740684508,-0.205760158513123,-0.715981972222366,0.504472942133186,-0.933221748034676,-0.753994939916711,0.621282204759455,-0.765325670641104,-0.562212262506228,0.243910506191385,-0.18722314471612,-0.0622969167703062,0.0963860089397601,0.928532546621442,-0.507291127983297,-0.630232748807505,-0.0529850241214187,-0.208055896467325,-0.472601167201975,0.957981449874736,0.723691603043618,-0.219370649330408,-0.669939321722994,0.621084763406392,-0.241891792922036,-0.451905434155806],"y":[0.320823614802014,-0.0478880386712512,0.477325399681831,-0.89197832487611,0.397348215427841,-0.376971291288428,0.621642106592727,-0.181708946039531,0.365252212575312,0.174879522170079,0.540958983227488,-0.340883640154498,-0.795994802509153,-0.660218836382031,-0.760783808034329,-0.886599261961266,-0.927258837269484,-0.881858602701351,-0.658148754505659,-0.661786881212473,-0.30238376542911,-0.36760281920543,-0.651374591246716,-0.430142440585986,-0.511362056847935,-0.43683673219542,0.0451386809035841,0.269941198371394,0.0473155361917101,0.319453198088061,0.447792244334034,0.264570805549879,0.183342904125188,0.784393439753486,0.897245291633196,1,0.782506366483784,0.936496973362628,0.745644126694574,0.942061592052742,0.724345924657612,0.978955293875768,0.836801781631213,0.996961329724645,0.69263253534032,0.82756793362754,0.936245392072623,-0.527956547921937,-0.420819876091258,-0.269818681699067,-0.488293929342932,-0.235980782802038,0.334357826677187,0.121782655174708,-0.113441255238103,0.531738115533031,0.713390688224296,0.479709256342257,0.217514428171753,-0.41524969401214,-0.662821662608121,-0.370968157258403,0.330965668277077,0.0458320543521034,0.27226872722711,-0.244074147907702,0.566391925061238,0.878366974564696,-0.765674712393238,0.326565355236694,0.339854198208423,0.371884774818352,0.234866532623506,-0.177708595591489,-0.79399423512255,-0.227961612758284,0.103542078090095,-0.407053440366185,0.152165996609355,-0.174848570447403,-0.236049900145605,0.148616593499717,-0.718300962123788,-0.747597449050371,0.58767474017208,0.599664650694844,-0.528029093286429,-0.81285362281936,-0.758974802777443,-1,-0.117305308026925,0.168916439828084,-0.481548525349525,-0.500393199796389,0.851046706807978,0.387746581989011,0.223981648147316,-0.383720280541913,0.050276448958323,0.650615234365739,-0.108319211522887,-0.670864465943555,0.0732882273309834]},"edges":{"from":[1,1,2,2,3,4,5,5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,28,29,29,29,29,30,31,31,31,32,32,33,33,33,34,34,34,35,35,36,37,38,38,39,40,40,41,41,41,42,42,42,42,43,43,43,43,43,44,44,44,44,45,46,47,47,47,47,47,48,48,49,49,49,49],"to":[50,51,50,52,50,50,53,54,50,50,55,55,55,55,55,56,56,56,56,56,56,56,56,56,57,57,58,58,59,60,61,62,63,64,65,66,67,68,62,69,64,70,71,72,73,74,75,75,76,77,78,75,79,80,81,82,79,79,83,84,79,85,80,79,86,79,79,87,88,87,87,88,89,90,91,90,92,91,93,90,94,95,91,96,95,94,90,91,97,97,98,99,100,101,97,98,97,97,102,103,99],"distance_kb":[411,73,23,197,25,52,37,37,88,99,501,316,124,8,14,104,99,68,39,17,12,8,24,45,10,86,23,265,13,52,52,71,9,71,52,39,292,378,224,292,224,241,378,132,292,378,79,75,36,36,361,68,913,1120,52,127,228,95,150,155,106,841,49,203,92,217,236,242,9,161,20,292,8,214,214,41,8,41,8,7,245,245,7,0,86,86,178,178,1156,18,9,116,146,89,8,81,17,35,131,131,32],"inv_dist":[2.43090753070844e-06,1.35712831648232e-05,4.31276146116358e-05,5.05985812157827e-06,3.96165121622692e-05,1.921561845468e-05,2.6576660376857e-05,2.6576660376857e-05,1.12737029604744e-05,1.00224502886466e-05,1.99513983935134e-06,3.15976731473494e-06,8.04453454323133e-06,0.000122744568552842,7.06114955514758e-05,9.58809542072563e-06,1.00568210388696e-05,1.45959831854274e-05,2.51806713166973e-05,5.73591832052312e-05,8.02632635042941e-05,0.000115754138210441,4.11742907728414e-05,2.2031769812069e-05,9.15834783405074e-05,1.16229064239804e-05,4.3232026285072e-05,3.76585424637725e-06,7.62311327946333e-05,1.89652556516462e-05,1.89652556516462e-05,1.3959461723156e-05,0.000104690117252931,1.3959461723156e-05,1.89652556516462e-05,2.53030034665115e-05,3.42327415633408e-06,2.64316714860414e-06,4.44974458466084e-06,3.42327415633408e-06,4.44974458466084e-06,4.13902145254819e-06,2.64316714860414e-06,7.54142125624995e-06,3.42327415633408e-06,2.64316714860414e-06,1.25161144974154e-05,1.3329245697986e-05,2.7590773645293e-05,2.7590773645293e-05,2.76568558579986e-06,1.46862287233261e-05,1.09451462107357e-06,8.92601316943983e-07,1.88935912938331e-05,7.86311882745172e-06,4.37623355083214e-06,1.04775675279227e-05,6.64332644641825e-06,6.42793322662964e-06,9.40035157314883e-06,1.18827928841083e-06,2.01507274412606e-05,4.91767806912288e-06,1.07965710090475e-05,4.59645429515672e-06,4.23121025306869e-06,4.11978692462026e-06,0.000103092783505155,6.18601342364913e-06,4.77874414603842e-05,3.41791534507273e-06,0.000122819945959224,4.66753171587801e-06,4.66753171587801e-06,2.40552308099396e-05,0.000116103564379426,2.40552308099396e-05,0.000116103564379426,0.000137230684781117,4.06955686595287e-06,4.06955686595287e-06,0.000137230684781117,0.5,1.1591111935368e-05,1.1591111935368e-05,5.60789591745177e-06,5.60789591745177e-06,8.6450224122206e-07,5.30701056095102e-05,0.000110314396028682,8.5652371285899e-06,6.81031899534174e-06,1.11997132873398e-05,0.000121550990640574,1.22884844612114e-05,5.58503211393465e-05,2.83494925440835e-05,7.59659065011623e-06,7.59659065011623e-06,3.08661028458547e-05],"rescaled_log_inv_dist":[2.03386603841288,3.75355740967234,4.90976488354472,2.76693982504792,4.82484738650113,4.10132479089019,4.42563478253509,4.42563478253509,3.56807422336192,3.45042898690226,1.83631552415368,2.29609977166554,3.2305943034736,5.95570689900906,5.40279433910799,3.40612364997709,3.4538524959761,3.8263477479303,4.37167807032702,5.19493433774818,5.53091340702084,5.89706982496082,4.86341543256748,4.23808687539631,5.66285226930564,3.59857922404941,4.9121829511334,2.47157610833937,5.47937132746602,4.08821003658715,4.08821003658715,3.78175891943308,5.79660610368732,3.78175891943308,4.08821003658715,4.37652448413431,2.37619883011053,2.11757925697719,2.6384480791809,2.37619883011053,2.6384480791809,2.56606077713546,2.11757925697719,3.16601204113611,2.37619883011053,2.11757925697719,3.67261835520391,3.7355619271472,4.46308280997775,4.46308280997775,2.16288993705953,3.83251161450613,1.23591237792544,1.03198612996955,4.08442416093634,3.207784706479,2.62178981589923,3.49483792719329,3.03921419024589,3.00625444128465,3.38634847133079,1.31810766592117,4.14884177039974,2.73843786329211,3.52482996590981,2.67088658222654,2.588089445328,2.5614028260543,5.78123077482823,2.96789222562066,5.01236425537871,2.37463219835555,5.95632081025678,2.68623177281226,2.68623177281226,4.32595385155288,5.90008397053053,4.32595385155288,5.90008397053053,6.06726472164942,2.54913549673802,2.54913549673802,6.06726472164942,14.2679647587598,3.59583997330867,3.59583997330867,2.86977697125995,2.86977697125995,1,5.11721516819721,5.84893581613075,3.29331319851475,3.06404034261388,3.56148955998415,5.94593523215078,3.65426398253789,5.16827665688717,4.4902105116894,3.17330092932831,3.17330092932831,4.57525996783962],"width":[4.11,0.73,0.23,1.97,0.25,0.52,0.37,0.37,0.88,0.99,5.01,3.16,1.24,0.08,0.14,1.04,0.99,0.68,0.39,0.17,0.12,0.08,0.24,0.45,0.1,0.86,0.23,2.65,0.13,0.52,0.52,0.71,0.09,0.71,0.52,0.39,2.92,3.78,2.24,2.92,2.24,2.41,3.78,1.32,2.92,3.78,0.79,0.75,0.36,0.36,3.61,0.68,9.13,11.2,0.52,1.27,2.28,0.95,1.5,1.55,1.06,8.41,0.49,2.03,0.92,2.17,2.36,2.42,0.09,1.61,0.2,2.92,0.08,2.14,2.14,0.41,0.08,0.41,0.08,0.07,2.45,2.45,0.07,0,0.86,0.86,1.78,1.78,11.56,0.18,0.09,1.16,1.46,0.89,0.08,0.81,0.17,0.35,1.31,1.31,0.32]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","physics":false},"manipulation":{"enabled":false},"groups":{"enhancer":{"color":"black","shape":"dot"},"gene":{"color":"darkred","shape":"star"},"useDefaultGroups":true,"known_gene":{"color":"green","shape":"star","shadow":{"enabled":true}}},"edges":{"color":{"color":"grey","highlight":"red"},"smooth":false},"physics":{"stabilization":false},"interaction":{"navigationButtons":true,"zoomSpeed":1}},"groups":["gene","known_gene","enhancer"],"width":"100%","height":"400","idselection":{"enabled":false,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":true,"hideColor":"rgba(200,200,200,0.5)","highlight":false,"variable":"sample","main":"Select by sample"},"main":{"text":"E-G Network in Liver and Intestine cells","style":"color:black;font-size:20px;text-align:center;"},"submain":{"text":"Network of genes involved in hemochromatosis and iron metabolism, their enhancers, and other genes regulated by the latter enhancers","style":"color:grey;font-size:14px;text-align:center;"},"footer":null,"background":"rgba(0, 0, 0, 0)","igraphlayout":{"type":"square"},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"nodes":[{"label":"gene","shape":"star","color":"darkred","size":15},{"label":"known\n gene","shape":"star","color":"green","size":15},{"label":"enhancer","shape":"dot","color":"black","size":15}],"nodesToDataframe":false},"highlight":{"enabled":true,"hoverNearest":false,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"}},"evals":[],"jsHooks":[]}</script>
```

Width proportional to inverse distance:





```{=html}
<div id="htmlwidget-2ea43febdeeb29adb7bd" style="width:100%;height:400px;" class="visNetwork html-widget"></div>
<script type="application/json" data-for="htmlwidget-2ea43febdeeb29adb7bd">{"x":{"nodes":{"label":["ACTL6B","ADPGK","AGFG2","ANKRD34A","ARL2BP","ATF1","BMP6","C1QTNF6","C7orf61","CCDC135","CCL22","CELA1","CERS5","chr1:144995955-144999800","chr1:145382770-145387982","chr1:145445045-145450003","chr1:145471844-145478745","chr1:145508505-145515028","chr1:145519579-145527792","chr12:50495719-50503727","chr12:51184932-51188867","chr12:51316042-51321932","chr12:51540821-51574422","chr12:51637790-51642516","chr12:51652001-51656616","chr12:51670781-51676690","chr15:73089434-73099467","chr15:73177101-73180543","chr15:73372310-73380150","chr16:57256998-57262124","chr16:57428584-57434799","chr16:57491859-57510132","chr16:57662892-57669586","chr2:171868841-171871980","chr2:172052913-172056719","chr2:172244003-172248890","chr2:172359194-172365051","chr2:172396460-172401984","chr2:190325771-190329372","chr2:190329373-190334233","chr2:190360995-190365156","chr2:190390528-190393955","chr2:190410593-190416234","chr2:190416235-190421209","chr2:190421210-190425029","chr2:190475143-190478818","chr2:190496245-190499388","chr22:36318565-36323139","chr22:37451969-37461031","chr22:37461032-37471647","chr22:37533344-37544456","chr22:37550713-37556000","chr3:195829327-195839503","chr3:195904445-195910895","chr6:26065333-26073782","chr6:7757920-7762092","chr6:8000333-8004738","chr7:100151130-100158146","chr7:100479857-100484394","chr8:22136614-22143599","chr8:22143600-22148473","chr8:22148474-22155405","CIAPIN1","CLDN15","COQ9","CSRNP2","CX3CL1","CYBRD1","DAZAP2","DOK4","EPO","FBXO24","FIS1","HFE","HFE2","IL2RB","KATNB1","KCTD17","LRCH4","METTL7A","MPST","NEO1","NR4A1","NUDT17","NYAP1","PCOLCE","PDE4DIP","PHYHIP","POLR3D","POLR3GL","RAC2","SAP25","SLC11A2","SLC39A14","SLC40A1","TFR2","TFRC","TMPRSS6","TRIM56","TSC22D4","TST","XPO7","ZNHIT1"],"id":[62,88,63,54,89,84,58,99,60,95,92,86,81,1,2,3,4,5,6,31,32,33,34,35,36,37,38,39,40,41,42,43,44,7,8,9,10,11,12,13,14,15,16,17,18,19,20,45,46,47,48,49,21,22,25,23,24,26,27,28,29,30,90,68,91,83,93,55,80,96,72,67,71,59,50,101,94,98,73,82,102,87,85,52,65,66,51,77,76,53,100,69,79,75,56,64,57,97,70,61,103,78,74],"group":["gene","gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","gene","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","known_gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","known_gene","known_gene","gene","gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","known_gene","known_gene","known_gene","known_gene","known_gene","known_gene","gene","gene","gene","gene","gene"],"d_in":["2","2","1","1","1","1","4","2","1","2","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","4","1","4","1","1","4","2","1","1","1","1","4","4","1","2","2","1","1","1","4","1","1","1","1","1","1","1","1","1","1","4","4","4","4","4","4","1","1","1","1","1"],"value":["2","2","1","1","1","1","4","2","1","2","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","4","1","4","1","1","4","2","1","1","1","1","4","4","1","2","2","1","1","1","4","1","1","1","1","1","1","1","1","1","1","4","4","4","4","4","4","1","1","1","1","1"],"x":[-0.539818751300899,-0.931913348437006,-0.802483689532863,0.042341379387264,0.131248551851264,-0.335719889281041,0.795838934090759,0.777394074358613,-0.741992069607261,0.612493147733946,0.535567835398878,-0.929201045039506,-0.42945760512469,0.250534949073871,0.496612543671098,0.43816573451323,0.31035503702018,0.170649906065797,0.388854514980547,-0.517947641690017,-0.674176452268129,-0.475134862481194,-0.660357120330651,-0.826765922304183,-0.507428915035857,-0.73212114324448,-1,-0.96428528441766,-0.857507540943602,0.24607872993314,0.420488778078205,0.541758050026986,0.453013530154901,0.261230038367272,0.504992470097493,0.357927667157214,0.449145444267993,0.235698453949501,-0.319774391792109,-0.0488209932200246,-0.0873347504653815,-0.147378088033297,-0.390949973292158,-0.255911094461514,-0.200828050812204,-0.0445512521836724,-0.352415607847606,0.687127519864245,0.580916252582676,0.871933029789944,0.828375049429453,0.649412498456972,0.950751305445568,0.983336487908838,0.140813713337155,0.862425816135643,0.730662655775512,-0.644113612725189,-0.344529171020632,0.0563763509191202,-0.12570176956381,-0.0796009229594368,0.399135737651386,-0.338669375801438,0.378819062252798,-0.414132657492003,0.350559219163499,0.361235414295681,-0.578465264438049,0.690364926426057,-0.211080009602393,-0.321083299405402,-0.146438099233644,0.224577497766095,0.333718394723932,1,0.520464098703783,0.912827740684508,-0.205760158513123,-0.715981972222366,0.504472942133186,-0.933221748034676,-0.753994939916711,0.621282204759455,-0.765325670641104,-0.562212262506228,0.243910506191385,-0.18722314471612,-0.0622969167703062,0.0963860089397601,0.928532546621442,-0.507291127983297,-0.630232748807505,-0.0529850241214187,-0.208055896467325,-0.472601167201975,0.957981449874736,0.723691603043618,-0.219370649330408,-0.669939321722994,0.621084763406392,-0.241891792922036,-0.451905434155806],"y":[0.320823614802014,-0.0478880386712512,0.477325399681831,-0.89197832487611,0.397348215427841,-0.376971291288428,0.621642106592727,-0.181708946039531,0.365252212575312,0.174879522170079,0.540958983227488,-0.340883640154498,-0.795994802509153,-0.660218836382031,-0.760783808034329,-0.886599261961266,-0.927258837269484,-0.881858602701351,-0.658148754505659,-0.661786881212473,-0.30238376542911,-0.36760281920543,-0.651374591246716,-0.430142440585986,-0.511362056847935,-0.43683673219542,0.0451386809035841,0.269941198371394,0.0473155361917101,0.319453198088061,0.447792244334034,0.264570805549879,0.183342904125188,0.784393439753486,0.897245291633196,1,0.782506366483784,0.936496973362628,0.745644126694574,0.942061592052742,0.724345924657612,0.978955293875768,0.836801781631213,0.996961329724645,0.69263253534032,0.82756793362754,0.936245392072623,-0.527956547921937,-0.420819876091258,-0.269818681699067,-0.488293929342932,-0.235980782802038,0.334357826677187,0.121782655174708,-0.113441255238103,0.531738115533031,0.713390688224296,0.479709256342257,0.217514428171753,-0.41524969401214,-0.662821662608121,-0.370968157258403,0.330965668277077,0.0458320543521034,0.27226872722711,-0.244074147907702,0.566391925061238,0.878366974564696,-0.765674712393238,0.326565355236694,0.339854198208423,0.371884774818352,0.234866532623506,-0.177708595591489,-0.79399423512255,-0.227961612758284,0.103542078090095,-0.407053440366185,0.152165996609355,-0.174848570447403,-0.236049900145605,0.148616593499717,-0.718300962123788,-0.747597449050371,0.58767474017208,0.599664650694844,-0.528029093286429,-0.81285362281936,-0.758974802777443,-1,-0.117305308026925,0.168916439828084,-0.481548525349525,-0.500393199796389,0.851046706807978,0.387746581989011,0.223981648147316,-0.383720280541913,0.050276448958323,0.650615234365739,-0.108319211522887,-0.670864465943555,0.0732882273309834]},"edges":{"from":[1,1,2,2,3,4,5,5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,28,29,29,29,29,30,31,31,31,32,32,33,33,33,34,34,34,35,35,36,37,38,38,39,40,40,41,41,41,42,42,42,42,43,43,43,43,43,44,44,44,44,45,46,47,47,47,47,47,48,48,49,49,49,49],"to":[50,51,50,52,50,50,53,54,50,50,55,55,55,55,55,56,56,56,56,56,56,56,56,56,57,57,58,58,59,60,61,62,63,64,65,66,67,68,62,69,64,70,71,72,73,74,75,75,76,77,78,75,79,80,81,82,79,79,83,84,79,85,80,79,86,79,79,87,88,87,87,88,89,90,91,90,92,91,93,90,94,95,91,96,95,94,90,91,97,97,98,99,100,101,97,98,97,97,102,103,99],"distance_kb":[411,73,23,197,25,52,37,37,88,99,501,316,124,8,14,104,99,68,39,17,12,8,24,45,10,86,23,265,13,52,52,71,9,71,52,39,292,378,224,292,224,241,378,132,292,378,79,75,36,36,361,68,913,1120,52,127,228,95,150,155,106,841,49,203,92,217,236,242,9,161,20,292,8,214,214,41,8,41,8,7,245,245,7,0,86,86,178,178,1156,18,9,116,146,89,8,81,17,35,131,131,32],"inv_dist":[2.43090753070844e-06,1.35712831648232e-05,4.31276146116358e-05,5.05985812157827e-06,3.96165121622692e-05,1.921561845468e-05,2.6576660376857e-05,2.6576660376857e-05,1.12737029604744e-05,1.00224502886466e-05,1.99513983935134e-06,3.15976731473494e-06,8.04453454323133e-06,0.000122744568552842,7.06114955514758e-05,9.58809542072563e-06,1.00568210388696e-05,1.45959831854274e-05,2.51806713166973e-05,5.73591832052312e-05,8.02632635042941e-05,0.000115754138210441,4.11742907728414e-05,2.2031769812069e-05,9.15834783405074e-05,1.16229064239804e-05,4.3232026285072e-05,3.76585424637725e-06,7.62311327946333e-05,1.89652556516462e-05,1.89652556516462e-05,1.3959461723156e-05,0.000104690117252931,1.3959461723156e-05,1.89652556516462e-05,2.53030034665115e-05,3.42327415633408e-06,2.64316714860414e-06,4.44974458466084e-06,3.42327415633408e-06,4.44974458466084e-06,4.13902145254819e-06,2.64316714860414e-06,7.54142125624995e-06,3.42327415633408e-06,2.64316714860414e-06,1.25161144974154e-05,1.3329245697986e-05,2.7590773645293e-05,2.7590773645293e-05,2.76568558579986e-06,1.46862287233261e-05,1.09451462107357e-06,8.92601316943983e-07,1.88935912938331e-05,7.86311882745172e-06,4.37623355083214e-06,1.04775675279227e-05,6.64332644641825e-06,6.42793322662964e-06,9.40035157314883e-06,1.18827928841083e-06,2.01507274412606e-05,4.91767806912288e-06,1.07965710090475e-05,4.59645429515672e-06,4.23121025306869e-06,4.11978692462026e-06,0.000103092783505155,6.18601342364913e-06,4.77874414603842e-05,3.41791534507273e-06,0.000122819945959224,4.66753171587801e-06,4.66753171587801e-06,2.40552308099396e-05,0.000116103564379426,2.40552308099396e-05,0.000116103564379426,0.000137230684781117,4.06955686595287e-06,4.06955686595287e-06,0.000137230684781117,0.5,1.1591111935368e-05,1.1591111935368e-05,5.60789591745177e-06,5.60789591745177e-06,8.6450224122206e-07,5.30701056095102e-05,0.000110314396028682,8.5652371285899e-06,6.81031899534174e-06,1.11997132873398e-05,0.000121550990640574,1.22884844612114e-05,5.58503211393465e-05,2.83494925440835e-05,7.59659065011623e-06,7.59659065011623e-06,3.08661028458547e-05],"rescaled_log_inv_dist":[2.03386603841288,3.75355740967234,4.90976488354472,2.76693982504792,4.82484738650113,4.10132479089019,4.42563478253509,4.42563478253509,3.56807422336192,3.45042898690226,1.83631552415368,2.29609977166554,3.2305943034736,5.95570689900906,5.40279433910799,3.40612364997709,3.4538524959761,3.8263477479303,4.37167807032702,5.19493433774818,5.53091340702084,5.89706982496082,4.86341543256748,4.23808687539631,5.66285226930564,3.59857922404941,4.9121829511334,2.47157610833937,5.47937132746602,4.08821003658715,4.08821003658715,3.78175891943308,5.79660610368732,3.78175891943308,4.08821003658715,4.37652448413431,2.37619883011053,2.11757925697719,2.6384480791809,2.37619883011053,2.6384480791809,2.56606077713546,2.11757925697719,3.16601204113611,2.37619883011053,2.11757925697719,3.67261835520391,3.7355619271472,4.46308280997775,4.46308280997775,2.16288993705953,3.83251161450613,1.23591237792544,1.03198612996955,4.08442416093634,3.207784706479,2.62178981589923,3.49483792719329,3.03921419024589,3.00625444128465,3.38634847133079,1.31810766592117,4.14884177039974,2.73843786329211,3.52482996590981,2.67088658222654,2.588089445328,2.5614028260543,5.78123077482823,2.96789222562066,5.01236425537871,2.37463219835555,5.95632081025678,2.68623177281226,2.68623177281226,4.32595385155288,5.90008397053053,4.32595385155288,5.90008397053053,6.06726472164942,2.54913549673802,2.54913549673802,6.06726472164942,14.2679647587598,3.59583997330867,3.59583997330867,2.86977697125995,2.86977697125995,1,5.11721516819721,5.84893581613075,3.29331319851475,3.06404034261388,3.56148955998415,5.94593523215078,3.65426398253789,5.16827665688717,4.4902105116894,3.17330092932831,3.17330092932831,4.57525996783962],"width":[4.88633787317293e-05,0.000272794724102501,0.000866902973454641,0.000101707597101171,0.000796326727101369,0.000386250826184984,0.000534214240983676,0.000534214240983676,0.00022661134185805,0.000201460062996039,4.01040649886731e-05,6.3514101237342e-05,0.000161702217439688,0.00246727375052078,0.00141935314542386,0.000192729148246268,0.000202150945295849,0.000293392095479518,0.000506153633457376,0.00115297001522845,0.00161336216754898,0.00232675995433415,0.000827639446843693,0.000442857944556892,0.00184090843900474,0.000233630638510092,0.000869001739894201,7.56969814625552e-05,0.00153231279505205,0.000381218313713639,0.000381218313713639,0.000280597454429236,0.00210436340509765,0.000280597454429236,0.000381218313713639,0.000508612617228632,6.88108204406875e-05,5.31299836797453e-05,8.94437785695529e-05,6.88108204406875e-05,8.94437785695529e-05,8.31979704121752e-05,5.31299836797453e-05,0.000151589198011273,6.88108204406875e-05,5.31299836797453e-05,0.000251584906135309,0.00026792955820872,0.000554598809333759,0.000554598809333759,5.5592711991163e-05,0.000295206110139266,2.2000706230626e-05,1.79420712862588e-05,0.000379777797111033,0.000158055602043568,8.79661421553509e-05,0.000210608319665271,0.000133536702687842,0.000129207109586573,0.000188955331837043,2.38854584862662e-05,0.000405047339005617,9.88496530356471e-05,0.000217020570118252,9.23927727443717e-05,8.50510463592245e-05,8.28113394889518e-05,0.00207225559231884,0.000124344308839768,0.000960569590246237,6.87031036226237e-05,0.0024687889026643,9.38214913953715e-05,9.38214913953715e-05,0.000483531289733053,0.00233378372756215,0.000483531289733053,0.00233378372756215,0.00275845742356152,8.18016711451845e-05,8.18016711451845e-05,0.00275845742356152,10.0504396227464,0.000232991541333821,0.000232991541333821,0.00011272363865799,0.00011272363865799,1.73772551582625e-05,0.00106675578440231,0.002217416353612,0.000172168797230797,0.00013689339974865,0.000225124084372959,0.00244328178503619,0.000247009342264925,0.001122640561044,0.000569849726299619,0.000152698151335426,0.000152698151335426,0.000620435806083485]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","physics":false},"manipulation":{"enabled":false},"groups":{"enhancer":{"color":"black","shape":"dot"},"gene":{"color":"darkred","shape":"star"},"useDefaultGroups":true,"known_gene":{"color":"green","shape":"star","shadow":{"enabled":true}}},"edges":{"color":{"color":"grey","highlight":"red"},"smooth":false},"physics":{"stabilization":false},"interaction":{"navigationButtons":true,"zoomSpeed":1}},"groups":["gene","known_gene","enhancer"],"width":"100%","height":"400","idselection":{"enabled":false,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":true,"hideColor":"rgba(200,200,200,0.5)","highlight":false,"variable":"sample","main":"Select by sample"},"main":{"text":"E-G Network in Liver and Intestine cells","style":"color:black;font-size:20px;text-align:center;"},"submain":{"text":"Network of genes involved in hemochromatosis and iron metabolism, their enhancers, and other genes regulated by the latter enhancers","style":"color:grey;font-size:14px;text-align:center;"},"footer":null,"background":"rgba(0, 0, 0, 0)","igraphlayout":{"type":"square"},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"nodes":[{"label":"gene","shape":"star","color":"darkred","size":15},{"label":"known\n gene","shape":"star","color":"green","size":15},{"label":"enhancer","shape":"dot","color":"black","size":15}],"nodesToDataframe":false},"highlight":{"enabled":true,"hoverNearest":false,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"}},"evals":[],"jsHooks":[]}</script>
```

Width proportional to rescaled and translated log inverse distance:





```{=html}
<div id="htmlwidget-bd1b85b429f7afb2e53f" style="width:100%;height:400px;" class="visNetwork html-widget"></div>
<script type="application/json" data-for="htmlwidget-bd1b85b429f7afb2e53f">{"x":{"nodes":{"label":["ACTL6B","ADPGK","AGFG2","ANKRD34A","ARL2BP","ATF1","BMP6","C1QTNF6","C7orf61","CCDC135","CCL22","CELA1","CERS5","chr1:144995955-144999800","chr1:145382770-145387982","chr1:145445045-145450003","chr1:145471844-145478745","chr1:145508505-145515028","chr1:145519579-145527792","chr12:50495719-50503727","chr12:51184932-51188867","chr12:51316042-51321932","chr12:51540821-51574422","chr12:51637790-51642516","chr12:51652001-51656616","chr12:51670781-51676690","chr15:73089434-73099467","chr15:73177101-73180543","chr15:73372310-73380150","chr16:57256998-57262124","chr16:57428584-57434799","chr16:57491859-57510132","chr16:57662892-57669586","chr2:171868841-171871980","chr2:172052913-172056719","chr2:172244003-172248890","chr2:172359194-172365051","chr2:172396460-172401984","chr2:190325771-190329372","chr2:190329373-190334233","chr2:190360995-190365156","chr2:190390528-190393955","chr2:190410593-190416234","chr2:190416235-190421209","chr2:190421210-190425029","chr2:190475143-190478818","chr2:190496245-190499388","chr22:36318565-36323139","chr22:37451969-37461031","chr22:37461032-37471647","chr22:37533344-37544456","chr22:37550713-37556000","chr3:195829327-195839503","chr3:195904445-195910895","chr6:26065333-26073782","chr6:7757920-7762092","chr6:8000333-8004738","chr7:100151130-100158146","chr7:100479857-100484394","chr8:22136614-22143599","chr8:22143600-22148473","chr8:22148474-22155405","CIAPIN1","CLDN15","COQ9","CSRNP2","CX3CL1","CYBRD1","DAZAP2","DOK4","EPO","FBXO24","FIS1","HFE","HFE2","IL2RB","KATNB1","KCTD17","LRCH4","METTL7A","MPST","NEO1","NR4A1","NUDT17","NYAP1","PCOLCE","PDE4DIP","PHYHIP","POLR3D","POLR3GL","RAC2","SAP25","SLC11A2","SLC39A14","SLC40A1","TFR2","TFRC","TMPRSS6","TRIM56","TSC22D4","TST","XPO7","ZNHIT1"],"id":[62,88,63,54,89,84,58,99,60,95,92,86,81,1,2,3,4,5,6,31,32,33,34,35,36,37,38,39,40,41,42,43,44,7,8,9,10,11,12,13,14,15,16,17,18,19,20,45,46,47,48,49,21,22,25,23,24,26,27,28,29,30,90,68,91,83,93,55,80,96,72,67,71,59,50,101,94,98,73,82,102,87,85,52,65,66,51,77,76,53,100,69,79,75,56,64,57,97,70,61,103,78,74],"group":["gene","gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","gene","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","enhancer","known_gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","known_gene","known_gene","gene","gene","gene","gene","gene","gene","known_gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","known_gene","known_gene","known_gene","known_gene","known_gene","known_gene","gene","gene","gene","gene","gene"],"d_in":["2","2","1","1","1","1","4","2","1","2","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","4","1","4","1","1","4","2","1","1","1","1","4","4","1","2","2","1","1","1","4","1","1","1","1","1","1","1","1","1","1","4","4","4","4","4","4","1","1","1","1","1"],"value":["2","2","1","1","1","1","4","2","1","2","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","4","1","4","1","1","4","2","1","1","1","1","4","4","1","2","2","1","1","1","4","1","1","1","1","1","1","1","1","1","1","4","4","4","4","4","4","1","1","1","1","1"],"x":[-0.539818751300899,-0.931913348437006,-0.802483689532863,0.042341379387264,0.131248551851264,-0.335719889281041,0.795838934090759,0.777394074358613,-0.741992069607261,0.612493147733946,0.535567835398878,-0.929201045039506,-0.42945760512469,0.250534949073871,0.496612543671098,0.43816573451323,0.31035503702018,0.170649906065797,0.388854514980547,-0.517947641690017,-0.674176452268129,-0.475134862481194,-0.660357120330651,-0.826765922304183,-0.507428915035857,-0.73212114324448,-1,-0.96428528441766,-0.857507540943602,0.24607872993314,0.420488778078205,0.541758050026986,0.453013530154901,0.261230038367272,0.504992470097493,0.357927667157214,0.449145444267993,0.235698453949501,-0.319774391792109,-0.0488209932200246,-0.0873347504653815,-0.147378088033297,-0.390949973292158,-0.255911094461514,-0.200828050812204,-0.0445512521836724,-0.352415607847606,0.687127519864245,0.580916252582676,0.871933029789944,0.828375049429453,0.649412498456972,0.950751305445568,0.983336487908838,0.140813713337155,0.862425816135643,0.730662655775512,-0.644113612725189,-0.344529171020632,0.0563763509191202,-0.12570176956381,-0.0796009229594368,0.399135737651386,-0.338669375801438,0.378819062252798,-0.414132657492003,0.350559219163499,0.361235414295681,-0.578465264438049,0.690364926426057,-0.211080009602393,-0.321083299405402,-0.146438099233644,0.224577497766095,0.333718394723932,1,0.520464098703783,0.912827740684508,-0.205760158513123,-0.715981972222366,0.504472942133186,-0.933221748034676,-0.753994939916711,0.621282204759455,-0.765325670641104,-0.562212262506228,0.243910506191385,-0.18722314471612,-0.0622969167703062,0.0963860089397601,0.928532546621442,-0.507291127983297,-0.630232748807505,-0.0529850241214187,-0.208055896467325,-0.472601167201975,0.957981449874736,0.723691603043618,-0.219370649330408,-0.669939321722994,0.621084763406392,-0.241891792922036,-0.451905434155806],"y":[0.320823614802014,-0.0478880386712512,0.477325399681831,-0.89197832487611,0.397348215427841,-0.376971291288428,0.621642106592727,-0.181708946039531,0.365252212575312,0.174879522170079,0.540958983227488,-0.340883640154498,-0.795994802509153,-0.660218836382031,-0.760783808034329,-0.886599261961266,-0.927258837269484,-0.881858602701351,-0.658148754505659,-0.661786881212473,-0.30238376542911,-0.36760281920543,-0.651374591246716,-0.430142440585986,-0.511362056847935,-0.43683673219542,0.0451386809035841,0.269941198371394,0.0473155361917101,0.319453198088061,0.447792244334034,0.264570805549879,0.183342904125188,0.784393439753486,0.897245291633196,1,0.782506366483784,0.936496973362628,0.745644126694574,0.942061592052742,0.724345924657612,0.978955293875768,0.836801781631213,0.996961329724645,0.69263253534032,0.82756793362754,0.936245392072623,-0.527956547921937,-0.420819876091258,-0.269818681699067,-0.488293929342932,-0.235980782802038,0.334357826677187,0.121782655174708,-0.113441255238103,0.531738115533031,0.713390688224296,0.479709256342257,0.217514428171753,-0.41524969401214,-0.662821662608121,-0.370968157258403,0.330965668277077,0.0458320543521034,0.27226872722711,-0.244074147907702,0.566391925061238,0.878366974564696,-0.765674712393238,0.326565355236694,0.339854198208423,0.371884774818352,0.234866532623506,-0.177708595591489,-0.79399423512255,-0.227961612758284,0.103542078090095,-0.407053440366185,0.152165996609355,-0.174848570447403,-0.236049900145605,0.148616593499717,-0.718300962123788,-0.747597449050371,0.58767474017208,0.599664650694844,-0.528029093286429,-0.81285362281936,-0.758974802777443,-1,-0.117305308026925,0.168916439828084,-0.481548525349525,-0.500393199796389,0.851046706807978,0.387746581989011,0.223981648147316,-0.383720280541913,0.050276448958323,0.650615234365739,-0.108319211522887,-0.670864465943555,0.0732882273309834]},"edges":{"from":[1,1,2,2,3,4,5,5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,28,29,29,29,29,30,31,31,31,32,32,33,33,33,34,34,34,35,35,36,37,38,38,39,40,40,41,41,41,42,42,42,42,43,43,43,43,43,44,44,44,44,45,46,47,47,47,47,47,48,48,49,49,49,49],"to":[50,51,50,52,50,50,53,54,50,50,55,55,55,55,55,56,56,56,56,56,56,56,56,56,57,57,58,58,59,60,61,62,63,64,65,66,67,68,62,69,64,70,71,72,73,74,75,75,76,77,78,75,79,80,81,82,79,79,83,84,79,85,80,79,86,79,79,87,88,87,87,88,89,90,91,90,92,91,93,90,94,95,91,96,95,94,90,91,97,97,98,99,100,101,97,98,97,97,102,103,99],"distance_kb":[411,73,23,197,25,52,37,37,88,99,501,316,124,8,14,104,99,68,39,17,12,8,24,45,10,86,23,265,13,52,52,71,9,71,52,39,292,378,224,292,224,241,378,132,292,378,79,75,36,36,361,68,913,1120,52,127,228,95,150,155,106,841,49,203,92,217,236,242,9,161,20,292,8,214,214,41,8,41,8,7,245,245,7,0,86,86,178,178,1156,18,9,116,146,89,8,81,17,35,131,131,32],"inv_dist":[2.43090753070844e-06,1.35712831648232e-05,4.31276146116358e-05,5.05985812157827e-06,3.96165121622692e-05,1.921561845468e-05,2.6576660376857e-05,2.6576660376857e-05,1.12737029604744e-05,1.00224502886466e-05,1.99513983935134e-06,3.15976731473494e-06,8.04453454323133e-06,0.000122744568552842,7.06114955514758e-05,9.58809542072563e-06,1.00568210388696e-05,1.45959831854274e-05,2.51806713166973e-05,5.73591832052312e-05,8.02632635042941e-05,0.000115754138210441,4.11742907728414e-05,2.2031769812069e-05,9.15834783405074e-05,1.16229064239804e-05,4.3232026285072e-05,3.76585424637725e-06,7.62311327946333e-05,1.89652556516462e-05,1.89652556516462e-05,1.3959461723156e-05,0.000104690117252931,1.3959461723156e-05,1.89652556516462e-05,2.53030034665115e-05,3.42327415633408e-06,2.64316714860414e-06,4.44974458466084e-06,3.42327415633408e-06,4.44974458466084e-06,4.13902145254819e-06,2.64316714860414e-06,7.54142125624995e-06,3.42327415633408e-06,2.64316714860414e-06,1.25161144974154e-05,1.3329245697986e-05,2.7590773645293e-05,2.7590773645293e-05,2.76568558579986e-06,1.46862287233261e-05,1.09451462107357e-06,8.92601316943983e-07,1.88935912938331e-05,7.86311882745172e-06,4.37623355083214e-06,1.04775675279227e-05,6.64332644641825e-06,6.42793322662964e-06,9.40035157314883e-06,1.18827928841083e-06,2.01507274412606e-05,4.91767806912288e-06,1.07965710090475e-05,4.59645429515672e-06,4.23121025306869e-06,4.11978692462026e-06,0.000103092783505155,6.18601342364913e-06,4.77874414603842e-05,3.41791534507273e-06,0.000122819945959224,4.66753171587801e-06,4.66753171587801e-06,2.40552308099396e-05,0.000116103564379426,2.40552308099396e-05,0.000116103564379426,0.000137230684781117,4.06955686595287e-06,4.06955686595287e-06,0.000137230684781117,0.5,1.1591111935368e-05,1.1591111935368e-05,5.60789591745177e-06,5.60789591745177e-06,8.6450224122206e-07,5.30701056095102e-05,0.000110314396028682,8.5652371285899e-06,6.81031899534174e-06,1.11997132873398e-05,0.000121550990640574,1.22884844612114e-05,5.58503211393465e-05,2.83494925440835e-05,7.59659065011623e-06,7.59659065011623e-06,3.08661028458547e-05],"rescaled_log_inv_dist":[2.03386603841288,3.75355740967234,4.90976488354472,2.76693982504792,4.82484738650113,4.10132479089019,4.42563478253509,4.42563478253509,3.56807422336192,3.45042898690226,1.83631552415368,2.29609977166554,3.2305943034736,5.95570689900906,5.40279433910799,3.40612364997709,3.4538524959761,3.8263477479303,4.37167807032702,5.19493433774818,5.53091340702084,5.89706982496082,4.86341543256748,4.23808687539631,5.66285226930564,3.59857922404941,4.9121829511334,2.47157610833937,5.47937132746602,4.08821003658715,4.08821003658715,3.78175891943308,5.79660610368732,3.78175891943308,4.08821003658715,4.37652448413431,2.37619883011053,2.11757925697719,2.6384480791809,2.37619883011053,2.6384480791809,2.56606077713546,2.11757925697719,3.16601204113611,2.37619883011053,2.11757925697719,3.67261835520391,3.7355619271472,4.46308280997775,4.46308280997775,2.16288993705953,3.83251161450613,1.23591237792544,1.03198612996955,4.08442416093634,3.207784706479,2.62178981589923,3.49483792719329,3.03921419024589,3.00625444128465,3.38634847133079,1.31810766592117,4.14884177039974,2.73843786329211,3.52482996590981,2.67088658222654,2.588089445328,2.5614028260543,5.78123077482823,2.96789222562066,5.01236425537871,2.37463219835555,5.95632081025678,2.68623177281226,2.68623177281226,4.32595385155288,5.90008397053053,4.32595385155288,5.90008397053053,6.06726472164942,2.54913549673802,2.54913549673802,6.06726472164942,14.2679647587598,3.59583997330867,3.59583997330867,2.86977697125995,2.86977697125995,1,5.11721516819721,5.84893581613075,3.29331319851475,3.06404034261388,3.56148955998415,5.94593523215078,3.65426398253789,5.16827665688717,4.4902105116894,3.17330092932831,3.17330092932831,4.57525996783962],"width":[2.03386603841288,3.75355740967234,4.90976488354472,2.76693982504792,4.82484738650113,4.10132479089019,4.42563478253509,4.42563478253509,3.56807422336192,3.45042898690226,1.83631552415368,2.29609977166554,3.2305943034736,5.95570689900906,5.40279433910799,3.40612364997709,3.4538524959761,3.8263477479303,4.37167807032702,5.19493433774818,5.53091340702084,5.89706982496082,4.86341543256748,4.23808687539631,5.66285226930564,3.59857922404941,4.9121829511334,2.47157610833937,5.47937132746602,4.08821003658715,4.08821003658715,3.78175891943308,5.79660610368732,3.78175891943308,4.08821003658715,4.37652448413431,2.37619883011053,2.11757925697719,2.6384480791809,2.37619883011053,2.6384480791809,2.56606077713546,2.11757925697719,3.16601204113611,2.37619883011053,2.11757925697719,3.67261835520391,3.7355619271472,4.46308280997775,4.46308280997775,2.16288993705953,3.83251161450613,1.23591237792544,1.03198612996955,4.08442416093634,3.207784706479,2.62178981589923,3.49483792719329,3.03921419024589,3.00625444128465,3.38634847133079,1.31810766592117,4.14884177039974,2.73843786329211,3.52482996590981,2.67088658222654,2.588089445328,2.5614028260543,5.78123077482823,2.96789222562066,5.01236425537871,2.37463219835555,5.95632081025678,2.68623177281226,2.68623177281226,4.32595385155288,5.90008397053053,4.32595385155288,5.90008397053053,6.06726472164942,2.54913549673802,2.54913549673802,6.06726472164942,14.2679647587598,3.59583997330867,3.59583997330867,2.86977697125995,2.86977697125995,1,5.11721516819721,5.84893581613075,3.29331319851475,3.06404034261388,3.56148955998415,5.94593523215078,3.65426398253789,5.16827665688717,4.4902105116894,3.17330092932831,3.17330092932831,4.57525996783962]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot","physics":false},"manipulation":{"enabled":false},"groups":{"enhancer":{"color":"black","shape":"dot"},"gene":{"color":"darkred","shape":"star"},"useDefaultGroups":true,"known_gene":{"color":"green","shape":"star","shadow":{"enabled":true}}},"edges":{"color":{"color":"grey","highlight":"red"},"smooth":false},"physics":{"stabilization":false},"interaction":{"navigationButtons":true,"zoomSpeed":1}},"groups":["gene","known_gene","enhancer"],"width":"100%","height":"400","idselection":{"enabled":false,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":true,"hideColor":"rgba(200,200,200,0.5)","highlight":false,"variable":"sample","main":"Select by sample"},"main":{"text":"E-G Network in Liver and Intestine cells","style":"color:black;font-size:20px;text-align:center;"},"submain":{"text":"Network of genes involved in hemochromatosis and iron metabolism, their enhancers, and other genes regulated by the latter enhancers","style":"color:grey;font-size:14px;text-align:center;"},"footer":null,"background":"rgba(0, 0, 0, 0)","igraphlayout":{"type":"square"},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"nodes":[{"label":"gene","shape":"star","color":"darkred","size":15},{"label":"known\n gene","shape":"star","color":"green","size":15},{"label":"enhancer","shape":"dot","color":"black","size":15}],"nodesToDataframe":false},"highlight":{"enabled":true,"hoverNearest":false,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"}},"evals":[],"jsHooks":[]}</script>
```



# Results

The list of the initial genes + all genes regulated by enhancers regulating the 12 initial genes (12 found among 13), can be found here: `/work2/project/regenet/workspace/thoellinger/shared/automne_2021/promoter_capture_hic/results/new_genes_v1.list`.


```{.r .badCode}
# write.table(genes_dist1, file='results/new_genes_v1.list', quote=FALSE,
# sep='\t', row.names=F, col.names=F)
```

The list all genes regulated by enhancers regulating the 12 initial genes, with useful info, can be found here: `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/results/new_genes_chic_v2_more_info.list`.


```{.r .badCode}
# write.table(genes_dist1.more,
# file='results/new_genes_chic_v2_more_info.list', quote=FALSE, sep='\t',
# row.names=F, col.names=T)
```


For each one of the 12 initial genes, the list of genes regulated by enhancers regulating that initial gene, can be found here:


```{.r .badCode}
# for (gene in genes){ current_gene_dist1 =
# unique(eg_dist1[eg_dist1$from==gene,'gene.symbol'])
# write.table(current_gene_dist1, file=paste('results/separate/',gene,'.list',
# sep=''), quote=FALSE, sep='\t', row.names=F, col.names=F) }
```


The list of the enhancers regulating the 12 initial genes, can be found here: `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/results/enhancers.list`


```{.r .badCode}
# write.table(enhancers, file='results/enhancers.list', quote=FALSE, sep='\t',
# row.names=F, col.names=F)
```

The list of E-G pairs involving the 12 initial genes, can be found here: `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/results/eg_dist0.bedpe`


```{.r .badCode}
# write.table(eg_dist0, file='results/eg_dist0.bedpe', quote=FALSE, sep='\t',
# row.names=F, col.names=F)
```

The list of the enhancers regulating the 54 initial+new genes, can be found here: `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/results/enhancers.new_genes.list`


```{.r .badCode}
# eg_dist2 = eg[eg$gene.symbol %in% genes_dist1,] enhancers.new_genes =
# unique(paste(eg_dist2$chrom1, ':', eg_dist2$start1, '-', eg_dist2$end1,
# sep='')) write.table(enhancers.new_genes,
# file='results/enhancers.new_genes.list', quote=FALSE, sep='\t', row.names=F,
# col.names=F)
```

The list of E-G pairs involving the 54 initial+new genes, can be found here: `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/results/eg_dist2.bedpe`


```{.r .badCode}
# write.table(eg_dist2, file='results/eg_dist2.bedpe', quote=FALSE, sep='\t',
# row.names=F, col.names=F)
```


# Appendix

## Tables

### Table of the initial ABC-predicted E-G pairs involving the initial genes



```{.r .badCode}
options(max.print = 1000)
eg_dist0
```

```
      chrom1    start1      end1 chrom2    start2      end2
1602    chr1 144995955 144999800   chr1 145413094 145417545
1613    chr1 145382770 145387982   chr1 145413094 145417545
1620    chr1 145445045 145450003   chr1 145413094 145417545
1622    chr1 145471844 145478745   chr1 145413094 145417545
1626    chr1 145508505 145515028   chr1 145413094 145417545
1627    chr1 145519579 145527792   chr1 145413094 145417545
3982    chr2 171868841 171871980   chr2 172378756 172414643
3983    chr2 172052913 172056719   chr2 172378756 172414643
3985    chr2 172244003 172248890   chr2 172378756 172414643
3986    chr2 172359194 172365051   chr2 172378756 172414643
3987    chr2 172396460 172401984   chr2 172378756 172414643
4074    chr2 190325771 190329372   chr2 190425304 190448484
4075    chr2 190329373 190334233   chr2 190425304 190448484
4076    chr2 190360995 190365156   chr2 190425304 190448484
4077    chr2 190390528 190393955   chr2 190425304 190448484
4078    chr2 190410593 190416234   chr2 190425304 190448484
4079    chr2 190416235 190421209   chr2 190425304 190448484
4080    chr2 190421210 190425029   chr2 190425304 190448484
4081    chr2 190475143 190478818   chr2 190425304 190448484
4085    chr2 190496245 190499388   chr2 190425304 190448484
6360    chr3 195829327 195839503   chr3 195754053 195809060
6361    chr3 195904445 195910895   chr3 195754053 195809060
9222    chr6   7757920   7762092   chr6   7727029   7881655
9223    chr6   8000333   8004738   chr6   7727029   7881655
9356    chr6  26065333  26073782   chr6  26087508  26098571
12009   chr7 100151130 100158146   chr7 100218038 100240402
12040   chr7 100479857 100484394   chr7 100218038 100240402
12830   chr8  22136614  22143599   chr8  22224761  22291642
12831   chr8  22143600  22148473   chr8  22224761  22291642
12835   chr8  22148474  22155405   chr8  22224761  22291642
20296  chr12  50495719  50503727  chr12  51373183  51422349
20338  chr12  51184932  51188867  chr12  51373183  51422349
20341  chr12  51316042  51321932  chr12  51373183  51422349
20362  chr12  51540821  51574422  chr12  51373183  51422349
20367  chr12  51637790  51642516  chr12  51373183  51422349
20369  chr12  51652001  51656616  chr12  51373183  51422349
20373  chr12  51670781  51676690  chr12  51373183  51422349
25072  chr15  73089434  73099467  chr15  73344050  73597547
25074  chr15  73177101  73180543  chr15  73344050  73597547
25078  chr15  73372310  73380150  chr15  73344050  73597547
26940  chr16  57256998  57262124  chr16  57462080  57481440
26952  chr16  57428584  57434799  chr16  57462080  57481440
26967  chr16  57491859  57510132  chr16  57462080  57481440
26985  chr16  57662892  57669586  chr16  57462080  57481440
38401  chr22  36318565  36323139  chr22  37461475  37505603
38476  chr22  37451969  37461031  chr22  37461475  37505603
38481  chr22  37461032  37471647  chr22  37461475  37505603
38484  chr22  37533344  37544456  chr22  37461475  37505603
38485  chr22  37550713  37556000  chr22  37461475  37505603
                                                      name score.contact
1602    chr1:144995955-144999800::ENSG00000168509.13::HFE2      2.691613
1613    chr1:145382770-145387982::ENSG00000168509.13::HFE2      3.500083
1620    chr1:145445045-145450003::ENSG00000168509.13::HFE2      2.747677
1622    chr1:145471844-145478745::ENSG00000168509.13::HFE2      2.173732
1626    chr1:145508505-145515028::ENSG00000168509.13::HFE2      3.404033
1627    chr1:145519579-145527792::ENSG00000168509.13::HFE2      3.329819
3982   chr2:171868841-171871980::ENSG00000071967.7::CYBRD1      3.161168
3983   chr2:172052913-172056719::ENSG00000071967.7::CYBRD1      2.180688
3985   chr2:172244003-172248890::ENSG00000071967.7::CYBRD1      2.113738
3986   chr2:172359194-172365051::ENSG00000071967.7::CYBRD1      2.961889
3987   chr2:172396460-172401984::ENSG00000071967.7::CYBRD1      2.013556
4074  chr2:190325771-190329372::ENSG00000138449.6::SLC40A1      4.556188
4075  chr2:190329373-190334233::ENSG00000138449.6::SLC40A1      2.672282
4076  chr2:190360995-190365156::ENSG00000138449.6::SLC40A1      3.338258
4077  chr2:190390528-190393955::ENSG00000138449.6::SLC40A1      2.310444
4078  chr2:190410593-190416234::ENSG00000138449.6::SLC40A1      2.833319
4079  chr2:190416235-190421209::ENSG00000138449.6::SLC40A1      2.538705
4080  chr2:190421210-190425029::ENSG00000138449.6::SLC40A1      2.206945
4081  chr2:190475143-190478818::ENSG00000138449.6::SLC40A1      2.561220
4085  chr2:190496245-190499388::ENSG00000138449.6::SLC40A1      2.306987
6360     chr3:195829327-195839503::ENSG00000072274.8::TFRC      2.880979
6361     chr3:195904445-195910895::ENSG00000072274.8::TFRC      2.461830
9222         chr6:7757920-7762092::ENSG00000153162.8::BMP6      3.182610
9223         chr6:8000333-8004738::ENSG00000153162.8::BMP6      2.373608
9356       chr6:26065333-26073782::ENSG00000010704.14::HFE      2.104796
12009    chr7:100151130-100158146::ENSG00000106327.8::TFR2      2.202777
12040    chr7:100479857-100484394::ENSG00000106327.8::TFR2      2.085438
12830  chr8:22136614-22143599::ENSG00000104635.9::SLC39A14      2.941420
12831  chr8:22143600-22148473::ENSG00000104635.9::SLC39A14      2.234129
12835  chr8:22148474-22155405::ENSG00000104635.9::SLC39A14      2.240686
20296 chr12:50495719-50503727::ENSG00000110911.10::SLC11A2      2.207329
20338 chr12:51184932-51188867::ENSG00000110911.10::SLC11A2      2.739063
20341 chr12:51316042-51321932::ENSG00000110911.10::SLC11A2      2.933383
20362 chr12:51540821-51574422::ENSG00000110911.10::SLC11A2      2.623970
20367 chr12:51637790-51642516::ENSG00000110911.10::SLC11A2      2.607919
20369 chr12:51652001-51656616::ENSG00000110911.10::SLC11A2      2.923081
20373 chr12:51670781-51676690::ENSG00000110911.10::SLC11A2      2.106838
25072    chr15:73089434-73099467::ENSG00000067141.12::NEO1      3.560153
25074    chr15:73177101-73180543::ENSG00000067141.12::NEO1      2.328833
25078    chr15:73372310-73380150::ENSG00000067141.12::NEO1      3.416054
26940 chr16:57256998-57262124::ENSG00000005194.10::CIAPIN1      2.036134
26952 chr16:57428584-57434799::ENSG00000005194.10::CIAPIN1      5.012442
26967 chr16:57491859-57510132::ENSG00000005194.10::CIAPIN1      2.158098
26985 chr16:57662892-57669586::ENSG00000005194.10::CIAPIN1      3.803306
38401 chr22:36318565-36323139::ENSG00000187045.12::TMPRSS6      2.429089
38476 chr22:37451969-37461031::ENSG00000187045.12::TMPRSS6      2.256193
38481 chr22:37461032-37471647::ENSG00000187045.12::TMPRSS6      3.763882
38484 chr22:37533344-37544456::ENSG00000187045.12::TMPRSS6      2.021822
38485 chr22:37550713-37556000::ENSG00000187045.12::TMPRSS6      2.426097
      strand1 strand2 tissue gene.symbol original.distance            gene.id
1602        .       .  liver        HFE2            411368 ENSG00000168509.13
1613        .       .  liver        HFE2             23186 ENSG00000168509.13
1620        .       .  liver        HFE2             25241 ENSG00000168509.13
1622        .       .  liver        HFE2             52040 ENSG00000168509.13
1626        .       .  liver        HFE2             88701 ENSG00000168509.13
1627        .       .  liver        HFE2             99775 ENSG00000168509.13
3982        .       .  liver      CYBRD1            501217  ENSG00000071967.7
3983        .       .  liver      CYBRD1            316478  ENSG00000071967.7
3985        .       .  liver      CYBRD1            124307  ENSG00000071967.7
3986        .       .  liver      CYBRD1              8146  ENSG00000071967.7
3987        .       .  liver      CYBRD1             14161  ENSG00000071967.7
4074        .       .  liver     SLC40A1            104295  ENSG00000138449.6
4075        .       .  liver     SLC40A1             99434  ENSG00000138449.6
4076        .       .  liver     SLC40A1             68511  ENSG00000138449.6
4077        .       .  liver     SLC40A1             39712  ENSG00000138449.6
4078        .       .  liver     SLC40A1             17433  ENSG00000138449.6
4079        .       .  liver     SLC40A1             12458  ENSG00000138449.6
4080        .       .  liver     SLC40A1              8638  ENSG00000138449.6
4081        .       .  liver     SLC40A1             24286  ENSG00000138449.6
4085        .       .  liver     SLC40A1             45388  ENSG00000138449.6
6360        .       .  liver        TFRC             10918  ENSG00000072274.8
6361        .       .  liver        TFRC             86036  ENSG00000072274.8
9222        .       .  liver        BMP6             23130  ENSG00000153162.8
9223        .       .  liver        BMP6            265543  ENSG00000153162.8
9356        .       .  liver         HFE             13117 ENSG00000010704.14
12009       .       .  liver        TFR2             71635  ENSG00000106327.8
12040       .       .  liver        TFR2            224731  ENSG00000106327.8
12830       .       .  liver    SLC39A14             79896  ENSG00000104635.9
12831       .       .  liver    SLC39A14             75022  ENSG00000104635.9
12835       .       .  liver    SLC39A14             68090  ENSG00000104635.9
20296       .       .  liver     SLC11A2            913646 ENSG00000110911.10
20338       .       .  liver     SLC11A2            228506 ENSG00000110911.10
20341       .       .  liver     SLC11A2             95441 ENSG00000110911.10
20362       .       .  liver     SLC11A2            106378 ENSG00000110911.10
20367       .       .  liver     SLC11A2            203347 ENSG00000110911.10
20369       .       .  liver     SLC11A2            217558 ENSG00000110911.10
20373       .       .  liver     SLC11A2            236338 ENSG00000110911.10
25072       .       .  liver        NEO1            242730 ENSG00000067141.12
25074       .       .  liver        NEO1            161654 ENSG00000067141.12
25078       .       .  liver        NEO1             20925 ENSG00000067141.12
26940       .       .  liver     CIAPIN1            214245 ENSG00000005194.10
26952       .       .  liver     CIAPIN1             41570 ENSG00000005194.10
26967       .       .  liver     CIAPIN1              7286 ENSG00000005194.10
26985       .       .  liver     CIAPIN1            178319 ENSG00000005194.10
38401       .       .  liver     TMPRSS6           1156734 ENSG00000187045.12
38476       .       .  liver     TMPRSS6             18842 ENSG00000187045.12
38481       .       .  liver     TMPRSS6              8226 ENSG00000187045.12
38484       .       .  liver     TMPRSS6             17904 ENSG00000187045.12
38485       .       .  liver     TMPRSS6             35273 ENSG00000187045.12
```

```{.r .badCode}
options(max.print = 75)
```


### List of the infered new genes


```{.r .badCode}
options(max.print = 1000)
genes_dist1
```

```
 [1] "HFE2"     "PDE4DIP"  "NUDT17"   "POLR3GL"  "ANKRD34A" "CYBRD1"  
 [7] "SLC40A1"  "TFRC"     "BMP6"     "HFE"      "C7orf61"  "TSC22D4" 
[13] "ACTL6B"   "AGFG2"    "TFR2"     "NYAP1"    "PCOLCE"   "FBXO24"  
[19] "CLDN15"   "SAP25"    "TRIM56"   "FIS1"     "EPO"      "LRCH4"   
[25] "ZNHIT1"   "SLC39A14" "POLR3D"   "PHYHIP"   "XPO7"     "SLC11A2" 
[31] "DAZAP2"   "CERS5"    "METTL7A"  "CSRNP2"   "ATF1"     "NR4A1"   
[37] "CELA1"    "NEO1"     "ADPGK"    "ARL2BP"   "CIAPIN1"  "COQ9"    
[43] "CCL22"    "CX3CL1"   "KATNB1"   "CCDC135"  "DOK4"     "TMPRSS6" 
[49] "KCTD17"   "C1QTNF6"  "RAC2"     "IL2RB"    "MPST"     "TST"     
```

```{.r .badCode}
length(genes_dist1)
```

```
[1] 54
```

```{.r .badCode}
options(max.print = 75)
```

