## Context

If GWAS have identified hundreds of variants associated to a broad range of traits and diseases, little insight has been gained about the biological mechanisms underlying these associations. Part of this gap in knowledge is likely explained by the fact that the large majority of these variants do not impact protein-coding regions but are rather located in regulatory regions active in disease-related cell types, and in particular in enhancers.

Therefore it is crucial to be able to better identify E-G relationships. In particular, identifying E-G relationhsips seems crucial for a better understanding of complex genetic disorders such as Parkinson’s [[1]](#1), and might be of great importance towards a better understanding of the severity of haemochromatosis.

Hence, we proposed to used E-G relationships obtained from [[2]](#2) in order to identify genes sharing enhancers with those known to cause haemochromatosis or to regulate iron metabolism in liver and intestine, aiming at identifying new genes potentially involved in the severity of haemochromatosis.

## Objectives

We aim at finding new genes potentially related to haemochromatosis / involved in the severity of haemochromatosis, hopefully different from the ones already known to be strongly involved in iron metabolism. To that purpose, we started from a list of up to 13 genes, 6 of which being the ones causally involved in haemochromatosis, and the other ones being involved in iron metabolism regulation in liver and intestine.

We worked with E-G pairs obtained with either ABC or CHiC data, and inferred new genes as follows:

1. Find all the enhancers E regulating the list of 13 initial genes (later one we will denote by "I" a gene from this list)
2. Find all genes G regulated by those very same enhancers E

In other words, if G is a gene inferred, it means it is regulated by an enhancer E that is also regulating a initial gene, so we have a relation of the form "I <- E -> G".

## Results

We finally obtained a list of 453 putative new genes potentially involved in haemochromatosis. The list is available in the main [GitHub repository of this project](https://github.com/hoellin/hoellin.github.io) in `data_and_results/haemochromatosis/merged_inferred_genes_v1.csv`.

Those 453 genes are ranked according to a custom confidence index:

- 1 gene at the maximum confidence index (namely 15) : *COQ9*
- 11 genes at a confidence index >= 11 :

*COQ9, CX3CL1, POLR3GL, POLR3D, NUDT17, KATNB1, DOK4*

- 17 genes at a confidence index  >= 10
- 33 genes at a confidence index >= 8
- 44 genes at a confidence index >= 7
- 88 genes at a confidence index >= 6
- 148 genes at a confidence index >= 5
- 245 genes at a confidence index >= 4
- 363 genes at a confidence index >= 3
- 453 genes at a confidence index >= 2


## References

<a id="1">[1]</a> 
Farh, KH., Marson, A., Zhu, J. et al. (2015).
Genetic and epigenetic fine mapping of causal autoimmune disease variants. 
Nature 518, 337–343.

<a id="2">[2]</a>
Nasser, S., Farh, K.H., Zhu, J. et al. (2021).
Genome-wide mapping of enhancer-gene relationships in human cells.
Nature 591, 1–7.