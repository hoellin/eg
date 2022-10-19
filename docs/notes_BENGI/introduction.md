# Context

We want to evaluate:

- the performance of ABC method [[1]](#1) to predict Enhancer-Gene relationships, over:
    - [the CRISPRi-FlowFISH set](../notes_ABC/K562/recap_investigations_ABC_over_CRiFF.html)
    - [the BENGI sets](../notes_ABC/BENGI/ABC_BENGI_results.html)
- the performance of Average Rank method [[2]](#2) to predict Enhancer-Gene relationships, over:
    - [the CRISPRi-FlowFISH set](CRISPRi_FlowFISH/avg_rank_method/summary_avg_rank.html)
    - [the BENGI sets](avg_rank_method/summary_avg_rank_method.html)

In particular, we needed to understand:

- how to use the BENGI sets
- [how to use the Average Rank method](avg_rank_method/avg_rank_method_with_code.html)

This is the subject of this section.

# Introduction

These pages gather guidelines to use Moore et al. BENGI benchmarks, mainly for GM12878.

For most of them, the code can be directly adapted to re-run the steps anywhere needed. For some others, one should use the dedicated notebook. All notebooks can be found at [scripts/notebooks/ipynb/](../scripts/notebooks/ipynb/).

# References

<a id="1">[1]</a> 
Moore, J.E., Pratt, H.E., Purcaro, M.J. et al.
A curated benchmark of enhancer-gene interactions for evaluating enhancer-target gene prediction methods.
Genome Biol 21, 17 (2020)

<a id="2">[2]</a>
Activity-by-Contact model of enhancer-promoter regulation from thousands of CRISPR perturbations.
Fulco CP, Nasser J, Jones TR, et al.
Nat Genet. 51(12):1664-1669 (2019)