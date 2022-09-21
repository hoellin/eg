# Context

We want to evaluate:

- the performance of Nasser et al 's ABC method to predict Enhancer-Gene relationships, over:
  - [the CRISPRi-FlowFISH set](../../notes_ABC/K562/recap_investigations_ABC_over_CRiFF/)
  - [the BENGI sets](../../notes_ABC/BENGI/ABC_BENGI_results/)
- the performance of Moore et al 's Average Rank method to predict Enhancer-Gene relationships, over:
  - [the CRISPRi-FlowFISH set](../CRISPRi_FlowFISH/avg_rank_method/summary_avg_rank/)
  - [the BENGI sets](../avg_rank_method/summary_avg_rank_method/)

In particular, we needed to understand:

- how to use the BENGI sets
- [how to use the Average Rank method](../avg_rank_method/avg_rank_method_with_code/)

This is the subject of this "BENGI" section.

# Introduction

These pages gather guidelines to use Moore et al. BENGI benchmarks, mainly for GM12878.

For most of them, the code can be directly adapted to re-run the steps anywhere needed. For some others, one should use the dedicated notebook. All notebooks can be found at [http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/](http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/).

