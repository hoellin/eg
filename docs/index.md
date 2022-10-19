# How to navigate through this site

These pages were first intented for my personal use to organize my notes, and to share them with my colleagues. Now they could be used as guidelines to help going through the data and scripts used to reproduce the results presented in [E/G relationships: need for more reliable genome-wide reference sets](https://doi.org/10.1101/2022.10.12.511908). Yet as they were not intented to be read by another audience, they are not very user-friendly. If you have any question or suggestion, please do not hesitate to [contact me](about_me).

## Context

If the billions of cells of a vertebrate organism all have the same genome, their functions are as diverse as the multitude of cell types these organisms are made of. These functional differences are due to the differential expression of genes in the cell types, which is due to the differential action of their regulatory elements (promoters, enhancers, insulators, etc).


Among those regulatory elements, enhancers appear as particularly interesting, not only because they are predominant and cover more genomic space [[1]](#1), but also because they appear to play important roles in human diseases [[2](#2),[3](#3)]. Enhancers, which, like promoters, are DNA elements bound by transcription factors (TF), are known to activate the expression of one or several genes by getting physically close to their promoters in the 3D space of the nucleus [[4](#4)-[5](#5)].


If the identification of all enhancers present in a given cell type is not a solved problem, the identification of the relationships between enhancers and genes, i.e. which genes are the targets of which enhancers, in a particular cell type is even more complex. 

Here we reviewed 2 of the most recent methods to identify enhancer-gene relationships.


## Organization

### This website

The GitHub Pages site is organized as follows:

- The [Guidebooks](guidebooks/introduction) section contains basic technical guidelines that are referenced in the other sections.
- The [ABC](notes_abc/introduction) section contains the scripts used to reproduce the results presented in the paper when using the ABC method.
- The [BENGI](notes_bengi/introduction) section contains the guidelines to:
    - use the Average Rank method
    - evaluate the ABC model and the Average Rank method over the BENGI sets
- In the [Haemochromatosis](haemochromatosis/introduction) section we make use of the ABC method to identify putative genes involved in the severity of haemochromatosis based on the enhancers they are regulated by.

### The GitHub repository

The present website is generated from the [GitHub repository](https://github.com/hoellin/eg) of the project. At the same time, the repository contains most (hopefully all - don't hesitate to [contact me](about_me) otherwise) of the resources required to reproduce the results presented in the paper. The repository itself is organized as follows:

- The `docs/` directory contains
    - the files used to generate the GitHub Pages site.
    - the scripts and notebooks used to reproduce the results presented in the paper (`docs/scripts/`))
    - some of the required data (`docs/data_and_results/`)
- The `.github/workflows/` directory contains the configuration files for the GitHub Actions workflows that are used to automatically generate the GitHub Pages site. It shall not be modified.
- The `mkdocs.yml` file contains the configuration and index of the GitHub Pages site. Any modification to the site structure should be done here (for instance if one adds a new file in `docs/`it should be referenced in `mkdocs.yml` to be accessible from the site - from then everything will be deployed automatically by GitHub Actions when pushing the commit).

## References

<a id="1">[1]</a> 
Pennacchio LA, Bickmore W, Dean A, Nobrega MA, Bejerano G.
Enhancers: five essential questions.
Nature Reviews Genetics. 2013 Apr;14(4):288-95

<a id="2">[2]</a>
Zhang G, Shi J, Zhu S, Lan Y, Xu L, Yuan H, Liao G, Liu X, Zhang Y, Xiao Y, Li X.
DiseaseEnhancer: a resource of human disease-associated enhancer catalog.
Nucleic acids research. 2018 Jan 4;46(D1):D78-84.

<a id="3">[3]</a>
Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR, Patwardhan TA, Jones TR, Nguyen TH, Ulirsch JC, Lekschas F, Mualim K. 
Genome-wide enhancer maps link risk variants to disease genes.
Nature. 2021 May;593(7858):238-43.

<a id="4">[4]</a>
Krivega I, Dean A.
Enhancer and promoter interactions—long distance calls.
Current opinion in genetics & development. 2012 Apr 1;22(2):79-85.

<a id="5">[5]</a>
Schoenfelder S, Fraser P.
Long-range enhancer–promoter contacts in gene expression control.
Nature Reviews Genetics. 2019 Aug;20(8):437-55.
