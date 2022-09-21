# Investigating reasons why we don't find the same predictions over the CRISPRi-FlowFISH dataset as Fulco et al, using similar accessions

In this markdown we investigate some possible reasons why we don't find the same predictions as Fulco et al. (namely, we find `AUPR=0.56` instead of `AUPR=0.66`).

## Recap

First of all let's recap what we've done to obtain those different predictions.

### Accessions

#### DNase OK

* dnase_rep1 = `ENCFF001DOX`  see https://www.encodeproject.org/files/ENCFF001DOX/

The name given in Supplementary Table 4 is actually `wgEncodeUwDnaseK562AlnRep1.bam` but it is mentioned in the link given above that `hg19/wgEncodeUwDnase/wgEncodeUwDnaseK562AlnRep1.bam` is nothing but the original name given to the replicate which is now referred to as `ENCFF001DOX`.

This replicate was added the 2011-01-11, was obtained through the [ENCODE2 Stamatoyannopoulos DNase Pipeline](https://www.encodeproject.org/pipelines/ENCPL104AAN/), and is mapped over the hg19 assembly (= GRCh37).

* dnase_rep2 = `wgEncodeUwDnaseK562AlnRep2`

Although I found the new accession name corresponding to the old name of replicate 1, I did not manage to find the new one corresponding to replicate 2,  `wgEncodeUwDnaseK562AlnRep2.bam`. Fortunately the replicate 2 was available for download, under its old name, here: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwDnase/

#### h3k27ac ChIP-seq OK

* h3k27_rep1 = `ENCFF384ZZM` # exp `ENCSR000AKP` (same as Fulco et al)
* h3k27_rep2 = `ENCFF070PWH` # exp `ENCSR000AKP` (same as Fulco et al)

Here is all the information given in Supplementary Table 4 : 

> ENCFF384ZZM, ENCFF070PWH
>
> ENCFF070PWH used for reproducibility analysis only (Supplementary Fig. 4). **Not** used to compute ABC score

Note that it is clearly written that **replicate 2 was not used** to compute the ABC score!

<span style="color:red">**We must ensure that we did not use the replicate 2 to compute the ABC score!**</span>

#### RNA-seq

At the time I computed the predictions, I did not manage to find the replicates used by Fulco et al! Instead, I used:

* rnaseq_rep1 = `ENCFF172GIN` # exp ENCSR000CPH (different from rnaseq used by Fulco et al)
* rnaseq_rep2 = `ENCFF768TKT` # exp ENCSR000CPH (different from Fulco et al rnaseq used by Fulco et al)

Fulco et al. give the [following GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87257) in Supplementary Table 4: `GSE87257`. The title of this serie in GEO browser is "Systematic mapping of functional enhancer-promoter connections with CRISPR interference". It has been published on Nov 11, 2016.

> | Experiment type | Expression profiling by high throughput sequencing           |
> | --------------- | ------------------------------------------------------------ |
> | Summary         | Gene expression in mammals is regulated  by noncoding elements that can impact physiology and disease, yet the  functions and target genes of most noncoding elements remain unknown. We present a high-throughput approach that uses CRISPR interference  (CRISPRi) to discover regulatory elements and identify their target  genes. We assess >1 megabase (Mb) of sequence in the vicinity of 2  essential transcription factors, MYC and GATA1, and identify 9 distal  enhancers that control gene expression and cellular proliferation.  Quantitative features of chromatin state and chromosome conformation  distinguish the 7 enhancers that regulate MYC from other elements that  do not, suggesting a strategy for predicting enhancer-promoter  connectivity. This CRISPRi-based approach can be applied to dissect  transcriptional networks and interpret the contributions of noncoding  genetic variation to human disease. |
> |                 |                                                              |
> | Overall design  | We examined the effects of using CRISPRi to inhibit a putative enhancer of GATA1, eHDAC6. We performed RNA  sequencing on K562 cells expressing individual sgRNAs targeting the  transcription start site of GATA1 (2 sgRNAs), eHDAC6 (2 sgRNAs) and  non-targeting, negative controls (4 sgRNAs). We generated paired-end RNA sequencing libraries from 3 biological replicates for each sgRNA. |

**Remark (digression - this has nothing to do with RNA-seq):** one should definitely pay attention to the third sentence of the summary, which is also mentioned in the first pages of Fulco et al. paper. That is, only sequences **in the vicinity** of target genes were assessed, ie although very reliable, the validation dataset based on CRISPR interference inherently misses highly distal regulatory elements.

**Remarque importante : en fait en me repenchant sur le fonctionnement du modèle ABC / le Github / le papier de Fulco et al, je vois qu'il n'est jamais nécessaire d'utiliser des données RNA-seq, sauf pour définir la liste des gènes exprimés ! Il ne faut pas confondre les données RNA-seq utilisées dans l'étape 2 seulement pour déterminer la liste des gènes exprimés (donc on n'utilise pas vraiment l'information quantitative, seulement un seuil), avec les données RNA-seq qui ont été utilisées seulement pour la constitution du dataset de validation par interférence CRISPR, et qui sont celles données dans le tableau Supplémentaire 4 ! (et c'est bien pour ça qu'elles proviennent, d'après les informations relatives à la référence GEO donnée dans le Tableau Supplémentaire 4, d'une série d'expériences de Fulco et al mais datant de 2016).**

L'origine de la confusion est la suivante : dans l'étape 2 d'ABC (cf Github) il faut effectivement utiliser une table d'expression des gènes... On voit que dans l'exemple sur le `chr22`, la table d'expression, au vu de son nom, est obtenue à partir du tsv `ENCFF934YBO`, issu de l'expérience `ENCSR000AEM` (on ne trouve pas le lien directement dans le résumé de l'expérience [ENCSR000AEM](https://www.encodeproject.org/experiments/ENCSR000AEM/) ; il faut plustôt se fier à la page de [ENCFF934YBO](https://www.encodeproject.org/files/ENCFF934YBO/) qui indique que ce tsv est obtenu à partir de [ENCFF213TSN](https://www.encodeproject.org/files/ENCFF213TSN/) et de [ENCFF826ONU](https://www.encodeproject.org/files/ENCFF826ONU/), le premier étant effectivement obtenu à partir de deux réplicats fastq issus de l'expérience `ENCSR000AEM`, le deuxième étant un fichier de génome de référence).

Ainsi visiblement, cette table d'expression ne sert qu'à savoir quels gènes sont exprimés, et quels gènes ne le sont pas (sans s'intéresser à leurs niveaux d'expression exacts). 

Sur le README.md du Github, voir aussi la section "Gene Expression":

> ## Gene Expression in ABC
>
> The ABC model is designed to predict the effect of enhancers on  expressed genes. If a gene is not expressed in a given cell type (or  cell state) then we assume it does not have any activating enhancers  (enhancers for which inhibition of the enhancer would lead to decrease  in gene expression). **Thus we typically only report enhancer-gene  connections for expressed genes.**
>
> In the absence of expression data, DNase-seq and H3K27ac ChIP-seq at  the gene promoter can be used as a proxy for expression. We suggest only considering enhancer-gene connections for genes with sufficiently  active promoters (**for instance in the top 60% of gene promoters in the  cell type**)

Est-ce qu'il faut comprendre que les auteurs ont considéré comme exprimés, les gènes dont l'expression fournie par le fichier de RNA-seq était dans le top 60% ?

Dans tous les cas, ce qui me parait le plus raisonnable, au vu de la forme du modèle ABC (score d'une relation E-G calculé dépendant de tous les éléments à moins de 5 Mb du gène G considéré), il me semble nécessaire de faire dans un premier temps des prédictions tout génome, puis ensuite seulement d'intersecter avec le dataset de validation, pour obtenir la performance. Et en effet, dans leur papier, Fulco et al. écrivent à propos de la méthode utilisée pour évaluer le pouvoir prédictif de leur modèle sur d'autres datasets de référence que CRISPRi-FlowFISH :

> We generated **genome-wide predictions** of functional enhancer–gene connections in each of these five cell types and compared them to the functional data in the corresponding cell type.

J'imagine que la démarche est la même pour leur dataset CRISPRi-FlowFISH : prédictions sur tout le génôme, puis intersection ensuite.

<span style="color:red">**Conclusion regarding RNA-seq: we should take a look at how it is used directly in the code, then recompute the ABC predictions with the same accession as for the example over chr22 found in Github, and NOT the accession given in Supplementary Table 4, because, on the one hand, it was used to create the CRISPRi-FlowFISH dataset, not to compute the ABC predictions, and, on the other and, we did not manage to find these files anyways.**</span>



Finalement sur Github on trouve : `ENCFF934YBO`

#### Gene annotation OK

* *gene_annotation*: `/work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf`
* *gnid_gname*: `/work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gnid.gnname.tsv`

#### Blacklist OK

Provided by Fulco et al. https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz

#### Ubiquitously expressed genes OK

Provided by Fulco et al.

### General remarks on what we have done and what we should have done

* We used `ENCFF172GIN` and `ENCFF768TKT` for gene quantifications (tsv, RNA-seq) (that we chose by ourselves based on... nothing, these replicates **were not** used by Fulco et al) but we should have used `ENCFF934YBO` (which is not mentioned in the paper, but on the example over chr22 on the Github). Gene quantification is not really essential for the ABC model (it is solely used to find the list of expressed genes), so it is not shocking that Fulco et al did not give the accession in their paper. Note that the RNA-seq experiment mentioned in Supplementary Table 4 has nothing to do with step 2 of the ABC method, nor with the ABC model in general ; but is solely about the generation of the CRISPRi-FlowFISH validation dataset!
* *We kept chromosomes only (we removed scaffolds) after calling peaks with macs2. I think that's ok.*
* We concatenated the counts from the 2 DNase-seq replicates. We must verify that Fulco et al did the same. Anyways, I think that we shall test used only replicate 1, using only replicate 2, and using both. I think it makes little sense if these 3 ways of doing provide very different results, as it would mean that the ABC model is not robust. But in their paper, Fulco et al wrote that they realized many tests and that all of them gave very similar results.
* At the end of the day, we obtained several distinct ABC scores for some enhancer-gene connection. We should test all the possibles ways of selecting 1 (take the one that maximize the score? average the scores? take the one that minimize the score?)



## Before we run the ABC model again

First of all, let's take a look at the code to understand what exactly are RNA-seq data used for.

### Step 2

In step 2 we run `run.neighborhoods.py` with the following args:

```bash
%%writetemplate $slurm_2
#!/bin/sh
python {run_neighborhoods} \
--candidate_enhancer_regions {peaks}{dnase_rep1}_{dnase_rep2}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes {annotations_dir}gene_ids.bed \
--H3K27ac {h3k27_file_rep1},{h3k27_file_rep1} \
--DHS {dnase_file_rep1},{dnase_file_rep2} \
--expression_table {gene_expression_table}  \
--chrom_sizes {genome_file} \
--ubiquitously_expressed_genes {ubiquitously_expressed_genes} \
--cellType {cell_type} \
--outdir {neighborhoods}
```

Now let's take a look at the argument parser in the code:

```python
...
    parser.add_argument('--expression_table', default="", nargs='?', help="Comma delimited string of gene expression files")  
...
```

Help for `expression_table` is "Comma delimited string of gene expression files". That does not help... + the other arguments were neither of any help (even though there are plenty unused arguments).

Now apart from the argument parsing, the code is as follows:

```python
def processCellType(args):
    params = parse_params_file(args)

    os.makedirs(args.outdir, exist_ok=True)

    #Setup Genes
    genes, genes_for_class_assignment = load_genes(file = args.genes, 
                                                ue_file = args.ubiquitously_expressed_genes,
                                                chrom_sizes = args.chrom_sizes,
                                                outdir = args.outdir, 
                                                expression_table_list = params["expression_table"], 
                                                gene_id_names = args.gene_name_annotations, 
                                                primary_id = args.primary_gene_identifier,
                                                cellType = args.cellType,
                                                class_gene_file = args.genes_for_class_assignment)

    if not args.skip_gene_counts:
        annotate_genes_with_features(genes = genes, 
                                        genome_sizes = args.chrom_sizes, 
                                        use_fast_count = (not args.use_secondary_counting_method),
                                        default_accessibility_feature = params['default_accessibility_feature'],
                                        features = params['features'],
                                        outdir = args.outdir)

    #Setup Candidate Enhancers
    load_enhancers(genes=genes_for_class_assignment, 
                    genome_sizes=args.chrom_sizes, 
                    candidate_peaks=args.candidate_enhancer_regions, 
                    skip_rpkm_quantile=args.skip_rpkm_quantile, 
                    qnorm = args.qnorm,
                    tss_slop_for_class_assignment=args.tss_slop_for_class_assignment,
                    use_fast_count = (not args.use_secondary_counting_method),
                    default_accessibility_feature = params['default_accessibility_feature'],
                    features = params['features'],
                    cellType = args.cellType,
                    class_override_file = args.enhancer_class_override,
                    outdir = args.outdir)

    print('Neighborhoods Complete! \n')

def main(args):
    processCellType(args)

if __name__ == '__main__':
    args = parseargs()
    main(args)
```

`expression_table` is only used in `load_genes`, which is imported from `neighborhoods.py`. So let's take a look at `neighborhoods.py`.

First, we recall that our expression table looks as follows:

> ```bash
> head CRISPRi_FlowFISH_K562/reference/expression/k562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt
> ```
>
> ```
> ENSG00000278842.1	0.38
> ENSG00000271876.1	0
> ENSG00000125975.13	0
> ENSG00000149313.10	42.285
> ENSG00000063515.2	0
> ENSG00000201567.1	0
> ENSG00000234601.1	0
> ENSG00000229745.1	0
> ENSG00000203950.6	22.135
> ENSG00000185522.8	1.435
> ```

also, looking at the argument parser of `run.neighborhoods.py ` and at the code of `neighborhoods.py`, we could have passed a list of expression tables instead of just one. So, even working with 2 files instead of 1 to compute gene expression, we could have created 2 separate files.

Of the following, the last line is important:

```python
name = os.path.basename(expression_table)
expr = pd.read_table(expression_table, names=[primary_id, name + '.Expression'])
expr[name + '.Expression'] = expr[name + '.Expression'].astype(float)
expr = expr.groupby(primary_id).max()

genes = genes.merge(expr, how="left", right_index=True, left_on='symbol')
```

It means that at this step, the expression of each gene is simply added to the list of genes passed as argument. This means that, starting from now, in `neighborhoods.py`, "genes" also contains gene expression. Let's determine what it is used for.

Remark: the following line,

```python
genes['Expression'] = genes[names_list].mean(axis = 1)
```

automatically averages the TPM expression values obtain for all expressions tables given as input, so it was really not necessary to average expression "by hand".

Now the following line is very important:

```python
genes['Expression.quantile'] = genes['Expression'].rank(method='average', na_option="top", ascending=True, pct=True)
```

That may be what we were searching for! A new column is created in the pandas dataframe `genes`, namely `Expression.quantile `, which contains the ranks of all genes according to their expression, in ascending order (meaning the first gene of the dataframe is the more expressed!)

Let's see how this column is used.

But first, let's note that:

* `expression_table_list` is never used again in the code
* The intermediates `epxr` are obviously not used again
* The column `'Expression'` of the pd dataframe `gene` is never used again
* Finally, only the `Expression.quantile ` column of `gene` keeps track of the expression table starting from now

well, and we reached the end of the `load_genes` function, meaning we can get back to `run.neighborhoods.py`. The next step is

```python
if not args.skip_gene_counts:
        annotate_genes_with_features(genes = genes, 
                                        genome_sizes = args.chrom_sizes, 
                                        use_fast_count = (not args.use_secondary_counting_method),
                                        default_accessibility_feature = params['default_accessibility_feature'],
                                        features = params['features'],
                                        outdir = args.outdir)
```

which takes us back to `run.neighborhoods.py`.

```python
def annotate_genes_with_features(genes, 
           genome_sizes,
           skip_gene_counts=False,
           features={},
           outdir=".",
           force=False,
           use_fast_count=True,
           default_accessibility_feature = "",
           **kwargs):

    #Setup files for counting
    bounds_bed = os.path.join(outdir, "GeneList.bed")
    tss1kb = make_tss_region_file(genes, outdir, genome_sizes)
    tss1kb_file = os.path.join(outdir, "GeneList.TSS1kb.bed")

    #Count features over genes and promoters
    genes = count_features_for_bed(genes, bounds_bed, genome_sizes, features, outdir, "Genes", force=force, use_fast_count=use_fast_count)
    tsscounts = count_features_for_bed(tss1kb, tss1kb_file, genome_sizes, features, outdir, "Genes.TSS1kb", force=force, use_fast_count=use_fast_count)
    tsscounts = tsscounts.drop(['chr','start','end','score','strand'], axis=1)

    merged = genes.merge(tsscounts, on="name", suffixes=['','.TSS1Kb'])

    access_col = default_accessibility_feature + ".RPKM.quantile.TSS1Kb"  

    if 'H3K27ac.RPKM.quantile.TSS1Kb' in merged.columns:
        merged['PromoterActivityQuantile'] = ((0.0001+merged['H3K27ac.RPKM.quantile.TSS1Kb'])*(0.0001+merged[access_col])).rank(method='average', na_option="top", ascending=True, pct=True)
    else:
        merged['PromoterActivityQuantile'] = ((0.0001+merged[access_col])).rank(method='average', na_option="top", ascending=True, pct=True)


    merged.to_csv(os.path.join(outdir, "GeneList.txt"),
             sep='\t', index=False, header=True, float_format="%.6f")

    return merged
```

Well nothing major is being done until the end of step 2.

The most important thing is that the `GeneList.txt` file generated during step 2, and used for step 3, now contains an `Expression.quantile ` column.

So let's take a look at step 3.

### Step 3

Let's investigate the content of `src/predict.py`. Well after the argument parser we see the following, which may answer our questions:

```python
def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()

    validate_args(args)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    write_params(args, os.path.join(args.outdir, "parameters.predict.txt"))
    
    print("reading genes")
    genes = pd.read_csv(args.genes, sep = "\t")
    genes = determine_expressed_genes(genes, args.expression_cutoff, args.promoter_activity_quantile_cutoff)
    genes = genes.loc[:,['chr','symbol','tss','Expression','PromoterActivityQuantile','isExpressed']]
    genes.columns = ['chr','TargetGene', 'TargetGeneTSS', 'TargetGeneExpression', 'TargetGenePromoterActivityQuantile','TargetGeneIsExpressed']
```

Let's investigate what the line

```python
genes = determine_expressed_genes(genes, args.expression_cutoff, args.promoter_activity_quantile_cutoff)
```

does, but it seems all clear! It explicitly cut genes list according to an expression cutoff!

> ```python
> parser.add_argument('--expression_cutoff', type=float, default=1, help="Make predictions for genes with expression higher than this value")
> ```

We see that by default, only genes with expression > 1 are considered, ie expressed genes. Given no expression cutoff is given as argument, this "1" threshold is indeed the one that is used.

`determine_expressed_genes` is from `tools.py`.

```python
def determine_expressed_genes(genes, expression_cutoff, activity_quantile_cutoff):
    #Evaluate whether a gene should be considered 'expressed' so that it runs through the model

    #A gene is runnable if:
    #It is expressed OR (there is no expression AND its promoter has high activity)

    genes['isExpressed'] = np.logical_or(genes.Expression >= expression_cutoff, np.logical_and(np.isnan(genes.Expression), genes.PromoterActivityQuantile >= activity_quantile_cutoff))

    return(genes)
```

> ```python
> parser.add_argument('--promoter_activity_quantile_cutoff', type=float, default=.4, help="Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data")
> ```



Now it's clear. And genes which are not expressed won't considered in establishing candidate E-G pairs. We shall note that genes with no TPM counts but the promoter of which has high activity (is in the 60% percentile of most expressed promoters), are considered as expressed.

We shall also note that the `Expression.quantile ` column of the `GeneList.txt` file is not used here, so it is given just for further analysis!



Now, let's ru-run the ABC model one very last time, making sure we're doing the right think at each single step!



## Gene annotation

This time we should use same gene annotation as Fulco et al... Well they DO NOT use GENCODE genes. They use RefSeq genes! https://en.wikipedia.org/wiki/RefSeq. Take a look at the name of `RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed ` in the small example over chr22 in Github, and at the Supplementary Table 5a. Their list contain 60445 genes.

> **Gene and TSS Annotation**
> We downloaded the UCSC RefSeq track (refGene, version 2017-03-08). This track contains multiple isoforms per gene symbol. We selected one TSS for each gene in the genome. To make this selection, we used the TSS used by the largest number of coding isoforms.
>
> ...
>
> When making genome-wide predictions, we removed genes corresponding to small RNAs (gene symbol contains ‘MIR’ or ‘RNU’, or gene body length <300 bp), as well as very long RNAs (gene body >2 Mb), which appear to correspond to artifactual UCSC transcript alignments.



## New run

https://www.sciencedirect.com/science/article/pii/S0022283619302530

Now let's try!

We should find 60 445 unique genes

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {split($9,parts,"\""); if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/) && ($3 ~ /(^transcript$)/)){if(300<=$5-$4 && $5-$4<=2000000){if(!(parts[2] ~ /(^MIR)|(^RNU)/)){print $1, $4, $5, parts[2]";"parts[4], $7}}}}' hg19.refGene.gtf > filtered_hg19.refGene.gtf
> ```
>
> ```bash
> wc -l filtered_hg19.refGene.gtf
> ```
>
> 70951
>
> ```bash
> bedtools sort -faidx ../ABC-Enhancer-Gene-Prediction/CRISPRi_FlowFISH_K562/reference/chr_sizes -i filtered_hg19.refGene.gtf > sorted_filtered_hg19.refGene.gtf
> ```

En fait rien que le choix des transcrits / des TSS ça va induire des différences notables...

On va faire deux choses : essayer avec tous les transcrits, et aussi essayer avec directement la liste donnée en table supplémentaire par FUlco et al.



```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); if(!_[$1$2$3,parts[1]]++){print $1, $2, $3, parts[1], 1, $5, parts[2]}}' sorted_filtered_hg19.refGene.gtf |head -n 25
```



> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); if(!_[parts[1]]++){print $0}}' sorted_filtered_hg19.refGene.gtf |wc -l
> ```
>
> 25556







```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if(!occurences_genes[gene]++){chr[gene]=$1}; if($5=="+"){TSS=$2} else {TSS=$3}; transcripts[gene]=trans";"transcripts[gene]; id=gene";"trans";"TSS; occurences[id]++; val[id]=$0;} END{for(u in transcripts){print u, transcripts[u]}}' sorted_filtered_hg19.refGene.gtf |head
```



```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2} else {TSS=$3}; transcripts[gene][trans]++; val[gene][trans]=$0;} END{for(u in transcripts){for(v in transcripts[u]){print u, v, transcripts[u][v]}}}' sorted_filtered_hg19.refGene.gtf |head

```



```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2} else {TSS=$3}; transcripts[gene][TSS]++; val[gene][TSS]=$0;} END{for(u in transcripts){for(v in transcripts[u]){print u, v, transcripts[u][v]}}}' sorted_filtered_hg19.refGene.gtf |head
```





```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2; TES=$3} else {TSS=$3; TES=$2}; transcripts[gene][TSS]++; TESs[gene][TSS]=TES";"TESs[gene][TSS]; val[gene][TSS][TES]=$0;} END{for(u in transcripts){for(v in transcripts[u]){print u, v, transcripts[u][v]}}}' sorted_filtered_hg19.refGene.gtf |head
```



```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2; TES=$3} else {TSS=$3; TES=$2}; transcripts[gene][TSS]++; TESs[gene][TSS]=TES";"TESs[gene][TSS]; val[gene][TSS][TES]=$0;} END{for(u in transcripts){for(v in transcripts[u]){print u, v, transcripts[u][v], TESs[u][v]}}}' sorted_filtered_hg19.refGene.gtf |sort -k1,1 -k3,3nr |head
```





```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2; TES=$3} else {TSS=$3; TES=$2}; transcripts[gene][TSS]++; TESs[gene][TSS]=TES";"TESs[gene][TSS]; val[gene][TSS][TES]=$0;} END{for(u in transcripts){for(v in transcripts[u]){split(substr(TESs[u][v], 1, length(TESs[u][v])-1),list_TESs,";"); for(t in list_TESs){print u, v, transcripts[u][v], t, list_TESs[t]}}}}' sorted_filtered_hg19.refGene.gtf |sort -k1,1 -k3,3nr |head
```





```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2; TES=$3} else {TSS=$3; TES=$2}; transcripts[gene][TSS]++; TESs[gene][TSS]=TES";"TESs[gene][TSS]; val[gene][TSS][TES]=$0;} END{for(u in transcripts){for(v in transcripts[u]){split(substr(TESs[u][v], 1, length(TESs[u][v])-1),list_TESs,";"); for(t in list_TESs){d=list_TESs[t]-v; print u, v, transcripts[u][v], d<0?-1*d:d}}}}' sorted_filtered_hg19.refGene.gtf |sort -k1,1 -k3,3nr -k4,4nr |head -n 30
```





C'est presque bon ! Il reste à ajouter un print des infos dont on a besoin (qu'il faut avoir enregistrées dans un tableau, $0 ne fonctionnera pas puisqu'on looop seulement lorsqu'on est arrivé à la dernière ligne !), puis à refaire passer ça dans un autre awk pour prendre les max. C'est simple comme  c'est déjà trié : on parcourt les lignes, on regarde si le gène est le même que le précédent, si oui on passe à la ligne suivante, sinon on affiche la ligne. Tout simplement !







```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; trans=parts[2]; if($5=="+"){TSS=$2; TES=$3} else {TSS=$3; TES=$2}; transcripts[gene][TSS]++; TESs[gene][TSS]=TES";"TESs[gene][TSS]; val[gene][TSS][TES]=$0;} END{for(u in transcripts){for(v in transcripts[u]){split(substr(TESs[u][v], 1, length(TESs[u][v])-1),list_TESs,";"); for(t in list_TESs){TES=list_TESs[t]; d=TES-v; print u, v, transcripts[u][v], d<0?-1*d:d, val[u][v][TES]}}}}' sorted_filtered_hg19.refGene.gtf |sort -k1,1 -k3,3nr -k4,4nr |head -n 30
```





### Well à refaire:

#### Gene annotation

We should do something like

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){if($3=="CDS"){split($9,fields,"\""); coding_genes[fields[2]]++;}; next}; if($3=="transcript"){split($9, fields, "\""); gene=fields[2]; if(coding_genes[gene]){if(genes[gene]){split(genes[gene],line,"\t"); if($4<line[2]){line[2]=$4}; if($5>line[3]){line[3]=$5}} else {line[2]=$4; line[3]=$5;} genes[gene]=$1"\t"line[2]"\t"line[3]"\t"gene"\t0\t"$7;}}} END{for(gene in genes){print genes[gene]}}' hg19.refGene.gtf hg19.refGene.gtf |wc -l
> ```
>
> 19433
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){if($3=="CDS"){split($9,fields,"\""); coding_genes[fields[2]]++;}; next}; if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/) && ($3=="transcript")){split($9, fields, "\""); gene=fields[2]; if(coding_genes[gene]){if(genes[gene]){split(genes[gene],line,"\t"); if($4<line[2]){line[2]=$4}; if($5>line[3]){line[3]=$5}} else {line[2]=$4; line[3]=$5;} genes[gene]=$1"\t"line[2]"\t"line[3]"\t"gene"\t0\t"$7;}}} END{for(gene in genes){print genes[gene]}}' hg19.refGene.gtf hg19.refGene.gtf |wc -l
> ```
>
> 19366

```bash
conda activate base && module load bioinfo/bedtools-2.27.1
srun --pty bash
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){if($3=="CDS"){split($9,fields,"\""); coding_genes[fields[2]]++;}; next}; if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/) && ($3=="transcript")){split($9, fields, "\""); gene=fields[2]; if(coding_genes[gene]){if(genes[gene]){split(genes[gene],line,"\t"); if($4<line[2]){line[2]=$4}; if($5>line[3]){line[3]=$5}} else {line[2]=$4; line[3]=$5;} genes[gene]=$1"\t"line[2]"\t"line[3]"\t"gene"\t0\t"$7;}}} END{for(gene in genes){print genes[gene]}}' hg19.refGene.gtf hg19.refGene.gtf |bedtools sort -faidx chr_sizes.bed -i > RefSeq_GRCh37_p13_coding_genes.sorted.bed
```



> ```bash
> wc -l RefSeq_GRCh37_p13_coding_genes.sorted.bed 
> ```
>
> 19366 RefSeq_GRCh37_p13_coding_genes.sorted.bed
>
> ```bash
> head RefSeq_GRCh37_p13_coding_genes.sorted.bed 
> ```
>
> ```
> chr1	69091	70008	OR4F5	0	+
> chr1	367659	180795226	OR4F3	0	+
> chr1	367659	180795226	OR4F16	0	-
> chr1	861111	879954	SAMD11	0	+
> chr1	879583	894636	NOC2L	0	-
> chr1	895964	901099	KLHL17	0	+
> chr1	901862	911245	PLEKHN1	0	+
> chr1	910578	917473	PERM1	0	-
> ```

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {l=$3-$2; if(l>300 && l<2000000){print $0}}' RefSeq_GRCh37_p13_coding_genes.sorted.bed |wc -l
> ```
>
> 19323

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {l=$3-$2; if(l>300 && l<2000000){print $0}}' RefSeq_GRCh37_p13_coding_genes.sorted.bed > RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed
```

#### OLD TSS annotation

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; $5=="+"?TSS=$2:TSS=$3; transcripts[gene][TSS]++; chr[gene][TSS]=$1; strand[gene][TSS]=$5;} END{for(u in transcripts){for(v in transcripts[u]){print chr[u][v], v, u, strand[u][v], transcripts[u][v]}}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed |sort -k3,3d -k5,5nr -u |awk 'BEGIN{FS="\t"; OFS="\t"} {if(!genes[$3]++){print $1, $2-250, $2+250, $3, 0, $4}}' |wc -l
> ```
>
> 19366

```bash
conda activate base && module load bioinfo/bedtools-2.27.1
srun --pty bash
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($4,parts,";"); gene=parts[1]; $5=="+"?TSS=$2:TSS=$3; transcripts[gene][TSS]++; chr[gene][TSS]=$1; strand[gene][TSS]=$5;} END{for(u in transcripts){for(v in transcripts[u]){print chr[u][v], v, u, strand[u][v], transcripts[u][v]}}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed |sort -k3,3d -k5,5nr -u |awk 'BEGIN{FS="\t"; OFS="\t"} {if(!genes[$3]++){print $1, $2-250, $2+250, $3, 0, $4}}' |bedtools sort -faidx chr_sizes.bed -i > RefSeq_GRCh37_p13_coding_genes_gt300bp_st2Mbp_max_nbtranscripts_per_TSS.TSS500bp.bed
```

### NEW TSS annotation

```bash
conda activate base && module load bioinfo/bedtools-2.27.1
srun --pty bash
```

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){curated[$4]++; next}; if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/) && ($3=="transcript")){split($9,parts,"\""); gene=parts[2]; if(curated[gene]){$7=="+"?TSS=$4:TSS=$5; transcripts[gene][TSS]++; chr[gene][TSS]=$1; strand[gene][TSS]=$7;}}} END{for(u in transcripts){for(v in transcripts[u]){print chr[u][v], v, u, strand[u][v], transcripts[u][v]}}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed hg19.refGene.gtf |sort -k3,3d -k5,5nr -u |awk 'BEGIN{FS="\t"; OFS="\t"} {if(!genes[$3]++){print $1, $2-250, $2+250, $3, 0, $4}}' |wc -l
> ```
>
> 19323

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){curated[$4]++; next}; if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/) && ($3=="transcript")){split($9,parts,"\""); gene=parts[2]; if(curated[gene]){$7=="+"?TSS=$4:TSS=$5; transcripts[gene][TSS]++; chr[gene][TSS]=$1; strand[gene][TSS]=$7;}}} END{for(u in transcripts){for(v in transcripts[u]){print chr[u][v], v, u, strand[u][v], transcripts[u][v]}}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed hg19.refGene.gtf |sort -k3,3d -k5,5nr -u |awk 'BEGIN{FS="\t"; OFS="\t"} {if(!genes[$3]++){print $1, $2-250, $2+250, $3, 0, $4}}' |bedtools sort -faidx chr_sizes.bed -i > RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_TSS_that_maximize_nb_of_distinct_coding_transcripts.TSS500bp.bed
```

> ```bash
> wc -l RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_TSS_that_maximize_nb_of_distinct_coding_transcripts.TSS500bp.bed
> ```
>
> 19323

It seems it worked!



## NOW:

Directly adding all candidate regulatory regions of the validation dataset in the whitelist, we obtained an AUPR of 0.63 and had no degeneracy in the predictions, yet we lacked predictions for a few hundreds validations pairs.

As we had no degeneracy in the predictions, hence no arbitrary choice to make, the predictions lacking are most likely due to the different RefSeq version of gene annotation we use. Let's check this.

There are 59 distinct genes in the CRISPRi-FlowFISH validation dataset:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {genes[$12]++} END{for(u in genes){print u}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
> 59



Yet only 54 of these 59 genes are found in the reference annotation we used:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]++; next}; if(genes[$12]){criff_genes[$12]++}} END{for(u in criff_genes){print u}}' /work2/project/regenet/workspace/thoellinger/RefSeq/RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
> 54



Even using Fulco et al 's annotation provided in their Supplementary Tables, only 58 out of these 59 genes are found:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){split($4,parts,";"); genes[parts[1]]++; next}; if(genes[$12]){criff_genes[$12]++}} END{for(u in criff_genes){print u}}' /work2/project/regenet/workspace/thoellinger/RefSeq/full_annotation_from_Fulco_et_al.gtf /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
>
> 58



> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){split($9,parts,"\""); genes[parts[2]]++; next}; if(genes[$12]){criff_genes[$12]++}} END{for(u in criff_genes){print u}}' /work2/project/regenet/workspace/thoellinger/RefSeq/hg19.refGene.gtf /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
>
> 56





### Intersection

The problem is that only gene names are given, so we need to find the loci of genes found in the CRiFF dataset, using the annotation given in the Supplementary Tables.

> ```
> chr1	11874	14409	DDX11L1;NR_046018	+
> chr1	14362	29370	WASH7P;NR_024540	-
> chr1	17369	17436	MIR6859-3;NR_107063	-
> chr1	17369	17436	MIR6859-2;NR_107062	-
> ```



```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/){split($4,fields,";"); gene=fields[1]; if(genes[gene]){split(genes[gene],line,"\t"); if($2<line[2]){line[2]=$2}; if($3>line[3]){line[3]=$3}} else {line[2]=$2; line[3]=$3}; genes[gene]=$1"\t"line[2]"\t"line[3]"\t"gene"\t0\t"$5;}} END{for(gene in genes){print genes[gene]}}' /work2/project/regenet/workspace/thoellinger/RefSeq/full_annotation_from_Fulco_et_al.gtf |bedtools sort -faidx chr_sizes.bed -i > gene_bodies_from_Fulco_et_al_annotation.bed
```

Now:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]=$0; next}; if(!unique[$12]++){if(genes[$12]){print genes[$12]}}}' /work2/project/regenet/workspace/thoellinger/RefSeq/gene_bodies_from_Fulco_et_al_annotation.bed /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
> 58

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]=$0; next}; if(!unique[$12]++){if(genes[$12]){print genes[$12]}}}' /work2/project/regenet/workspace/thoellinger/RefSeq/gene_bodies_from_Fulco_et_al_annotation.bed /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |bedtools sort -faidx chr_sizes.bed -i > /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/Fulco_annotation_for_genes_in_CRiFF.bed
```

Now:

```bash
bedtools intersect -sorted -wo -a /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/Fulco_annotation_for_genes_in_CRiFF.bed -b /work2/project/regenet/workspace/thoellinger/RefSeq/RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed -g chr_sizes.bed > /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/intersection_Fulco_and_newRefSeq_annotation_of_CRiFF_genes.bedpe
```

Finally:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"; k=0} {if(!OK[$4]){if($4==$10){OK[$4]=1}}; lines[k]=$0; k++} END{for(line in lines){split(lines[line],fields,"\t"); if(OK[fields[4]]){if(fields[4]==fields[10]){print lines[line]}} else {print lines[line]}}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/intersection_Fulco_and_newRefSeq_annotation_of_CRiFF_genes.bedpe |bedtools sort -faidx chr_sizes.bed -i |wc -l
> ```
>
> 56

This means that the 2 genes we are still missing with respect to the 58 ones found using directly Fulco et al 's annotation (before any filter), is likely because these 2 genes were considered as coding genes in the old RefSeq annotaiton, and as non-coding in the new RefSeq annotation we use.



So now, once and for all, I will recompute the ABC predictions 2 times with

* The 2 old gene names in the CRISPRi-FlowFISH validation dataset replaced by the 2 new gene names. So this one is easy: we don't have to re-run the ABC model, but only to reintersect the predictions
* The old RefSeq annotation from Fulco et al. For this one, we need to rerun the whole ABC model. But that will be easy as we do not have much things to change.
* Also, we should verify if we don't miss genes at each step of the ABC model because of the annotation.



```bash
awk 'BEGIN{FS="\t"; OFS="\t"; k=0} {if(!OK[$4]){if($4==$10){OK[$4]=1}}; lines[k]=$0; k++} END{for(line in lines){split(lines[line],fields,"\t"); if(OK[fields[4]]){if(fields[4]==fields[10]){print lines[line]}} else {print lines[line]}}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/intersection_Fulco_and_newRefSeq_annotation_of_CRiFF_genes.bedpe |bedtools sort -faidx chr_sizes.bed -i |awk 'BEGIN{FS="\t"; OFS="\t"} {print $4, $10}' > /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/old_RefSeq_names_new_RefSeq_name.txt
```



> ```bash
> wc -l /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames
> ```
>
> 3863
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){new[$1]=$2; next}; split($7,parts,"::"); if(new[$12]){print $1, $2, $3, $4, $5, $6, parts[1]"::"new[$12], $8, $9, $10, $11, new[$12]}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/old_RefSeq_names_new_RefSeq_name.txt /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
> 3854

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){new[$1]=$2; next}; split($7,parts,"::"); if(new[$12]){print $1, $2, $3, $4, $5, $6, parts[1]"::"new[$12], $8, $9, $10, $11, new[$12]}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/old_RefSeq_names_new_RefSeq_name.txt /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames > /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames.new
```



So now we re-intersect our ABC predictions:

```bash
cd /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/april_K562_candidates_in_whitelist
```

```bash
bedtools intersect -sorted -wo -a ABC_output/Predictions/AllPredictions.bedpe.sorted -b /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames.new -g reference/chr_sizes > ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.new.bedpe
```

keep only the lines for which the gene names match:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $0}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.new.bedpe |wc -l
> ```
>
> 3851

Now, in case there are multiple predictions for single validation pair, we arbitrarily keep the one that maximize the ABC score. We hope that such case does not happen, or not often ; let's see:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $23, $13, $14, $15, $24, $8, $25}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.new.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {uniq=$1"\t"$2"\t"$3"\t"$4"\t"$5; if($6>val[uniq]){val[uniq]=$6}} END{for(u in val){print u, val[u]}}' |wc -l
> ```

3851! So this case did not happen, and we lack only 12 predictions from the validation dataset, this time! Thats nice!

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $23, $13, $14, $15, $24, $8, $25}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.new.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {uniq=$1"\t"$2"\t"$3"\t"$4"\t"$5; if($6>val[uniq]){val[uniq]=$6}} END{for(u in val){print u, val[u]}}' > predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.new.bedpe
```

```bash
conda activate base && module load bioinfo/bedtools-2.27.1 && srun --mem=16G --pty bash
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $2, $3, $4, $5, $6, $1}' ABC_output/Predictions/predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.new.bedpe |bedtools sort -faidx reference/chr_sizes.bed -i |awk 'BEGIN{FS="\t"; OFS="\t"} {print $6, $1, $2, $3, $4, $5}' > ABC_output/Predictions/predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.new.bedpe.sorted
```

> ```bash
> wc -l ABC_output/Predictions/predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.new.bedpe.sorted
> ```
>
> 3851



Now we can perform further analyses with R.



Bon finalement ce n'était pas terrible. Il n'y a des prédictions que pour 99 des 109 positifs.



Ultimement, on va essayer de remplacer à la main les 3 gènes manquants.

## remplacement à la main des 3 derniers gènes manquants



> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$1]++; next}; if(!genes[$12]){criff_exclusive_gene[$12]++}} END{for(u in criff_exclusive_gene){print u, criff_exclusive_gene[u]}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/old_RefSeq_names_new_RefSeq_name.txt /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |head
> ```
>
> HBBP1	1 # Non codant !
> PVT1-TSS1	3 # Non codant !
> CCDC26	5 # Non codant !



Les 3 gènes en question sont non codant donc on les retirera définitivement du dataset de validation dans la suite ! Ici, on va les ajouter "à la main", même si c'est tout-à-fait artificiel et que ça introduit un certain biais, afin de vérifier qu'on obtient un AUPR proche de celui obtenu par Fulco et al. Encore une fois, on n'aura pas besoin de refaire tourner le modèle ABC, mais seulement de refaire l'intersection.



> ```
> (base) thoellinger@node120 /work2/project/regenet/workspace/thoellinger/RefSeq $ awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="transcript"){split($9,parts,"\""); if(parts[2] ~ /(^PVT1)/){print $0}}}' hg19.refGene.gtf |head
> ```
>
> chr8	refGene	transcript	128806779	129113499	.	+	.	gene_id "PVT1"; transcript_id "NR_003367";  gene_name "PVT1";
>
> ```
> (base) thoellinger@node120 /work2/project/regenet/workspace/thoellinger/RefSeq $ awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="transcript"){split($9,parts,"\""); if(parts[2] ~ /(^HBBP1)/){print $0}}}' hg19.refGene.gtf |head
> ```
>
> chr11	refGene	transcript	5263185	5264822	.	-	.	gene_id "HBBP1"; transcript_id "NR_001589";  gene_name "HBBP1";
>
> ```
> (base) thoellinger@node120 /work2/project/regenet/workspace/thoellinger/RefSeq $ awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="transcript"){split($9,parts,"\""); if(parts[2] ~ /(^CCDC26)/){print $0}}}' hg19.refGene.gtf |head
> ```
>
> chr8	refGene	transcript	130363940	130587264	.	-	.	gene_id "CCDC26"; transcript_id "NR_130919";  gene_name "CCDC26";
> chr8	refGene	transcript	130363940	130587264	.	-	.	gene_id "CCDC26"; transcript_id "NR_130918";  gene_name "CCDC26";
> chr8	refGene	transcript	130363940	130692485	.	-	.	gene_id "CCDC26"; transcript_id "NR_130917";  gene_name "CCDC26";
> chr8	refGene	transcript	130363940	130587264	.	-	.	gene_id "CCDC26"; transcript_id "NR_130920";  gene_name "CCDC26";



Ce qui nous donne (faisons ça vraiment à la main...) :

> chr8	128806779	129113499	PVT1	0	+
>
> chr11	5263185	5264822	HBBP1	0	-
>
> chr8	130363940	130692485	CCDC26	0	-



```bash
cp RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed temp.bed && echo -e "chr8\t128806779\t129113499\tPVT1\t0\t+
chr11\t5263185\t5264822\tHBBP1\t0\t-
chr8\t130363940\t130692485\tCCDC26\t0\t-" >> temp.bed && bedtools sort -faidx chr_sizes.bed -i temp.bed > RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed && rm temp.bed
```

Now:

```bash
cp /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/old_RefSeq_names_new_RefSeq_name.txt /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/all_old_RefSeq_names_new_RefSeq_name.txt
```

```bash
echo -e "PVT1-TSS1\tPVT1\nHBBP1\tHBBP1\nCCDC26\tCCDC26" >> /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/all_old_RefSeq_names_new_RefSeq_name.txt
```

> ```bash
> wc -l /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/all_old_RefSeq_names_new_RefSeq_name.txt 
> ```
>
> 59



> ```bash
> wc -l /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames
> ```
>
> 3863
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){new[$1]=$2; next}; split($7,parts,"::"); if(new[$12]){print $1, $2, $3, $4, $5, $6, parts[1]"::"new[$12], $8, $9, $10, $11, new[$12]}}' /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/all_old_RefSeq_names_new_RefSeq_name.txt /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames |wc -l
> ```
>
> 3863

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){new[$1]=$2; next}; split($7,parts,"::"); if(new[$12]){print $1, $2, $3, $4, $5, $6, parts[1]"::"new[$12], $8, $9, $10, $11, new[$12]}}' /work2/project/regenet/workspace /thoellinger/CRISPRi_FlowFISH/k562/all_old_RefSeq_names_new_RefSeq_name.txt /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames > /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames.byhand
```







So now we re-intersect our ABC predictions:

```bash
cd /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/april_K562_candidates_in_whitelist
```

```bash
bedtools intersect -sorted -wo -a ABC_output/Predictions/AllPredictions.bedpe.sorted -b /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames.byhand -g reference/chr_sizes > ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.byhand.bedpe
```

keep only the lines for which the gene names match:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $0}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.byhand.bedpe |wc -l
> ```
>
> 3851

Well unfortunately it changes nothing wrt the previous case! So we have to rerun all ABC predictions, in order to add 500bp TSS-centered regions for the 3 new genes in the whitelist.



## ABC 3rd RUN

* Include all the 59 (56 coding + 3 non-coding) genes in the whitelist!

```bash
conda activate base && module load bioinfo/bedtools-2.27.1
srun --pty bash
```

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){curated[$4]++; next}; if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/) && ($3=="transcript")){split($9,parts,"\""); gene=parts[2]; if(curated[gene]){$7=="+"?TSS=$4:TSS=$5; transcripts[gene][TSS]++; chr[gene][TSS]=$1; strand[gene][TSS]=$7;}}} END{for(u in transcripts){for(v in transcripts[u]){print chr[u][v], v, u, strand[u][v], transcripts[u][v]}}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed hg19.refGene.gtf |sort -k3,3d -k5,5nr -u |awk 'BEGIN{FS="\t"; OFS="\t"} {if(!genes[$3]++){print $1, $2-250, $2+250, $3, 0, $4}}' |wc -l
> ```
>
> 19326

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){curated[$4]++; next}; if(($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)/) && ($3=="transcript")){split($9,parts,"\""); gene=parts[2]; if(curated[gene]){$7=="+"?TSS=$4:TSS=$5; transcripts[gene][TSS]++; chr[gene][TSS]=$1; strand[gene][TSS]=$7;}}} END{for(u in transcripts){for(v in transcripts[u]){print chr[u][v], v, u, strand[u][v], transcripts[u][v]}}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed hg19.refGene.gtf |sort -k3,3d -k5,5nr -u |awk 'BEGIN{FS="\t"; OFS="\t"} {if(!genes[$3]++){print $1, $2-250, $2+250, $3, 0, $4}}' |bedtools sort -faidx chr_sizes.bed -i > RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_TSS_that_maximize_nb_of_distinct_coding_transcripts_with_3_more_genes_added_by_hand.TSS500bp.bed
```

> ```bash
> wc -l RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_TSS_that_maximize_nb_of_distinct_coding_transcripts_with_3_more_genes_added_by_hand.TSS500bp.bed
> ```
>
> 19326

It seems it worked!

### Problèmes à résoudre

* Lors de la construction du tableau d'expression des gènes, on a des gene ids, qu'on doit convertir en gene names. Verifier qu'on ne perd pas de gènes au passage. Eh si, on en perd. Beaucoup. On a seulement 791 des 847 gènes ubiquitaires qui se retrouvent dans notre annotation. Même en utilisant notre annotation complète sans filtre, on en a seulement 791. En utilisant l'annotation complète fournie en tableau supplémentaire par Fulco et al, on en aurait tout de même 831.

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]++; next}; if(genes[$1]){print $0}}' RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt |wc -l
> ```
>
>
> 791

Ce qu'on va faire est déjà récupérer ces 831, leurs coordonnées, et les réintersecter avec notre annotation, pour essayer d'en avoir un peu plus que 791...



> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]++; next}; if(genes[$1]){print $0}}' gene_bodies_from_Fulco_et_al_annotation.bed /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt |wc -l
> ```
>
> 831
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]++; next}; if(genes[$1]){print $0}}' curated_annotation_from_Fulco_et_al_table_5b.gtf /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt |wc -l
> ```
>
> 831
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]++; f13[$4]=$1"\t"$2"\t"$3; strand[$4]=$5; next}; if(genes[$1]){print f13[$1], $1, strand[$1]}}' curated_annotation_from_Fulco_et_al_table_5b.gtf /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt |head
> ```



```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$4]++; f13[$4]=$1"\t"$2"\t"$3; strand[$4]=$5; next}; if(genes[$1]){print f13[$1], $1, strand[$1]}}' curated_annotation_from_Fulco_et_al_table_5b.gtf /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19.txt > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b.bed
```

Now we can do the intersection.



```bash
cd RefSeq
conda activate base && module load bioinfo/bedtools-2.27.1
srun --pty bash
```

```bash
bedtools sort -faidx chr_sizes.bed -i /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b.bed > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b.sorted.bed && rm /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b.bed
```



```bash
bedtools intersect -sorted -wo -a /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b.sorted.bed -b RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed -g chr_sizes.bed > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.bedpe
```



> ```bash
> wc -l /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.bedpe
> ```
>
> 993

Now we keep only the appropriate lines:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"; k=0} {if(!OK[$4]){if($4==$9){OK[$4]=1}}; lines[k]=$0; k++} END{for(line in lines){split(lines[line],fields,"\t"); if(OK[fields[4]]){if(fields[4]==fields[9]){print lines[line]}} else {print lines[line]}}}' /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.bedpe |bedtools sort -faidx chr_sizes.bed -i |wc -l
> ```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"; k=0} {if(!OK[$4]){if($4==$9){OK[$4]=1}}; lines[k]=$0; k++} END{for(line in lines){split(lines[line],fields,"\t"); if(OK[fields[4]]){if(fields[4]==fields[9]){print lines[line]}} else {print lines[line]}}}' /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {print $6, $7, $8, $9, $10, $11}' |bedtools sort -faidx chr_sizes.bed -i > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.bed
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $4}' /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.bed > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.txt
```

OK, now we have 844 ubiquitous regions! with our custom RefSeq with 3 more genes. Sounds nice.





## ABC 4th RUN

* Include the 56 (not only 54 found previously) coding genes in the whitelist!



We also have to recreate the list of Ubiquitously expressed genes

```bash
cd RefSeq
conda activate base && module load bioinfo/bedtools-2.27.1
srun --pty bash
```

```bash
bedtools intersect -sorted -wo -a /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b.sorted.bed -b RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp.sorted.bed -g chr_sizes.bed > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.bedpe
```



> ```bash
> wc -l /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.bedpe
> ```
>
> 993

Now we keep only the appropriate lines:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"; k=0} {if(!OK[$4]){if($4==$9){OK[$4]=1}}; lines[k]=$0; k++} END{for(line in lines){split(lines[line],fields,"\t"); if(OK[fields[4]]){if(fields[4]==fields[9]){print lines[line]}} else {print lines[line]}}}' /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.bedpe |bedtools sort -faidx chr_sizes.bed -i |wc -l
> ```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"; k=0} {if(!OK[$4]){if($4==$9){OK[$4]=1}}; lines[k]=$0; k++} END{for(line in lines){split(lines[line],fields,"\t"); if(OK[fields[4]]){if(fields[4]==fields[9]){print lines[line]}} else {print lines[line]}}}' /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {print $6, $7, $8, $9, $10, $11}' |bedtools sort -faidx chr_sizes.bed -i > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.bed
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $4}' /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.bed > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation.txt
```

OK, now we have 844 blacklisted regions! with our custom RefSeq. Sounds nice.

