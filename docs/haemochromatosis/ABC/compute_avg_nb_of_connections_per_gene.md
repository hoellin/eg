# Are the 10 genes of interest more connected to other genes, than genes are in general?

To answer this question, we need to compute:

* A table indicating for each gene of our annotation, the list of enhancers regulating it (easy with Awk)

* Another table indicating for each enhancer of our annotation, the list of genes regulated by it (easy with Awk)

Then, we compute, with R:

* The distribution of the number of connections made by our 10 genes, with R (well, with 10 genes only, better compute only mean, min and max)
* The distribution of the number of connections made by all genes

## All genes

### Compute tables

```bash
cd .../networks_hemochromatosis/data
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,parts,"::"); e=parts[1]; g=parts[2]; if(!new[e,g]++){if(reg[g]){reg[g]=reg[g]","e} else {reg[g]=e}}} END{for(g in reg){print g, reg[g]}}' <(head -n 1000 Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe) |head
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,parts,"::"); e=parts[1]; g=parts[2]; if(!new[e,g]++){if(reg[g]){reg[g]=reg[g]","e} else {reg[g]=e}}} END{for(g in reg){print g, reg[g]}}' Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe |sort -dk1,1 > g_to_regE.tsv
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,parts,"::"); e=parts[1]; g=parts[2]; if(!new[e,g]++){if(reg[e]){reg[e]=reg[e]","g} else {reg[e]=g}}} END{for(e in reg){print e, reg[e]}}' <(head -n 1000 Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe) |head
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,parts,"::"); e=parts[1]; g=parts[2]; if(!new[e,g]++){if(reg[e]){reg[e]=reg[e]","g} else {reg[e]=g}}} END{for(e in reg){print e, reg[e]}}' Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe > e_to_regG.tsv
```



By the way, we also computed a table giving, for each enhancer, the number of genes it is connected to.

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(!new[$7]++){split($7,parts,"::"); nb[parts[1]]++}} END{for(e in nb){print e, nb[e]}}' <(head -n 1000 Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe) |head
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if(!new[$7]++){split($7,parts,"::"); nb[parts[1]]++}} END{for(e in nb){print e, nb[e]}}' Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe > e_to_nbG.tsv
```

### Compute distributions

```bash
cd .../shared/automne_2021/networks_hemochromatosis/
```

We use the script `compute_avg_nb_of_connections_per_gene.R`, the content of which is:

```R

rm(list = ls())

options <- commandArgs(trailingOnly = TRUE)

table1_file = options[1]
table2_file = options[2]

print("Loading input...")

if (length(options)==3) {
  gene_source_file = options[3]
  gene.source = read.table(gene_source_file, sep = "\t")[,c(1)]
} else{
  gene.source = NULL
}

gE = as.data.frame(read.table(table1_file, sep = "\t"))
eG = as.data.frame(read.table(table2_file, sep = "\t"))
print("Done.")

colnames(gE) = c("gene", "enhancer.list")
colnames(eG) = c("enhancer", "gene.list")

# For development purpose only, reduce the dimensions:
#gE = gE[1:10000,]
#eG = eG[1:10000,]

if(!is.null(gene.source)){
  gE = gE[gE$gene %in% gene.source,]
}

print("Computing the list of targeted genes for each gene...")
gE$target.genes.list = lapply(gE$enhancer.list, function(liste){
  list_of_list_of_genes = eG[eG$enhancer %in% unlist(strsplit(liste,',')),]$gene.list
  list_of_genes = Reduce(function(x,y) {
    paste(unlist(sort(unique(c(strsplit(x,',')[[1]],strsplit(y,',')[[1]])))), collapse=',')
  }, list_of_list_of_genes)
  return(list_of_genes)
})
gE$target.genes.list = as.character(gE$target.genes.list)
print("Done.")

print("Counting the number of targeted genes for each gene...")
gE$target.genes.nb = lapply(gE$target.genes.list, function(liste){
  if(liste=="NULL"){
    res = 0}
  else{
    res = length(unlist(strsplit(liste, ',')))}
  return(res)
})
gE$target.genes.nb = as.numeric(gE$target.genes.nb)
print("Done...")

print("Summary statistics of number of target genes per gene:")
summary(gE$target.genes.nb)
```

Well I should have done this with Python and `pandas` rather than R! The execution is very slow:

```bash
Rscript compute_avg_nb_of_connections_per_gene.R data/g_to_regE.tsv data/e_to_regG.tsv 10genes.list
```

> ```R
>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
>    1.00   36.00   53.00   61.84   82.00  237.00
> ```

## 10 initial genes

```bash
Rscript compute_avg_nb_of_connections_per_gene.R data/g_to_regE.tsv data/e_to_regG.tsv 10genes.list
```

> ```R
>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
>    28.0    40.0    72.0    77.7   101.8   160.0
> ```

## Coding genes

### Extract full lists of coding / non-coding genes from RefSeq annotation

We go back on Genotoul.

```bash
cd /work2/project/regenet/results/multi/abc.model/Nasser2021/
```

> ```bash
> wc -l all_genes_in_eg_pairs_131_biosamples.list
> ```
>
> 23,220

Our RefSeq annotation is in `/work2/project/regenet/workspace/thoellinger/data/GRCh37/GRCh37.p13/genes.gtf` (private). We copy it here as `genes.refseq.gtf`

```bash
cp /work2/project/regenet/workspace/thoellinger/data/GRCh37/GRCh37.p13/genes.gtf genes.refseq.gtf
```

which contains 48,531 lines.

Here are the types of genes it contains:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if($9~/(gene_biotype)/){split($9,fields,/(gene_biotype)/); split(fields[2],subf,"\""); biotypes[subf[2]]++}} END{for(u in biotypes){print u, biotypes[u]}}' genes.refseq.gtf |sort -nrk2,2
```

> ```bash
> protein_coding	21663
> pseudogene	15537
> lncRNA	5706
> miRNA	2091
> transcribed_pseudogene	1268
> tRNA	660
> snoRNA	593
> V_segment	307
> V_segment_pseudogene	280
> J_segment	111
> snRNA	76
> D_segment	60
> C_region	32
> guide_RNA	31
> other	27
> rRNA	23
> antisense_RNA	19
> misc_RNA	13
> J_segment_pseudogene	11
> C_region_pseudogene	7
> Y_RNA	4
> vault_RNA	4
> scRNA	4
> telomerase_RNA	1
> RNase_P_RNA	1
> RNase_MRP_RNA	1
> ncRNA_pseudogene	1
> ```

So let's extract the list of protein coding genes only: 

```bash
awk -F "\t" '$9~/(gene_biotype)/ {split($9,fields,/(gene_biotype)/); split(fields[2],sub1,"\""); if(sub1[2]=="protein_coding"){if($9~/(gene_id)/){split($9,fields,/(gene_id)/); split(fields[2],sub2,"\""); print sub2[2]}}}' genes.refseq.gtf |sort -dk1,1 |uniq > genes.protein_coding.refseq.list 
```

Contains the 21,663 unique protein coding genes, OK.

Lets also extract the non-protein coding genes while we're here:

```bash
awk -F "\t" '$9~/(gene_biotype)/ {split($9,fields,/(gene_biotype)/); split(fields[2],sub1,"\""); if(sub1[2]!="protein_coding"){if($9~/(gene_id)/){split($9,fields,/(gene_id)/); split(fields[2],sub2,"\""); print sub2[2]}}}' genes.refseq.gtf |sort -dk1,1 |uniq > genes.not_protein_coding.refseq.list 
```

Contains the 26,868 unique non -protein-coding genes.

And finally let's simply extract the list of all genes of the annotation

```bash
awk -F "\t" '$9~/(gene_id)/ {split($9,fields,/(gene_id)/); split(fields[2],subf,"\""); print subf[2]}' genes.refseq.gtf |sort -dk1,1 |uniq > genes.refseq.list 
```

Contains the 48,531 unique genes in our RefSeq annotation.

### Keep Nasser2021 genes that are found in our RefSeq annotation

Now let's extract from all the 23,220 genes in the 131 biosamples, the ones contained in the annotation:

```bash
awk -F "\t" 'NR==FNR {genes[$1]++; next}; genes[$1]' genes.refseq.list all_genes_in_eg_pairs_131_biosamples.list > all_genes_in_eg_pairs_131_biosamples.in_annotation.list
```

Contains 21,845 genes.

### Extract lists of coding / non-coding genes both in our RefSeq annotation and in Nasser2021 list of genes

```bash
awk -F "\t" 'NR==FNR {genes[$1]++; next}; genes[$1]' genes.protein_coding.refseq.list all_genes_in_eg_pairs_131_biosamples.list > all_genes_in_eg_pairs_131_biosamples.in_annotation.protein_coding.list
```

Contains 17,665 genes.

```bash
awk -F "\t" 'NR==FNR {genes[$1]++; next}; genes[$1]' genes.not_protein_coding.refseq.list all_genes_in_eg_pairs_131_biosamples.list > all_genes_in_eg_pairs_131_biosamples.in_annotation.non_protein_coding.list
```



### Compute statistics on number of targeted genes

TODO

## Non coding genes

TODO



