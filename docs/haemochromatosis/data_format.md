# Compute a list of reliable putative new genes

## Context

We aim at finding new genes potentially related to haemochromatosis / involved in the severity of haemochromatosis, hopefully different  from the ones already known to be strongly involved in iron metabolism. To that purpose, we started from a list of up to 13 genes, 6 of which being the ones causally involved in haemochromatosis, and the other ones being involved in iron metabolism regulation in liver and intestine.

We worked with E-G pairs obtained with either ABC or CHiC data, and inferred new genes as follows:

1. Find all the enhancers E regulating the list of 13 initial genes (later one we will denote by "I" a gene from this list)
2. Find all genes G regulated by those very same enhancers E

In other words, if G is a gene inferred, it means it is regulated by an enhancer E that is also regulating a initial gene, so we have a relation of the form "I <- E -> G".

So we worked with ABC and CHiC data, and obtained two separate list of inferred genes (of respective lengths 444 and 42). For more details on how those list have been obtained, one can refer to:

- http://genoweb.toulouse.inra.fr/~thoellinger/2022/ABC_inferred_genes/preliminary_analysis_v9.html for ABC data
- http://genoweb.toulouse.inra.fr/~thoellinger/2022/CHiC_ingerred_genes/preliminary_analysis_chic_v2.html for CHiC data

Now, we want to merge these list together, and to attribute a confidence label to each gene, according to multiple factors. The detailed process is explained below.

## Results

At the end of the day, one can directly use the following file:

`/work2/project/regenet/workspace/thoellinger/shared/2022/merged_inferred_genes_v1.csv`

## Remarks

- The field 2 (gene ids) is completed using a table then "by hand" at the very end
- The field 4 (source gene symbol) was added by hand at the very end (because I did not planned to add it at first), but could have been handled since the very beginning (as I added it in the input files, namely `new_genes_abc_v9_more_info.list` and `new_genes_chic_v2_more_info.list`)

## Data format

12 columns:

1. **gene** Gene symbol

2. **gene.id** ENSEMBL id (if applicable, "." otherwise)

3. **score.adj** The sum of all numerical $5-$9 labels (total "score" of the inferred gene) except that $9 (the number of enhancers) is divided by 4. See details of the following fields to understand why this score is expected to be relevant.

4. **source.gene** Gene symbol of the "source" gene (the initial gene I that shares an enhancer E with the current gene G)

5. **score.tot** The sum of all numerical $5-$9 labels (total "score" of the inferred gene, other approach). See details of the following fields to understand why this score is expected to be relevant.

6. **abc.enhancer.intersect.label** Label between 0 and 2. Let us denote G the current gene, and I <- E -> G the way it has been obtained.

   - 2 (the best) if E is a putative enhancer of the ABC data, and intersects >= 1 ccRE-ELS
   - 1 if E is a putative enhancer of the ABC data, and intersects 0 ccRE-ELS but >=1 ccRE
   - 0 otherwise (if the enhancer does not intersect any ccRE **or** if the gene G is not inferred from ABC element-Gene pairs)

7. **ABC.product.label** Label between 0 and 3. The higher the better. 0 if not inferred from ABC data, 1 if small ABC scores, 2 medium ABC score, 3 high ABC score ; according to the product of the ABC scores of:

   - the initial gene - enhancer pair (I-E)
   - the enhancer - gene pair (E-G)

   See the dedicated R markdown for more details ("small", "medium" and "high" have been defined in an arbitrary manner).

8. **chic.enhancer.intersect.label** Label between 0 and 2. Let us denote G the current gene, and I <- E -> G the way it has been obtained.

   - 2 (the best) if E is a putative enhancer of the CHiC data, and intersects >= 1 ccRE-ELS
   - 1 if E is a putative enhancer of the CHiC data, and  intersects 0 ccRE-ELS but >=1 ccRE
   - 0 otherwise (if the enhancer does not intersect any ccRE **or** if the gene G is not inferred from CHiC element-Gene pairs)

   Remark: as a matter of fact, there are no "1" in the csv file, because, to define the list of genes inferred with CHiC data, we restricting ourselves to E-G pairs in which the enhancers intersects >= 1 ccRE-ELS. We did so because putative enhancers in CHiC data are pretty big, and it would be surprising not to intersect any ccRE-ELS for an actual enhancer that is so large.

9. **contact.product.label** Label between 0 and 3. The higher the better. 0 if not inferred from CHiC data, 1 if small contact probability, 2 medium contact probability, 3 high contact probability ; according to the product of the contact probabilities of:

   - the initial gene - enhancer pair (I-E)
   - the enhancer - gene pair (E-G)

   See the dedicated R markdown for more details ("small", "medium" and "high" have been defined in an arbitrary manner).

10. **count** Label by the total number of E-G pairs involving the current gene (irrespective of if they are obtained from ABC or CHiC data). Note that this is nothing but the number of enhancers in $12.

11. **datatype** Indicates whether the gene comes from ABC data, CHiC data, or both "ABC,CHiC" data.

12. **sources** The list of enhancer coordinates of the I-E-G tuple (initial gene, enhancer, inferred gene).

The file is sorted by value of column $3.

## Steps to obtain this file

### Compute list of ABC-inferred genes

In order to do so, we started with a dedicated R markdown. Details here: http://genoweb.toulouse.inra.fr/~thoellinger/2022/ABC_inferred_genes/preliminary_analysis_v9.html

It produces an output file looking as follows

```
<gene_symbol>	<list_of_enhancers_ids>	<label_product_ABC_scores>	<nb_of_enhancers>
```

Now, we want to add relevant information to this list of genes (whether the associated enhancers intersect ccRE, etc). The steps are given below.

First, we compute:

  - a table which indicates, for every ABC enhancers id in our ABC data, whether the enhancer intersects >= 1 ccRE-ELS

    In `/work2/project/regenet/results/multi/abc.model/Nasser2021/`, we should find the 3 following files:

    - `list_all_enhancers.merged.bed`
    - `ccRE-ELS.bed`
    - `hg19-cCREs.bed`

    We copy them in our current working directory, namely, `/work2/project/regenet/workspace/thoellinger/shared/2022/data`.

    ```bash
    conda activate base && module load bioinfo/bedtools-2.27.1
    ```

    > ```bash
    > bedtools intersect -c -a list_all_enhancers.merged.bed -b ccRE-ELS.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2	
    > ```

    ```bash
    bedtools intersect -c -a list_all_abc_enhancers.merged.bed -b ccRE-ELS.bed |awk '$4>=1 {print $1":"$2"-"$3}' > list_all_abc_enhancers_intersecting_ccRE_ELS.bed
    ```

  - a table which indicates, for every ABC enhancers id in our ABC data, whether the enhancer intersects >= 1 ccRE

```bash
bedtools intersect -c -a list_all_abc_enhancers.merged.bed -b hg19-cCREs.bed |awk '$4>=1 {print $1":"$2"-"$3}' > list_all_abc_enhancers_intersecting_ccRE.bed
```

  - starting from those 2 tables, a table which associates to each ABC enhancer id:
    - 2 if it intersects >= 1 ccRE-ELS
    - 1 if it intersects 0 ccRE-ELS but >= 1 ccRE
    - 0 otherwise

Well, let us compute this table only for enhancers actually found in our data:

```bash
awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){els[$1]++; next}; if(els[$1]){print $1, 2} else{print $1, $2}}' data/list_all_abc_enhancers_intersecting_ccRE_ELS.bed <(awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){ccre[$1]++; next}; if(ccre[$1]){print $1, 1} else{print $1, 0}}' data/list_all_abc_enhancers_intersecting_ccRE.bed networks_hemochromatosis/results/enhancers.new_genes.list) > networks_hemochromatosis/results/enhancers.new_genes.intersection_label.list
```

> ```bash
> cd /work2/project/regenet/workspace/thoellinger/shared/2022/networks_hemochromatosis/results
> ```
>
> ```bash
> awk -F "\t" '{nb[$2]++} END{for(u in nb){print u, nb[u]}}' enhancers.new_genes.intersection_label.list
> ```
>
> 0 51
> 
> 1 501
> 
> 2 666

Now, in awk, add a column ($4) to our output files, which contains the list of scores of all enhancers id in column $2. Then, replace this column with the max in the list.

> Note: let us remark that there are only 101 distinct enhancers in our data:
> 
> ```bash
>awk -F "\t" 'BEGIN{OFS="\t"} {split($2,liste,","); for(i in liste){e=liste[i]; enhancers[e]++}} END{for(u in enhancers){print u, enhancers[u]}}' <(tail -n+2 new_genes_abc_v9_more_info.list) |wc -l
> ```
>
> 101

```bash
awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, "enhancer.intersect", $3, $4}' <(head -n 1 new_genes_abc_v9_more_info.list) && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){label[$1]=$2; next}; split($2,liste,","); labs=""; for(i in liste){e=liste[i]; labs=labs""label[e]}; print $1, $2, labs, $3, $4}' enhancers.new_genes.intersection_label.list <(tail -n+2 new_genes_abc_v9_more_info.list) |head
```

> ```bash
> awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, "enhancer.intersect", $3, $4}' <(head -n 1 new_genes_abc_v9_more_info.list) && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){label[$1]=$2; next}; split($2,liste,","); labs=""; for(i in liste){e=liste[i]; labs=labs""label[e]}; print $1, $2, labs, $3, $4}' enhancers.new_genes.intersection_label.list <(tail -n+2 new_genes_abc_v9_more_info.list) |awk -F "\t" 'BEGIN{OFS="\t"} {split($3,liste,""); max=0; for(i in liste){e=liste[i]; if(e>max){max=e}}; print $1, $2, max, $4, $5}' |head
> ```

```bash
(awk -F "\t" 'BEGIN{OFS="\t"} {print $1, "gene.id", $2, "enhancer.intersect", $3, $4}' <(head -n 1 new_genes_abc_v9_more_info.list) && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){label[$1]=$2; next}; split($2,liste,","); labs=""; for(i in liste){e=liste[i]; labs=labs""label[e]}; print $1, $2, labs, $3, $4}' enhancers.new_genes.intersection_label.list <(tail -n+2 new_genes_abc_v9_more_info.list) |awk -F "\t" 'BEGIN{OFS="\t"} {split($3,liste,""); max=0; for(i in liste){e=liste[i]; if(e>max){max=e}}; print $1, ".", $2, max, $4, $5}') > abc_eg.details.list
```

- Now, I should have an output file looking as follows

  ```
  <gene symbol>	<list_of_enhancers_ids>	<max_validity_score_of_enhaner>	<label_product_ABC_scores>	<nb_of_enhancers>
  ```

> So at the end of the day, how many Genes inferred are associated to at least one enhancer overlapping a ccRE-ELS or a ccRE?
>
> ```bash
> tail -n+2 abc_eg.details.list |awk -F "\t" '{nb[$4]++} END{for(u in nb){print u, nb[u]}}'
> ```
>
> ```bash
> 1 138 # at least 1 enhancer overlapping a ccRE
> 2 306 # at least 1 enhancer overlapping a ccRE-ELS
> ```

It appears that all the 444 inferred genes were obtain through at least 1 enhancer overlapping a ccRE.

### Compute list of CHiC-inferred genes

In order to do so, we started with a dedicated R markdown. Details here: http://genoweb.toulouse.inra.fr/~thoellinger/2022/CHiC_ingerred_genes/preliminary_analysis_chic_v2.html

It produces an output file looking as follows:

```
<gene_symbol>	<gene_id>	<list_of_enhancers_ids>	<label_product_CHiC_scores>	<nb_of_enhancers>
```

Now, we want to add relevant information to this list of genes (whether the associated enhancers intersect ccRE, etc). The steps are given below.

> Well, it was not really necessary to check whether the associated enhancers intersect ccRE or ccRE-ELS, because they all do: we restricted ourselves to CHiC data for which they do, since the beginning.

- We add a column ($4) to our output files, which contain the list of scores of all enhancers id in column $2. Then, we replace this column with the max in the list.

  > ```bash
  > awk -F "\t" 'BEGIN{OFS="\t"} {split($3,liste,","); for(i in liste){e=liste[i]; enhancers[e]++}} END{for(u in enhancers){print u, enhancers[u]}}' <(tail -n+2 new_genes_chic_v2_more_info.list) |wc -l
  > ```
  >
  > There are only 20 distinct enhancers involved.
  >
  > They must overlap with ccRE as with already restricted ourselves to data for which they do, but still, let us verify that they overlap ccRE-ELS:
  >
  > ```bash
  > awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){enhancers[$1":"$2"-"$3]++; next}; if(enhancers[$1]){print $0}}' ../data/list_all_enhancers.overlapping_ccRE-ELS.bed <(awk -F "\t" 'BEGIN{OFS="\t"} {split($3,liste,","); for(i in liste){e=liste[i]; enhancers[e]++}} END{for(u in enhancers){print u, enhancers[u]}}' <(tail -n+2 new_genes_chic_v2_more_info.list)) |wc -l
  > ```
  >
  > 20
  >
  > They do.
  
  All enhancers involved overlap >= 1 ccRE-ELS, so we add a score 2 everywhere, in the appropriate column.

  ```bash
  (awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, "enhancer.intersect", $4, $5}' <(head -n 1 new_genes_chic_v2_more_info.list) && awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, 2, $4, $5}' <(tail -n+2 new_genes_chic_v2_more_info.list)) > chic_eg.details.list
  ```

- Now, I should have an output file looking as follows

  ```
  <gene symbol>	<gene_id>	<list_of_enhancers_ids>	<max_validity_score_of_enhaner = 2> <label_product_CHiC_scores>	<nb_of_enhancers>
  ```

### Compute a table giving, for each gene symbol, its corresponding ENSEMBL id, whenever known

TODO.

### Combine everything

First, transform the file `abc_eg.details.list`, ie the one looking like:

> ```
> <gene symbol>	.	<list_of_enhancers_ids>	<label_product_ABC_scores>	<max_validity_score_of_enhaner>	<nb_of_enhancers>
> ```

so that it looks like:

> ```
> <gene symbol>	.	<list_of_enhancers_ids>	.	<max_validity_score_of_ABC_enhaner>	<label_product_ABC_scores>	0	0	<nb_of_enhancers>	ABC
> ```

So, go to `/work2/project/regenet/workspace/thoellinger/shared/2022/networks_hemochromatosis/results/` and do:

```bash
(echo -e "gene\tgene.id\tsources\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount" && awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, 0, 0, $6, "ABC"}' <(tail -n+2 abc_eg.details.list)) > abc_part.list
```

Then, transform the file looking like:

> ```
> <gene symbol>	<list_of_enhancers_ids>	<label_product_CHiC_scores>	<max_validity_score_of_enhaner>	<nb_of_enhancers>
> ```

so that it looks like:

> ```
> <gene symbol>	.	<list_of_enhancers_ids>	.	0	0	<max_validity_score_of_ChiC_enhancer>	<label_product_CHiC_scores>	<nb_of_enhancers>	CHiC
> ```

So, go to `/work2/project/regenet/workspace/thoellinger/shared/2022/promoter_capture_hic/results` and do:

```bash
(echo -e "gene\tgene.id\tsources\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount" && awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $3, 0, 0, $4, $5, $6, "CHiC"}' <(tail -n+2 chic_eg.details.list)) > chic_part.list
```

Then, concatenate the two files obtained above.

```bash
cd /work2/project/regenet/workspace/thoellinger/shared/2022
```

```bash
(cat networks_hemochromatosis/results/abc_part.list && tail -n+2 promoter_capture_hic/results/chic_part.list) |wc -l
```

Now, we have to merge rows for which the gene symbol (column $1) is the same. There should not be more than 2 occurrences of the same gene symbol because we started with only 2 lists, of unique gene symbols. So here are the steps (the code is given afterwards):

In awk, parse the file a first time and maintain, for each new gene, a table containing all column values. Whenever the current row corresponds to a gene already encountered in the above rows, do the following (note that the file is not sorted at the moment, so there should be first all ABC-inferred genes, then all CHiC-inferred genes, so whenever there is a gene already encountered, it means that the current row corresponds to a CHiC prediction):

- let us denote G the gene, and `info[G]` the information saved for this gene when previously encountered (which necessarily correspond to a row of ABC inference). Let $n denote the n-th column of the current row - that is, the information corresponding to the CHiC inference. We have to do the following:

  ```awk
  split(info[G],parts,"\t") # get the fields corresponding to the previously encountered occurence of the current gene
  parts[2] = $2 # replace the gene id by the current one (this won't change anything at the end, but we do so because the gene id is already known in CHiC data whereas it is not in ABC data ; but at the end of the day we will add all missing ids so it really changes nothing)
  parts[3] = parts[3]","$3 # add the current enhancers to the list
  parts[6] = $6 # add the CHiC enhancer label (which will always be 2 here because all CHiC enhancers we kept do intersect at least one ccRE-ELS) 
  parts[7] = $7 # add the CHiC contact product label (see definition in "Data format")
  parts[8] = parts[8]+$8 # sum the number of enhancers (ABC + CHiC)
  parts[9] = "ABC,CHiC" # if the gene is already encountered it means that it comes from both ABC and CHiC data
  ```

- when all the file is parsed - that is, at the END in awk, print all rows in `info`

- parse the resulting file once again, and simply compute, for each row, the sum of columns $5-$9, and put it in $4

- sort the result (descending) according to values of $4

So here is the code

> ```bash
> (echo -e "gene\tgene.id\tsources\tscore.sum\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount\tdatatype" && (tail -n+2 networks_hemochromatosis/results/abc_part.list && tail -n+2 promoter_capture_hic/results/chic_part.list) |awk -F "\t" 'BEGIN{OFS="\t"} {G=$1; if(!info[G]){info[G]=$0} else {split(info[G],parts,"\t"); parts[2] = $2; parts[3] = parts[3]","$3; parts[6] = $6; parts[7] = $7; parts[8] = parts[8]+$8; parts[9] = "ABC,CHiC"; n = length(parts); new=parts[1]; for(i=2;i<=n;i++){new=new""parts[i]}; info[G] = new}} END{for(g in info){print info[g]}}') |wc -l
> ```
>
> 454
>
> Instead of 444+42=486 inferred genes, we have 453 distinct inferred genes, which means we have 33 redundancies (in other words: 78% of the genes inferred with CHiC data, are inferred using ABC data too). Which is good!

```bash
(echo -e "gene\tgene.id\tscore.adj\tscore.tot\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount\tdatatype\tsources" && ((tail -n+2 networks_hemochromatosis/results/abc_part.list && tail -n+2 promoter_capture_hic/results/chic_part.list) |awk -F "\t" 'BEGIN{OFS="\t"} {G=$1; if(!info[G]){info[G]=$0} else {split(info[G],parts,"\t"); parts[2]=$2; parts[3] = parts[3]","$3; parts[6] = $6; parts[7] = $7; parts[8] = parts[8]+$8; parts[9] = "ABC,CHiC"; n = length(parts); new=parts[1]; for(i=2;i<=n;i++){new=new"\t"parts[i]}; info[G] = new}} END{for(g in info){print info[g]}}' |awk -F "\t" 'BEGIN{OFS="\t"} {print $1, $2, $4+$5+$6+$7+int($8/4), $4+$5+$6+$7+$8, $4, $5, $6, $7, $8, $9, $3}' |sort -nrk3,3)) > merged_inferred_genes_v1.list
```

Now we already have a good looking file. Lastly, we will add, whenever possible, the ensembl id corresponding to the gene symbol, using the table computed above for that purpose.

Here is the distribution of the scores obtained:

> ```bash
> tail -n+2 merged_inferred_genes_v1.list |awk -F "\t" '{nb_genes[$3]++} END{for(u in nb_genes){print u, nb_genes[u]}}' |sort -rnk1,1
> ```
>
> 15 1
>
> 12 1
>
> 11 6
>
> 10 9
>
> 9 6
>
> 8 10
>
> 7 11
>
> 6 44
>
> 5 60
>
> 4 97
>
> 3 118
>
> 2 90

Ok. Note that as I thought it would be better to give less importance to the number of enhancers, in the total score, I can divide it by 4 for the total score provided in the column $3 `score.adj`, whereas I kept it untouched for the total score provided in the column $4 `score.tot`. This is the only difference between these 2 columns.

### Add gene ids whenever possible

We use the table `/work2/project/regenet/results/multi/abc.model/Nasser2021/gname_gnid_full.tsv`

> ```bash
> awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){id[$1]=$2; next} split($0,parts,"\t"); n = length(parts); if(id[$1]){parts[2] = id[$1]} new=parts[1]; for(i=2;i<=n;i++){new=new"\t"parts[i]} print new}' /work2/project/regenet/results/multi/abc.model/Nasser2021/gname_gnid_full.tsv merged_inferred_genes_v1.list |head
> ```

```bash
mv merged_inferred_genes_v1.list backup.merged_inferred_genes_v1.list
awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){id[$1]=$2; next} split($0,parts,"\t"); n = length(parts); if(id[$1]){parts[2] = id[$1]} new=parts[1]; for(i=2;i<=n;i++){new=new"\t"parts[i]} print new}' /work2/project/regenet/results/multi/abc.model/Nasser2021/gname_gnid_full.tsv backup.merged_inferred_genes_v1.list > merged_inferred_genes_v1.list
mv merged_inferred_genes_v1.list merged_inferred_genes_v1.csv
```

OK.

Now we modify the csv by hand, to add ensembl id where missing.

### Add source gene symbol

One could also keep this information all along. I had to add it at the end because at first I did not think it would be useful.

```bash
cd .../2022/
```

> ```bash
> (echo -e "gene\tgene.id\tscore.adj\tgene.source\tscore.tot\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount\tdatatype\tsources" && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$5; next} source="."; if(from[$1]){source=from[$1]}; print $1, $2, $3, source, $4, $5, $6, $7, $8, $9, $10, $11}' networks_hemochromatosis/results/new_genes_abc_v9_more_info.list <(tail -n+2 merged_inferred_genes_v1.csv)) |head
> ```

> ```bash
> (echo -e "gene\tgene.id\tscore.adj\tgene.source\tscore.tot\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount\tdatatype\tsources" && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$6; next} source="."; if(from[$1]){if($4=="." || from[$1]==$4){source=from[$1]} else{source=$4","from[$1]}} else {source=$4}; print $1, $2, $3, source, $5, $6, $7, $8, $9, $10, $11, $12}' promoter_capture_hic/results/new_genes_chic_v2_more_info.list <(awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$5; next} source="."; if(from[$1]){source=from[$1]}; print $1, $2, $3, source, $4, $5, $6, $7, $8, $9, $10, $11}' networks_hemochromatosis/results/new_genes_abc_v9_more_info.list <(tail -n+2 merged_inferred_genes_v1.csv))) |head
> ```

> ```bash
> (echo -e "gene\tgene.id\tscore.adj\tgene.source\tscore.tot\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount\tdatatype\tsources" && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$6; next} source="."; if(from[$1]){if($4=="." || from[$1]==$4){source=from[$1]} else{source=$4","from[$1]}} else {source=$4}; print $1, $2, $3, source, $5, $6, $7, $8, $9, $10, $11, $12}' promoter_capture_hic/results/new_genes_chic_v2_more_info.list <(awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$5; next} source="."; if(from[$1]){source=from[$1]}; print $1, $2, $3, source, $4, $5, $6, $7, $8, $9, $10, $11}' networks_hemochromatosis/results/new_genes_abc_v9_more_info.list <(tail -n+2 merged_inferred_genes_v1.csv))) |awk -F "\t" 'BEGIN{OFS="\t"} {source[$4]++} END{for(u in source){print u, source[u]}}' |sort -nrk2,2
> ```

```bash
mv merged_inferred_genes_v1.csv temp && (echo -e "gene\tgene.id\tscore.adj\tgene.source\tscore.tot\tabc.enhancer.intersect.label\tABC.product.label\tchic.enhancer.intersect.label\tcontact.product.label\tcount\tdatatype\tsources" && awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$6; next} source="."; if(from[$1]){if($4=="." || from[$1]==$4){source=from[$1]} else{source=$4","from[$1]}} else {source=$4}; print $1, $2, $3, source, $5, $6, $7, $8, $9, $10, $11, $12}' promoter_capture_hic/results/new_genes_chic_v2_more_info.list <(awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){from[$1]=$5; next} source="."; if(from[$1]){source=from[$1]}; print $1, $2, $3, source, $4, $5, $6, $7, $8, $9, $10, $11}' networks_hemochromatosis/results/new_genes_abc_v9_more_info.list <(tail -n+2 temp))) > merged_inferred_genes_v1.csv && rm temp
```

