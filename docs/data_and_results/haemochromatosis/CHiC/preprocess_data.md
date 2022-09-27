# Explore and prepare data for further analysis of E-G network with genes involved in hemochromatosis, with E-G pairs based on CHiC data

## Requirements

### Data availability

All data can be found here:

`/work2/project/regenet/results/phic/homo_sapiens/hg19/jung.ren.2019/LI11/score2/pall/LI11.pall.score2.gninfo.bedpe`

```bash
cd /work2/project/regenet/results/phic/homo_sapiens/hg19/jung.ren.2019/LI11/score2/pall/
mkdir networks && cd networks
ln -s ../LI11.pall.score2.gninfo.bedpe raw_data_liver.full.bedpe
```

### `bedtools`

One can load the appropriate modules with:

```bash
conda activate base && module load bioinfo/bedtools-2.27.1	
```

### `Summary statistics with R`

We need to load a module containing R to be able to launch the `compute_summary_stats_enhancer_list.R` script:

```bash
module load system/R-4.1.1_gcc-9.3.0
```

```bash
cp /work2/project/regenet/results/multi/abc.model/Nasser2021/compute_summary_stats_enhancer_list.R .
```

We can execute the R script with:

```bash
Rscript compute_summary_stats_enhancer_list.R <relative_path_to_enhancers_list>
```

## Extracting data of interest

Note that the `LI11` repository already corresponds to putative E-G pairs in `liver` cells. The list of all available cell types can be found here; `/work2/project/regenet/results/phic/homo_sapiens/hg19/jung.ren.2019/metadata.tsv`.

There are 12 columns:

> ```bash
> head -n 1 raw_data_liver.full.bedpe
> ```
>
> ```bash
> chr1	943049	965801	chr1	1183838	1209657	chr1:943049:965801,chr1:1183838:1209657	2.76719343781528..	pp	11837	22511	7	242322	1.9015	3.6813	1.1317	2.1959	2	chr1:948802:949920:+:ENSG00000187608.5:ISG15,chr1:955502:991496:+:ENSG00000188157.9:AGRN,	1	chr1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,	elt2B	218037
> ```
>
> ```bash
> 1 chr1 # promoter of the gene
> 2 start1 # 0-based
> 3 end1 # 1-based
> 4 chr2 # elt2 (putative enhancer)
> 5 start2 # 0-based
> 6 end2 # 1-based
> 7 short name <chr1:start1:end1, chr2:start2:end2> # <prom_coord, elt2_coord>
> 8 strand1 # promoter
> 9 strand2 # elt2
> 10 score (does not matter because all p-values are small enough in those files)
> 11 type of relation (pp or po)  # I think that
> 								# pp stands for promoter-promoter
> 								# po stands for promoter-other
> 12-19 does not matter # cf mail "sens des champs des *gninfo.bedpe.gz" for details
> 20 number of genes whose tss+-500bp overlaps the prom parts
> 21 list of such genes in the form <chr1:start1:end1:strand1:gene_id:gene_symbol,chr2:...,...>
> 22 number of genes whose tss+-500bp overlaps the elt2 part
> 23 the list of such genes in the same form as for prom (21)
> 24 type of elt2 # two possible values:
> 				# elt2A if the element 2 is never a promoter anywhere else in the file
> 				# elt2B otherwise, that is if the element 2 also appears as a promoter somewhere else in the file
> 25 interfragment distance (from end to beg)
> ```

> ```bash
> awk -F "\t" '{types[$24]++} END{for(u in types){print u, types[u]}}' raw_data_liver.full.bedpe |sort -nrk2,2
> ```
>
> ```bash
> elt2A 33553
> elt2B 4706
> ```

### Exploration

#### Are there self-promoters ? If so, we have to remove them for comparison purpose with ABC data

> ```bash
> awk -F "\t" '$24~/(^elt2B$)/ {if($1":"$2":"$3==$4":"$5":"$6){print $0}}' raw_data_liver.bedpe |wc -l
> ```
>
> ```
> 0
> ```

=> No. OK.

### Main filters

#### Remove entries for which no gene TSS+-500 bp overlaps the promoter

This is the case for 46 entires:

> ```bash
> awk -F "\t" '{nb[$20]++} END{for(u in nb){print u, nb[u]}}' raw_data_liver.full.bedpe |sort -nrk2,2
> ```
>
> ```bash
> 1 30065
> 2 6672
> 3 1114
> 4 286
> 0 46
> 5 39
> 6 30
> 7 7
> ```

```bash
awk -F "\t" '$20>0' raw_data_liver.full.bedpe > temp.bedpe
mv temp.bedpe raw_data_liver.full.bedpe
```

> ```bash
> awk -F "\t" '{nb[$20]++} END{for(u in nb){print u, nb[u]}}' raw_data_liver.full.bedpe |sort -nrk2,2
> ```
>
> ```bash
> 1 30065
> 2 6672
> 3 1114
> 4 286
> 5 39
> 6 30
> 7 7
> ```

#### Deal with `pp` type of relation

Well there should be no problem keeping them:

> ```bash
> awk -F "\t" '{type[$11]++} END{for(u in type){print u, type[u]}}' raw_data_liver.full.bedpe |sort -nrk2,2
> ```
>
> ```bash
> po 33506
> pp 4707
> ```

#### Duplicates entries for which the E-Prom pair is associated to multiple E-G pairs, and modify file structure

```bash
awk -F "\t" 'BEGIN{OFS="\t"} {split($21,genes,","); i=0; n=length(genes)-1; while(i<n){i=i+1; gene=genes[i]; split(gene,parts,":"); entry[$4"\t"$5"\t"$6"\t"parts[1]"\t"parts[2]"\t"parts[3]"\t"parts[6]"\t"parts[5]]=$4"\t"$5"\t"$6"\t"parts[1]"\t"parts[2]"\t"parts[3]"\t"$4":"$5"-"$6"::"parts[5]"::"parts[6]"\t"$8"\t"$10"\t"$9"\tliver\t"parts[6]"\t"$25"\t"parts[5]}} END{for(u in entry){print entry[u]}}' raw_data_liver.full.bedpe > liver.all_putative_enhancers.bedpe
```

OK.

**Warning**:

- in `raw_data_liver.full.bedpe`, element 1 is the promoter and element 2 is the putative enhancer
- in `liver.all_putative_enhancers.bedpe`, element 1 is the putative enhancer and element 2 is the gene

Also, let us note that there are some duplicates, ie in a few cases (279), multiple entries in the original raw file, correspond to the same "id" (coordinates of the putative enhancer + coordinates, id and name of the gene) in the new file:

> ```bash
> awk -F "\t" 'BEGIN{OFS="\t"} {split($21,genes,","); i=0; n=length(genes)-1; while(i<n){i=i+1; gene=genes[i]; split(gene,parts,":"); entry[$4":"$5":"$6"::"parts[1]":"parts[2]":"parts[3]"::"parts[6]"::"parts[5]]++}} END{for(u in entry){if(entry[u]>1){print u, entry[u]}}}' raw_data_liver.full.bedpe |sort -nrk2,2 |wc -l
> ```
>
> 279

More precisely, there are 2 new entries that both match 3 old entries, and 277 new entries matching 2 old entries each.

The reason must be (I did not verify though, but this is the only explanation I see) that for some given putative enhancer, let us say for E, the are multiple promoters in the original file, let's say P1 and P2, each one overlapping the same TSS+- 500 bp of a gene G.

Well, for now we don't care.

```bash
cp /work2/project/regenet/results/multi/abc.model/Nasser2021/chr_sizes chr_sizes
```

```bash
bedtools sort -faidx chr_sizes -i liver.all_putative_enhancers.bedpe > liver.all_putative_enhancers.sorted.bedpe
```

#### Extract enhancers

Now we extract the list of all enhancers in our liver biosample:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' liver.all_putative_enhancers.sorted.bedpe |uniq > list_all_enhancers.bed
```

OK. Contains 31,749 enhancers.

And we define the list of merged enhancers, such that none are overlapping:

```bash
bedtools merge -i list_all_enhancers.bed > list_all_enhancers.merged_liver.bed
```

It did not change anything, so there was no overlap, which is expected as we are focused on only 1 biosample.

```bash
rm list_all_enhancers.merged_liver.bed
```

Here are the summary statistics on those enhancers:

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1111    4247    5596    7007    7756 1553733
```

**Note that they are pretty large!**

#### Overlap between enhancers and ccRE-ELS

We take the ccRE from `/work2/project/regenet/results/multi/bengi/Annotations/hg19-cCREs.bed`.

```bash
ln -s /work2/project/regenet/results/multi/bengi/Annotations/hg19-cCREs.bed .
```

```bash
awk -F "\t" '$6~/(^Enhancer-like$)/' hg19-cCREs.bed > ccRE-ELS.bed
```

989,712 ccRE-ELS

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   50.0   247.0   352.0   423.4   519.0 16633.0
```

```bash
bedtools intersect -c -a list_all_enhancers.bed -b ccRE-ELS.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2	
```

> ```bash
> 0	6149
> 1	4355
> 2	4238
> 3	3869
> 4	3302
> 5	2650
> 6	1975
> 7	1402
> 8	1040
> 9	744
> 10	568
> 11	363
> 12	247
> 13	176
> 14	148
> 15	105
> 16	90
> 17	73
> 18	48
> 19	34
> 20	31
> 22	19
> 23	18
> 21	18
> 25	13
> 24	13
> 26	9
> 34	7
> 30	7
> 28	7
> 27	5
> 32	4
> 29	4
> 33	3
> 49	2
> 39	2
> 37	2
> 36	2
> 31	2
> 51	1
> 42	1
> 41	1
> 40	1
> 35	1
> ```

6,149 (19%) of putative enhancers, do not match any ccRE-ELS. Only 4,355 (14%) of putative enhancers match exactly one ccRE-ELS. 

#### Overlap between enhancers and all ccRE

```bash
bedtools intersect -c -a list_all_enhancers.bed -b hg19-cCREs.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2
```

> ```bash
> 3	3962
> 2	3805
> 4	3629
> 1	3524
> 0	3184
> 5	3010
> 6	2443
> 7	1757
> 8	1397
> 9	998
> 10	785
> 11	589
> 12	435
> 13	351
> 14	292
> 15	237
> 16	198
> 17	162
> 19	131
> 18	120
> 20	93
> 21	76
> 22	75
> 24	63
> 23	56
> 26	36
> 28	34
> 25	34
> 29	32
> 27	25
> 30	24
> 33	22
> 31	22
> 32	16
> 36	13
> 35	13
> 38	12
> 37	12
> 34	12
> 41	11
> 40	6
> 39	5
> 53	4
> 52	4
> 51	4
> 45	4
> 42	4
> 50	3
> 49	3
> 44	3
> 43	3
> 48	2
> 46	2
> 79	1
> 76	1
> 74	1
> 69	1
> 64	1
> 61	1
> 60	1
> 58	1
> 57	1
> 55	1
> 54	1
> 47	1
> ```

3,184 (10%) of all putative enhancers, do not match any ccRE.

#### Intersect enhancers with Nasser2021 enhancers

TODO. Not essential.

#### Keep only putative enhancers that overlap at least one ccRE-ELS

At the end of the day we chose to cast out putative enhancers that do not overlap any ccRE-ELS, in order to remove element-Gene pairs in which the element is highly suspected not to be an enhancer. Note that we have not done so for Nasser2021 enhancers (instead, for Nasser2021 enhancers, we will create a confidence label on the resulting predictions according to whether the enhancer overlaps a ccRE-ELS or not), but for the ones here it is necessary given the important size of elements: if they do not overlap any ccRE-ELS despite their large size, they are very unlikely to be actual enhancers.

Doing so, we lost 19% of our enhancers, resulting in a list of 31,749−6,149 = 25,600 enhancers:

```bash
bedtools intersect -wa -a list_all_enhancers.bed -b ccRE-ELS.bed |uniq > list_all_enhancers.overlapping_ccRE-ELS.bed
```

> ```bash
> wc -l list_all_enhancers.overlapping_ccRE-ELS.bed
> ```
>
> 25,600

Then, over 48,038 initial putative E-G pairs in liver, there are 39,252 remaining pairs for which the enhancer overlaps at least one ccRE-ELS.

```bash
awk -F "\t" 'BEGIN{OFS="\t"} {if(NR==FNR){found[$1":"$2"-"$3]++; next}; if(found[$1":"$2"-"$3]){print $0}}' list_all_enhancers.overlapping_ccRE-ELS.bed liver.all_putative_enhancers.sorted.bedpe > liver.all_putative_enhancers.overlapping_ccRE-ELS.sorted.bedpe
```

> ```bash
> wc -l liver.all_putative_enhancers.overlapping_ccRE-ELS.sorted.bedpe
> ```
>
> 39,252

The file `liver.all_putative_enhancers.overlapping_ccRE-ELS.sorted.bedpe` is the one we will use from now.

