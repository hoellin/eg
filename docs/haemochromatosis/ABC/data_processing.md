# Explore and prepare data for further analysis of E-G network with genes involved in hemochromatosis

## Requirements

### Data availability

All data are downloaded from [Nasser et al. 2021](https://www.engreitzlab.org/resources/) and available in `/work2/project/regenet/results/multi/abc.model/Nasser2021` on Genotoul.

`AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt` contains [predictions made with the ABC model over 131 biosamples from Nasser et al. 2021](https://www.engreitzlab.org/resources/).

> ```bash
> wget ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz -P /work2/project/regenet/results/multi/abc.model/Nasser2021/
> ```

List of biosamples linked to liver:

```bash
hepatocyte-ENCODE
HepG2-Roadmap # human liver cancer cell line
liver-ENCODE
```

List of biosamples linked to intestine:

```bash
large_intestine_fetal-Roadmap
small_intestine_fetal-Roadmap
```

List of biosamples linked to colon:

```bash
#HCT116-ENCODE # colon
#HT29 # colon
#sigmoid_colon-ENCODE # sigmoid colon
#transverse_colon-ENCODE # transverse colon
```

### Output

At the end of the day, we want to compute `bedpe` files the columns of which should be as follows:

```R
"chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "ABC.score", "strand1", "strand2", "biosample", "gene", "distance", "tissue"(, ...)
```

where "1" refers to enhancers and "2" refers to genes (either its TSS or the most upstream base of its most probable promoter according to what we want to do with the bedpe).

### `bedtools`

One can load the appropriate bedtools module using:

```bash
conda activate base && module load bioinfo/bedtools-2.27.1	
```

## Main remarks

### On merges of enhancers

* When, for a single biosample, we replace the enhancers by the "merged" enhancers (whatever they are) and we merge the resulting E-G pairs in which the merged enhancers overlap, we will have to merge ABC predictions eventually. This raises the following question: how to combine the ABC scores? Well, strictly speaking, there are no easy rigorous way to combine two ABC scores such that the score obtained is an ABC score too, even if we take the overlap into account (see next paragraph) because it could happen that the denominator in the two ABC scores are different, and would be different also if one looked at the merged enhancers from the beginning.

  ***Here we simply take the mean of the two ABC scores.*** 

A slightly better approach could be to do as follows (do a drawing to re-understand when reading again): let E1 and E3 be the two enhancers, such that the intersection of them is F2 ; and denote F1 = E1-F2, F3 = E3-F2 and F = F1+F2+F3 = union of E1 and E3. Consider the following ratios, the sum of which is 1: r1 = F1/F, r2 = F2/F, r3 = F3/F. Then, the new ABC score could read as follows: ABC = r1 x ABC1 + r3 x ABC3 + r2 x (ABC1+ABC3)/2.


## Extracting data of interest

> ```bash
> head -n 1 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt
> ```
>
> ```bash
> 1 chr
> 2 start
> 3 end
> 4 name
> 5 class
> 6 activity_base
> 7 TargetGene
> 8 TargetGeneTSS
> 9 TargetGeneExpression
> 10 TargetGenePromoterActivityQuantile
> 11 TargetGeneIsExpressed
> 12 distance
> 13 isSelfPromoter
> 14 hic_contact
> 15 powerlaw_contact
> 16 powerlaw_contact_reference
> 17 hic_contact_pl_scaled
> 18 hic_pseudocount
> 19 hic_contact_pl_scaled_adj
> 20 ABC.Score.Numerator
> 21 ABC.Score
> 22 powerlaw.Score.Numerator
> 23 powerlaw.Score
> 24 CellType
> ```
>
> ```bash
> wc -l AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt
> ```
>
> 7,717,393

### Main filters

First we verified that all genes in `AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt` are expressed.

> ```bash
> tail -n+2 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt |awk -F "\t" '{_[$11]++} END{for(u in _){print u, _[u]}}'
> ```
>
> ```bash
> True	7717392
> ```

OK.

Now we verify that all ABC scores are above 0.015 as expected:

> ```bash
> tail -n+2 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt |awk -F "\t" '{if($21>=0.015){print $0}}' |wc -l
> ```
>
> 7,717,392

OK.

We see that in a lot (> 1 million) of those element-gene pairs, the element is actually the promoter of the considered gene:

> ```bash
> tail -n+2 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt |awk 'BEGIN{FS="\t"; OFS="\t"} {_[$13]++} END{for(u in _){print u, _[u]}}'
> ```
>
> ```bash
> True	1177021
> False	6540371
> ```

Let's keep only the rows corresponding to enhancers and store the results, with the fields of interest, in a bedpe file.

```bash
conda activate base && module load bioinfo/bedtools-2.27.1	
```

```bash
tail -n +2 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt |awk 'BEGIN{FS="\t"; OFS="\t"} {if($13=="False"){print $1, $2, $3, $1, $8, $8, $4"::"$7, $21, ".", ".", $24, $7, $12}}' |bedtools sort -faidx chr_sizes -i stdin > Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe
```

> ```bash
> wc -l Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.bedpe
> ```
>
> 6,540,371

OK, this is the expected length.

In the resulting `Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe` bedpe file, there are no more ***self***-promoters, and the candidate elements are distributed as follows (a promoter of a gene G1 can be an enhancer for gene G2):

> ```bash
> tail -n +2 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt |awk -F "\t" '{if($13=="False"){_[$5]++}} END{for(u in _){print u, _[u]}}'
> ```
>
> ```bash
> intergenic 2758413
> promoter 224350
> genic 3557608
> ```

Now we extract the list of all enhancers involved in the predictions of interest for all 131 biosamples:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe |uniq > list_all_enhancers.bed
```

Contains 2,463,310 enhancers.

And we define the list of merged enhancers among those 131 biosamples, so that there are no overlaps in the resulting list:

```bash
bedtools merge -i list_all_enhancers.bed > list_all_enhancers.merged.bed
```

Contains 269,254 enhancers.

Now let's replace all enhancers coordinates in the result, by the coordinates of the corresponding merged enhancers in `list_all_enhancers.merged.bed`.

WARNING: it is **absolutely necessary** to have **at least** 32 GB of memory to run the code below. If the memory is not sufficient, the arrays indexed by `$11":::"$7` are not going to be stored but no error message will be yielded, which could make it pretty hard to debug the code.

```bash
srun --mem=32G --pty bash
conda activate base && module load bioinfo/bedtools-2.27.1
```

Note that the "." we add at the very last column is for compatibility with scrips used in downstream analysis, as when we select particular biosamples, this fields will contain the list of tissues of interest.

```bash
bedtools intersect -wb -a Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe -b list_all_enhancers.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,l1,":"); print $14, $15, $16, $4, $5, $6, $14":"$15"-"$16"::"l1[4], $8, $9, $10, $11, $12, $13}' |awk 'BEGIN{FS="\t"; OFS="\t"} {lines[$11":::"$7]=$0; score[$11":::"$7]+=$8; nb[$11":::"$7]++} END{for(u in lines){split(u,parts,":::"); split(lines[u],l,"\t"); print l[1], l[2], l[3], l[4], l[5], l[6], parts[2], score[u]/nb[u], l[9], l[10], l[11], l[12], l[13], "."}}' |bedtools sort -faidx chr_sizes -i stdin > Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

After a pretty long time:

> ```bash
> wc -l Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.new.bedpe
> ```
>
> 6,086,345

> To ensure that the Awk array has not been "cropped out" for insufficient memory reasons, we executed  the same code with twice as much max memory - eg 64 GB - to see if there is a difference (as it simply gives nothing with a small max memory) and it gave the same result => OK
>

Note that by construction, `Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe` contains almost duplicated rows when identical E-G are found in different biosamples, the only difference being the ABC score, the biosample and the distance fields.

We will also need the corresponding a similar files where all unique E-G correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

For that purpose, we use the `pandas` library from Python, which is more efficient. We use the script `reduce_df.py`, the content of which is:

```python
import pandas as pd
import numpy as np
import argparse
import sys

# Create the parser
parser = argparse.ArgumentParser(description="<input file>, <output file>")

# Add the arguments
parser.add_argument("-i",
                       metavar="path",
                       type=str,
                       required=True,
                       help="path to input file")
parser.add_argument("-o",
                       metavar="path",
                       type=str,
                       required=True,
                       help="path to output file")

args = parser.parse_args()

input_path = args.i
output_path = args.o

colnames = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name",
            "ABC.score", "strand1", "strand2", "biosample", "gene", "original_distance", "tissue"]
print("Loading input...")
eg_pairs = pd.read_csv(input_path, sep='\t', header=None, dtype='str', engine='c', names = colnames)
print("Done.")
eg_pairs[["start1", "end1", "start2", "end2", "ABC.score", "original_distance"]] = eg_pairs[["start1", "end1", "start2", "end2", "ABC.score", "original_distance"]].apply(pd.to_numeric)
eg_pairs[["ABC.score", "original_distance"]].astype(float, copy=False)

# For testing purposes we keep only a few rows
#eg_pairs = eg_pairs.iloc[:10000,:]

column_identical_for_unique_EG_pair = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name",
                                       "strand1", "strand2", "gene"]
print("Reducing input (this may take a while)...")
eg_reduced = (
    eg_pairs.groupby(column_identical_for_unique_EG_pair)
    .agg({"ABC.score": list, "biosample": list, "original_distance": list, "tissue": list})
    .reset_index()
    )
print("Done.")
eg_reduced.rename({"ABC.score": "ABC.scores", "biosample": "biosamples", "original_distance": "original_distances", "tissue": "tissues"}, axis=1, inplace=True)

print("Adding more columns...")
eg_reduced["ABC.mean"] = eg_reduced["ABC.scores"].apply(np.mean)
eg_reduced["ABC.min"] = eg_reduced["ABC.scores"].apply(np.min)
eg_reduced["ABC.max"] = eg_reduced["ABC.scores"].apply(np.max)

eg_reduced["original_distance.mean"] = eg_reduced["original_distances"].apply(np.mean)
eg_reduced["original_distance.min"] = eg_reduced["original_distances"].apply(np.min)
eg_reduced["original_distance.max"] = eg_reduced["original_distances"].apply(np.max)
print("Done.")

print("Preparing output...")
eg_reduced["ABC.scores"] = eg_reduced["ABC.scores"].apply(lambda x: ','.join(list(map(str, x))))
eg_reduced["original_distances"] = eg_reduced["original_distances"].apply(lambda x: ','.join(list(map(str, x))))
eg_reduced["biosamples"] = eg_reduced["biosamples"].apply(lambda x: ','.join(x))
eg_reduced["tissues"] = eg_reduced["tissues"].apply(lambda x: ','.join(x))
print("Done")

print("Writing output... Don't interrupt...")
eg_reduced.to_csv(output_path, sep='\t', header=True, index=False)
print("Done.")
sys.exit(0)
```

```bash
# requires a lot of memory, better use at least 32GB
python reduce_df.py -i Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.uniques_eg.bedpe
```

> ```bash
> wc -l Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe
> ```
>
> 6,086,345
>
> ```bash
> wc -l Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.uniques_eg.bedpe
> ```
>
> 1,041,155



### Select biosamples

#### All liver and intestine biosamples

The five biosamples of interest are:

> ```bash
> hepatocyte-ENCODE
> HepG2-Roadmap # human liver cancer cell line
> liver-ENCODE
> large_intestine_fetal-Roadmap
> small_intestine_fetal-Roadmap
> ```

We extract corresponding data:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {bool=0; if($11~/^(hepatocyte-ENCODE|HepG2-Roadmap|liver-ENCODE)$/){bool=1; tissue="liver"} else if ($11~/^(large_intestine_fetal-Roadmap|small_intestine_fetal-Roadmap)$/){bool=1; tissue="intestine"}; if(bool){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, tissue}}' Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

> ```bash
> awk -F "\t" '{tissues[$14]++} END{for(u in tissues){print u, tissues[u]}}' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe
> ```
>
> ```bash
> liver 145132
> intestine 101488
> ```

Note that by construction, `Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe` contains almost duplicated rows when identical E-G are found in different biosamples, the only difference being the ABC score, the biosample and the distance fields.

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.uniques_eg.bedpe
```

> ```bash
> wc -l Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe
> ```
>
> 246,620
>
> ```bash
> wc -l Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.uniques_eg.bedpe
> ```
>
> 147,748

We do the same with unmerged (original) enhancers:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {bool=0; if($11~/^(hepatocyte-ENCODE|HepG2-Roadmap|liver-ENCODE)$/){bool=1; tissue="liver"} else if ($11~/^(large_intestine_fetal-Roadmap|small_intestine_fetal-Roadmap)$/){bool=1; tissue="intestine"}; if(bool){print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, tissue}}' Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe
```

> ```bash
> awk -F "\t" '{tissues[$14]++} END{for(u in tissues){print u, tissues[u]}}' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe
> ```
>
> ```bash
> liver 157901
> intestine 108221
> ```

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.uniques_eg.bedpe
```

> ```bash
> wc -l Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe
> ```
>
> 266,122
>
> ```bash
> wc -l Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.uniques_eg.bedpe
> ```
>
> 265,975

And finally, we do the same again with enhancers merged only over the selected biosamples, ie the 5 liver and intestine biosamples:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.liver_and_intestine.bed
```

Contains 102,996 enhancers.

```bash
bedtools merge -i list_enhancers.liver_and_intestine.bed > list_enhancers.liver_and_intestine.merged.bed
```

Contains 57,398 enhancers.

So we replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.liver_and_intestine.merged.bed`.

```bash
bedtools intersect -wb -a Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe -b list_enhancers.liver_and_intestine.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,l1,":"); print $15, $16, $17, $4, $5, $6, $15":"$16"-"$17"::"l1[4], $8, $9, $10, $11, $12, $13, $14}' |awk 'BEGIN{FS="\t"; OFS="\t"} {lines[$11":::"$7]=$0; score[$11":::"$7]+=$8; nb[$11":::"$7]++} END{for(u in lines){split(u,parts,":::"); split(lines[u],l,"\t"); print l[1], l[2], l[3], l[4], l[5], l[6], parts[2], score[u]/nb[u], l[9], l[10], l[11], l[12], l[13], l[14]}}' |bedtools sort -faidx chr_sizes -i stdin > Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.bedpe
```

> ```bash
> awk -F "\t" '{tissues[$14]++} END{for(u in tissues){print u, tissues[u]}}' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.bedpe
> ```
>
> ```bash
> liver 153404
> intestine 105800
> ```

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.uniques_eg.bedpe
```

> ```bash
> wc -l Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.bedp
> ```
>
> 259,204
>
> ```bash
> wc -l Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.uniques_eg.bedpe
> ```
>
> 159,675

#### All liver biosamples

```bash
awk '$14=="liver"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.uniques_egbedpe
```



We do the same with unmerged (original) enhancers:

```bash
awk '$14=="liver"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.bedpe
```

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.uniques_eg.bedpe
```



And finally, we do the same again with enhancers merged only over the selected biosamples, ie the 3 liver biosamples:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.liver.bed
```

Contains 61,838 enhancers.

```bash
bedtools merge -i list_enhancers.liver.bed > list_enhancers.liver.merged.bed
```

Contains 44,979 enhancers.

So we replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.liver.merged.bed`.

```bash
bedtools intersect -wb -a Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.bedpe -b list_enhancers.liver.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,l1,":"); print $15, $16, $17, $4, $5, $6, $15":"$16"-"$17"::"l1[4], $8, $9, $10, $11, $12, $13, $14}' |awk 'BEGIN{FS="\t"; OFS="\t"} {lines[$11":::"$7]=$0; score[$11":::"$7]+=$8; nb[$11":::"$7]++} END{for(u in lines){split(u,parts,":::"); split(lines[u],l,"\t"); print l[1], l[2], l[3], l[4], l[5], l[6], parts[2], score[u]/nb[u], l[9], l[10], l[11], l[12], l[13], l[14]}}' |bedtools sort -faidx chr_sizes -i stdin > Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_liver_enhancers.sorted.bedpe
```

> ```bash
> awk -F "\t" '{tissues[$14]++} END{for(u in tissues){print u, tissues[u]}}' Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_liver_enhancers.sorted.bedpe
> ```
>
> ```bash
> liver 155293
> ```

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_liver_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_liver_enhancers.sorted.uniques_eg.bedpe
```

#### All intestine biosamples

```bash
awk '$14=="intestine"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_enhancers.sorted.uniques_eg.bedpe
```



```bash
awk '$14=="intestine"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.bedpe
```

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.uniques_eg.bedpe
```



And finally, we do the same again with enhancers merged only over the selected biosamples, ie the 2 intestine biosamples:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.intestine.bed
```

Contains 41,102 enhancers.

```bash
bedtools merge -i list_enhancers.intestine.bed > list_enhancers.intestine.merged.bed
```

Contains 27,102 enhancers.

So we replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.liver.merged.bed`.

```bash
bedtools intersect -wb -a Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.bedpe -b list_enhancers.intestine.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {split($7,l1,":"); print $15, $16, $17, $4, $5, $6, $15":"$16"-"$17"::"l1[4], $8, $9, $10, $11, $12, $13, $14}' |awk 'BEGIN{FS="\t"; OFS="\t"} {lines[$11":::"$7]=$0; score[$11":::"$7]+=$8; nb[$11":::"$7]++} END{for(u in lines){split(u,parts,":::"); split(lines[u],l,"\t"); print l[1], l[2], l[3], l[4], l[5], l[6], parts[2], score[u]/nb[u], l[9], l[10], l[11], l[12], l[13], l[14]}}' |bedtools sort -faidx chr_sizes -i stdin > Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_intestine_enhancers.sorted.bedpe
```

> ```bash
> awk -F "\t" '{tissues[$14]++} END{for(u in tissues){print u, tissues[u]}}' Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_intestine_enhancers.sorted.bedpe
> ```
>
> ```bash
> intestine 106883
> ```

Just as before, we use `reduce_df.py` to compute a bedpe in which all E-G pairs correspond to unique lines with $11 is the list of biosamples, $13 the list of initial E-TSS distances (before replacing the enhancers by merged enhancers) and $14 the list of tissues + other additional fields are defined, see header of the result.

```bash
python reduce_df.py -i Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_intestine_enhancers.sorted.bedpe -o Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_intestine_enhancers.sorted.uniques_eg.bedpe
```



#### All 5 ***separate*** biosamples of interest

##### hepatocyte-ENCODE


```bash
awk '$11=="hepatocyte-ENCODE"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.hepatocyte-ENCODE.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

```bash
awk '$11=="hepatocyte-ENCODE"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.hepatocyte-ENCODE.all_putative_enhancers.sorted.bedpe
```

And finally, we do the same again with enhancers merged only over the selected biosample:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.hepatocyte-ENCODE.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.hepatocyte-ENCODE.bed
```

Contains 23,436 enhancers.

```bash
bedtools merge -i list_enhancers.hepatocyte-ENCODE.bed > list_enhancers.hepatocyte-ENCODE.merged.bed
```

Contains 23,436 enhancers => indeed, it seems natural that enhancers in one biosample only do not overlap each other.

So we do not need to replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.hepatocyte-ENCODE.merged.bed`.

##### liver-ENCODE


```bash
awk '$11=="liver-ENCODE"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.liver-ENCODE.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

```bash
awk '$11=="liver-ENCODE"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.liver-ENCODE.all_putative_enhancers.sorted.bedpe
```

And finally, we do the same again with enhancers merged only over the selected biosample:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.liver-ENCODE.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.liver-ENCODE.bed
```

Contains 18,623 enhancers.

```bash
bedtools merge -i list_enhancers.liver-ENCODE.bed > list_enhancers.liver-ENCODE.merged.bed
```

Contains 18,623 enhancers => indeed, it seems natural that enhancers in one biosample only do not overlap each other.

So we do not need to replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.liver-ENCODE.merged.bed`.

##### HepG2-Roadmap


```bash
awk '$11=="HepG2-Roadmap"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.HepG2-Roadmap.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

```bash
awk '$11=="HepG2-Roadmap"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.HepG2-Roadmap.all_putative_enhancers.sorted.bedpe
```

And finally, we do the same again with enhancers merged only over the selected biosample:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.HepG2-Roadmap.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.HepG2-Roadmap.bed
```

Contains 19,766 enhancers.

```bash
bedtools merge -i list_enhancers.HepG2-Roadmap.bed > list_enhancers.HepG2-Roadmap.merged.bed
```

Contains 19,766 enhancers => indeed, it seems natural that enhancers in one biosample only do not overlap each other.

So we do not need to replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.HepG2-Roadmap.merged.bed`.

##### large_intestine_fetal-Roadmap


```bash
awk '$11=="large_intestine_fetal-Roadmap"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.large_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

```bash
awk '$11=="large_intestine_fetal-Roadmap"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.large_intestine_fetal-Roadmap.all_putative_enhancers.sorted.bedpe
```

And finally, we do the same again with enhancers merged only over the selected biosample:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.large_intestine_fetal-Roadmap.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.large_intestine_fetal-Roadmap.bed
```

Contains 20,517 enhancers.

```bash
bedtools merge -i list_enhancers.large_intestine_fetal-Roadmap.bed > list_enhancers.large_intestine_fetal-Roadmap.merged.bed
```

Contains 20,517 enhancers => indeed, it seems natural that enhancers in one biosample only do not overlap each other.

So we do not need to replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.large_intestine_fetal-Roadmap.merged.bed`.

##### small_intestine_fetal-Roadmap


```bash
awk '$11=="small_intestine_fetal-Roadmap"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe > Nasser2021ABCPredictions.small_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers.sorted.bedpe
```

OK.

```bash
awk '$11=="small_intestine_fetal-Roadmap"' Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe > Nasser2021ABCPredictions.small_intestine_fetal-Roadmap.all_putative_enhancers.sorted.bedpe
```

And finally, we do the same again with enhancers merged only over the selected biosample:

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3}' Nasser2021ABCPredictions.small_intestine_fetal-Roadmap.all_putative_enhancers.sorted.bedpe |uniq > list_enhancers.small_intestine_fetal-Roadmap.bed
```

Contains 20,555 enhancers.

```bash
bedtools merge -i list_enhancers.small_intestine_fetal-Roadmap.bed > list_enhancers.small_intestine_fetal-Roadmap.merged.bed
```

Contains 20,555 enhancers => indeed, it seems natural that enhancers in one biosample only do not overlap each other.

So we do not need to replace all enhancers coordinates in the result, by that of the corresponding merged enhancers in `list_enhancers.small_intestine_fetal-Roadmap.merged.bed`.


#### All biosamples

We already did it in section [Main filters](#main-filters).

* `Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe` for merged enhancers
* `Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.uniques_eg.bedpe` for merged enhancers (with 2 rows = 2 distinct unique E-G pairs)
* `Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe` for unmerged (original) enhancers

## Investigate putative regulatory elements

### Small R script to compute summary statistics

See `compute_summary_stats_enhancer_list.R`. Content:

> ```R
> options <- commandArgs(trailingOnly = TRUE)
> 
> enhancers = as.data.frame(read.table(options[1], sep = "\t"))
> enhancers$length = abs(enhancers$V3-enhancers$V2)
> summary(enhancers$length)
> ```
>

We need to load a module containing R to be able to launch the script:

```bash
module load system/R-4.1.1_gcc-9.3.0
```

We can execute the R script with:

```bash
Rscript compute_summary_stats_enhancer_list.R <relative_path_to_enhancers_list>
```

### ccRE

We take the ccRE from `/work2/project/regenet/results/multi/bengi/Annotations/hg19-cCREs.bed`.

```bash
cp /work2/project/regenet/results/multi/bengi/Annotations/hg19-cCREs.bed .
```

```bash
awk -F "\t" '$6~/(^Enhancer-like$)/' hg19-cCREs.bed > /work2/project/regenet/results/multi/abc.model/Nasser2021/ccRE-ELS.bed
```

989,712 ccRE-ELS

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   50.0   247.0   352.0   423.4   519.0 16633.0
```

### Merged enhancers of all 131 biosamples

#### Original enhancers

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   200.0   308.0   481.7   637.0  6991.0
```

#### Merged enhancers

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   288.0   631.0   807.8  1078.0 11616.0
```

#### Overlap with ccRE-ELS

> ```bash
> wc -l list_all_enhancers.merged.bed
> ```
>
> 269,254

```bash
bedtools intersect -c -a list_all_enhancers.merged.bed -b ccRE-ELS.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2	
```

> ```bash
> 1	112356
> 0	85937
> 2	44230
> 3	16616
> 4	6135
> 5	2319
> 6	934
> 7	388
> 8	164
> 9	78
> 10	43
> 11	24
> 12	14
> 13	7
> 14	4
> 17	3
> 16	2
> ```

85,937 (32%) merged putative enhancers across all 131 biosamples, do not match any ccRE-ELS.

> ```bash
> wc -l list_all_enhancers.bed
> ```
>
> 2,463,310

```bash
bedtools intersect -c -a list_all_enhancers.bed -b ccRE-ELS.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2	
```

> ```bash
> 0	1117560
> 1	1054002
> 2	221538
> 3	52355
> 4	12610
> 5	3471
> 6	1025
> 7	422
> 8	180
> 9	69
> 10	38
> 11	26
> 12	9
> 13	2
> 17	1
> 16	1
> 14	1
> ```

#### Overlap with all ccRE

```bash
bedtools intersect -c -a list_all_enhancers.merged.bed -b hg19-cCREs.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2
```

> ```bash
> 1	144690
> 2	58313
> 0	25709
> 3	23341
> 4	9498
> 5	4016
> 6	1854
> 7	899
> 8	422
> 9	235
> 10	114
> 11	63
> 12	41
> 13	24
> 14	11
> 15	6
> 17	5
> 16	5
> 20	3
> 18	3
> 22	1
> 19	1
> ```

Only 25,709 (9%) merged putative enhancers across all 131 biosamples, do not match any ccRE (vs 32% against ccRE-ELS), suggesting that a lot of those "putative enhancers" are actually other types of regulatory element as defined by Moore/ENCODE.

### Merged enhancers of all liver and intestine biosamples

#### Original enhancers

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   200.0   353.0   490.7   644.0  4626.0
```

#### Merged enhancers (merged across 131 biosamples)

```bash
bedtools intersect -wb -a list_enhancers.liver_and_intestine.bed -b list_all_enhancers.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {print $4, $5, $6}' |bedtools sort -faidx chr_sizes -i stdin |uniq > list_enhancers.liver_and_intestine.merged131.bed
```

> ```bash
> wc -l list_enhancers.liver_and_intestine.bed
> 102996
> wc -l list_enhancers.liver_and_intestine.merged131.bed
> 50881
> ```

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    200     520     982    1231    1668   11616
```

#### Merged enhancers (merged across 5 liver and intestine biosamples)

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   200.0   423.0   574.6   756.0  6298.0 
```

#### Overlap with ccRE-ELS

```bash
bedtools intersect -c -a list_enhancers.liver_and_intestine.merged.bed -b ccRE-ELS.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2	
```

> ```bash
> 0	26768
> 1	22443
> 2	6037
> 3	1564
> 4	398
> 5	135
> 6	33
> 7	10
> 9	6
> 10	2
> 8	1
> 16	1
> ```

47% do not overlap any ccRE-ELS !

### Merged enhancers of all liver biosamples

#### Original enhancers

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   200.0   296.0   455.3   585.0  4626.0
```

#### Merged enhancers (merged across 131 biosamples)

```bash
bedtools intersect -wb -a list_enhancers.liver.bed -b list_all_enhancers.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {print $4, $5, $6}' |bedtools sort -faidx chr_sizes -i stdin |uniq > list_enhancers.liver.merged131.bed
```

> ```bash
> wc -l list_enhancers.liver.bed
> 61838
> wc -l list_enhancers.liver.merged131.bed
> 39477
> ```

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    200     528    1027    1292    1780   11616
```

#### Merged enhancers (across all 3 liver biosamples)

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   200.0   348.0   507.1   650.0  4737.0
```

#### Overlap with ccRE-ELS

```bash
bedtools intersect -c -a list_enhancers.liver.merged.bed -b ccRE-ELS.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {nbs[$4]++} END{for(u in nbs){print u, nbs[u]}}' |sort -nrk2,2
```

> ```bash
> 0	22717
> 1	17340
> 2	3802
> 3	851
> 4	184
> 5	64
> 6	11
> 9	5
> 7	3
> 8	1
> 10	1
> ```

51% do not overlap any ccRE-ELS !

### Merged enhancers of all intestine biosamples

#### Original enhancers

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    200     200     428     543     724    3679
```

#### Merged enhancers (merged across 131 biosamples)

```bash
bedtools intersect -wb -a list_enhancers.intestine.bed -b list_all_enhancers.merged.bed |awk 'BEGIN{FS="\t"; OFS="\t"} {print $4, $5, $6}' |bedtools sort -faidx chr_sizes -i stdin |uniq > list_enhancers.intestine.merged131.bed
```

> ```bash
> wc -l list_enhancers.intestine.bed
> 41102
> wc -l list_enhancers.intestine.merged131.bed
> 24492
> ```

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    200     710    1209    1469    1978   11073
```

#### Merged enhancers (across all 2 intestine biosamples)

```R
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  200.0   221.0   488.0   606.8   814.0  4461.0
```


