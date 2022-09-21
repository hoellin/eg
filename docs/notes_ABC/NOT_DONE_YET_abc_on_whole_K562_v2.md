# ABC model on K562

Adapted from [ABC model's README on github](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

The Activity-by-contact model predicts which enhancers regulate which genes on a cell type specific basis.

## Characteristics for this run

* No whitelist whenever possible (although `genes` argument was required for step 2 : quantifying enhancer activity ; and the input in the original example for this argument was the whitelist. So we used the whitelist here).
* Always used the 2 replicates, but had to concatenate them by hand at the end of step 1 (`makeCandidateRegions.py`). 

## Requirements

For each cell-type, the inputs to the ABC model are:

 * Required Inputs

   * bam file for DNase-Seq or ATAC-Seq (indexed and sorted)

     The DNase-Seq file will contain the alignment of the genome-wide sequencing of regions sensitive to cleavage by DNase I, over the GRCh38 human's genome assembly.

   * bam file for H3K27ac ChIP-Seq (indexed and sorted)

     (Wiki) H3K27ac is an epigenetic modification to the DNA packaging protein Histone H3 (one of the five main histones involved in the structure of chromatin in eukaryotic cells). It is a mark that indicates the acetylation at the 27th lysine residue of the histone H3 protein. H3K27ac is associated with the higher activation of transcription and therefore defined as an *active* enhancer mark. H3K27ac is found at both proximal and distal regions of transcription start site (TSS).

 * Optional Inputs

   * Hi-C data

     These Hi-C data serve to compute the contact frequency between regions. If not available, they can be estimated as it has been proven that the contact frequency approximately follows a power law of the distance.

   * A measure of gene expression

In addition the following (non-cell-type specific) genome annotation files are required

 * bed file containing gene annotations (may change across cell types if using cell-type specific TSS's)
 * bed file containing chromosome annotations

## Data acquisition

### Chromatin accessibility (DNase-seq)

Experiment: `ENCSR000EOT` [K562 DNase-seq whole genome on GRCh38 genome assembly](https://www.encodeproject.org/search/?type=Experiment&status=released&perturbed=false&assay_title=DNase-seq&biosample_ontology.term_name=K562&assembly=GRCh38&files.file_type=bam&lab.title=John+Stamatoyannopoulos%2C+UW)

Replicate 1: `ENCFF156LGK` ; Replicate 2:  `ENCFF134DLD`

```bash
cd K562/
wget https://www.encodeproject.org/files/ENCFF156LGK/@@download/ENCFF156LGK.bam -P reference/dnase/
wget https://www.encodeproject.org/files/ENCFF156LGK/@@download/ENCFF134DLD.bam -P reference/dnase/
```

### Blacklist

We use the ENCODE Blacklist for Human GRCh38 downloaded [here](https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz).

### Histone mark H3K27ac ChIP-seq

Experiment `ENCSR000AKP` https://www.encodeproject.org/experiments/ENCSR000AKP/

Replicate 1: `ENCFF301TVL` ; Replicate 2: `ENCFF879BWC`

```bash
wget https://www.encodeproject.org/files/ENCFF301TVL/@@download/ENCFF301TVL.bam -P reference/H3K27ac/
wget https://www.encodeproject.org/files/ENCFF879BWC/@@download/ENCFF879BWC.bam -P reference/H3K27ac/
```

### Gene expression (polyA plus RNA-seq)

Experiment: `ENCSR000CPH` (Thomas Gingeras) https://www.encodeproject.org/experiments/ENCSR000CPH/

Replicate 1: `ENCFF172GIN` ; Replicate 2: `ENCFF768TKT`

```bash
wget https://www.encodeproject.org/files/ENCFF172GIN/@@download/ENCFF172GIN.tsv -P reference/expression/
wget https://www.encodeproject.org/files/ENCFF768TKT/@@download/ENCFF768TKT.tsv -P reference/expression/
```

### Genes annotations

We extract them from the header of `ENCFF156LGK.bam`.

```bash
cd K562/
samtools view -H reference/dnase/ENCFF156LGK.bam | grep SQ | cut -f 2,3 |cut -c 4- |awk 'BEGIN{FS="\t"} {split($2,locus,":"); if($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/){print($1"\t"locus[2])}}' > reference/chr_sizes

awk '{print $1"\t"0"\t"$2}' reference/chr_sizes > reference/chr_sizes.bed
```

## Data processing

At the end of the day we should have the following in `K562/reference/`:

> ```bash
> reference/
> ├── blacklist
> │   └── hg38-blacklist.v2.bed.gz
> ├── chr_sizes
> ├── chr_sizes.bed
> ├── dnase
> │   ├── ENCFF134DLD.bam # 30,981,578 lines (associated sam). 198 for header.
> │   └── ENCFF156LGK.bam # 355,014,590 lines (associated sam). 198 for header.
> ├── expression
> │   ├── ENCFF172GIN.tsv
> │   ├── ENCFF768TKT.tsv
> │   └── K562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt # mean TPM expression table
> ├── genes_names
> │   └── gencode.v35.annotation.gnid.gnname.tsv
> └── H3K27ac
>  ├── ENCFF301TVL.bam
>  └── ENCFF879BWC.bam
>    ```

### Create expression table

We first compared the resulting expression table using the 1st or the 2nd replicate. This resulted in significant differences in TPM values, but all 55,777 listed genes where the same.

```bash
cd K562/
awk 'BEGIN{FS="\t";} {if(NR==FNR){split($1,parts,"."); name[parts[1]] = $2; next;}; if(FNR>1 && $1 ~ /^ENS/){split($1,parts,"."); if(length(name[parts[1]]) != 0){if(!_[name[parts[1]]]++){print(name[parts[1]]"\t"$6)}}}}' reference/genes_names/gencode.v35.annotation.gnid.gnname.tsv reference/expression/ENCFF172GIN.tsv > reference/expression/v1_K562.ENCFF172GIN_ENCFF768TKT.TPM.txt
```

```bash
cd K562/
awk 'BEGIN{FS="\t";} {if(NR==FNR){split($1,parts,"."); name[parts[1]] = $2; next;}; if(FNR>1 && $1 ~ /^ENS/){split($1,parts,"."); if(length(name[parts[1]]) != 0){if(!_[name[parts[1]]]++){print(name[parts[1]]"\t"$6)}}}}' reference/genes_names/gencode.v35.annotation.gnid.gnname.tsv reference/expression/ENCFF768TKT.tsv > reference/expression/v2_K562.ENCFF172GIN_ENCFF768TKT.TPM.txt
```

Hence we decided to take as reference TPM value, the mean value for the 2 replicates.

<span style="color:red">**Is that OK?**</span>

To that purpose, write the following in `compute_mean_expression.awk`.

```bash
#!/bin/bash
BEGIN{
	FS="\t";
	token=0;
	count=0;
}
{
	if(NR==FNR){
		split($1,parts,".");
		name[parts[1]] = $2; # name[id] is the name associated to gene id
		next;
	}
	if(FNR == 1){token++; next;}
	if(token==1 && $1 ~ /^ENS/){
		split($1,parts,".");
		if(length(name[parts[1]]) != 0 && !_1[name[parts[1]]]++){
			# "!_1[name[parts[1]]]++" ensures uniqueness (checks whether it's the first time that name
			# is encountered in the current)
			tabTPM[name[parts[1]]] = $6;
			count++;
			name_found[count] = name[parts[1]];
		}
	}
	if(token == 2 && $1 ~ /^ENS/){
		split($1,parts,".");
		if(length(name[parts[1]]) != 0 && !_2[name[parts[1]]]++){
			# "!_2[name[parts[1]]]++" ensures uniqueness (checks whether it's the first time that name
			# is encountered in the current)
			tabTPM[name[parts[1]]] += $6;
			if(!_1[name[parts[1]]]){
				count++;
				name_found[count] = name[parts[1]];
				#print("well we should not be here");
			}
		}
	}
}
END{
	 for(i=1;i<=count;i++){
	 	print(name_found[i]"\t"tabTPM[name_found[i]]/2)
	 }
}
```

Then:

```bash
awk -f compute_mean_expression.awk reference/genes_names/gencode.v35.annotation.gnid.gnname.tsv reference/expression/ENCFF172GIN.tsv reference/expression/ENCFF768TKT.tsv > reference/expression/K562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt
```

## Running the ABC model

### Step 1: define candidate elements

> 'Candidate elements' are the set of putative enhancer elements for which ABC Scores will be computed. A typical way to define candidate elements is by calling peaks on a DNase-Seq or ATAC-Seq bam file. In this implementation we first call peaks using MACS2 and then process these peaks using `makeCandidateRegions.py`.

#### Call peaks with `macs2`

```bash
conda activate py2
```

```bash
macs2 callpeak \
-t reference/dnase/ENCFF156LGK.bam reference/dnase/ENCFF134DLD.bam \
-n ENCFF156LGK_ENCFF134DLD.macs2 \
-f BAM \
-g hs \
-p .1 \
--call-summits \
--outdir ABC_output/Peaks/
```

Seems to have worked. It took 23 minutes. Output `foo.macs2_peaks.narrowPeak` contains peaks of various lengths.

> ```bash
> ├── ABC_output
> │   ├── Neighborhoods
> │   └── Peaks
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_model.r # 17 lines
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak # 1,321,235 lines
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls # 1,321,263 lines
> │       └── ENCFF156LGK_ENCFF134DLD.macs2_summits.bed # 1,321,235 lines
> ```

> ```
> ...
> INFO  @ Wed, 02 Dec 2020 16:58:50: #2 Build Peak Model... 
> INFO  @ Wed, 02 Dec 2020 16:58:50: #2 looking for paired plus/minus strand peaks... 
> INFO  @ Wed, 02 Dec 2020 16:59:00: #2 number of paired peaks: 100801 
> INFO  @ Wed, 02 Dec 2020 16:59:00: start model_add_line... 
> INFO  @ Wed, 02 Dec 2020 16:59:01: start X-correlation... 
> INFO  @ Wed, 02 Dec 2020 16:59:01: end of X-cor 
> INFO  @ Wed, 02 Dec 2020 16:59:01: #2 finished! 
> INFO  @ Wed, 02 Dec 2020 16:59:01: #2 predicted fragment length is 60 bps 
> INFO  @ Wed, 02 Dec 2020 16:59:01: #2 alternative fragment length(s) may be 60 bps 
> INFO  @ Wed, 02 Dec 2020 16:59:01: #2.2 Generate R script for model : ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_model.r 
> WARNING @ Wed, 02 Dec 2020 16:59:01: #2 Since the d (60) calculated from paired-peaks are smaller than 2*tag length, it may be influenced by unknown sequencing problem! 
> WARNING @ Wed, 02 Dec 2020 16:59:01: #2 You may need to consider one of the other alternative d(s): 60 
> WARNING @ Wed, 02 Dec 2020 16:59:01: #2 You can restart the process with --nomodel --extsize XXX with your choice or an arbitrary number. Nontheless, MACS will continute computing. 
> ...
> ```

```bash
conda deactivate
conda activate base && conda activate abcmodel
```

#### Use ABC `makeCandidateRegions.py` to define candidate regions

> Now `makeCandidateRegions.py` will take as input the narrowPeak file produced by MACS2 and then perform the following processing steps:
> 
> 1. Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
> 2. Resize each of these N peaks to be a fixed number of base pairs centered on the peak summit
> 3. Remove any blacklisted regions and include any whitelisted regions
> 4. Merge any overlapping regions

**This time we do not use any whitelist.** Although initially there were ~ 25 000 white-listed regions. It seems that they basically correspond to curated gene annotations across the human genome.  Later we shall verify whether this has much influence or not on selected candidate regions.

As there are peaks over 112 scaffolds listed in `ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak`, we modify it (1,321,235 lines initially) and `ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls` (1,321,263 lines) to keep chromosomes only (otherwise `bedtools sort` below would not work). **But we have to keep in mind that `ENCFF156LGK_ENCFF134DLD.macs2_model.r` will still corresponds to DNase peaks over the 112 scaffolds**.

Maybe we better have kept whole `ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak` but filtered using another more complete "chrom_sizes" list, because later we see that `makeCandidateRegions.py` performs such a filter by itself.

```bash
awk '$1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/ {print $0}' ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak
```

> ```bash
> $ wc -l ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak 
> 1318093
> ```

```bash
awk 'NR <= 28 || $1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/ {print $0}' ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
```

> ```bash
> $ wc -l ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls 
> 1318121
> ```

Finally there are 1,321,235 - 1,318,093 = 3,142 DNase peaks on others scaffolds.

Now:

```bash
#First sort narrowPeak file
bedtools sort -faidx reference/chr_sizes -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
```

> ```bash
> $ wc -l ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
> 1318093
> ```

> ```bash
> K562
> ├── ABC_output
> │   ├── Neighborhoods
> │   └── Peaks
> │       ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak # 1318093
> │       ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_model.r
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak # 1321235
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted # 1318093
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls
> │       └── ENCFF156LGK_ENCFF134DLD.macs2_summits.bed
> ```

`../src/makeCandidateRegions.py` does not support multiple `bam` input, so we first need to concatenate the bam files associated with the 2 replicates.  Indeed, `makeCandidateRegions.py` will use the bam input to count DNase reads on each peak of `ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted` found with `macs2`.

```bash
samtools view -h -@ 3 reference/dnase/ENCFF156LGK.bam > reference/dnase/ENCFF156LGK_ENCFF134DLD.sam

samtools view reference/dnase/ENCFF134DLD.bam >> reference/dnase/ENCFF156LGK_ENCFF134DLD.sam
# now 385,995,970-lines long

samtools view -S -b reference/dnase/ENCFF156LGK_ENCFF134DLD.sam > reference/dnase/ENCFF156LGK_ENCFF134DLD.bam
```

Not enough space on local disk so we deported the computation of the intermediate `.sam` files on genotoul and copied the resulting compressed `.bam` onto local disk.

```bash
#Then
python ../src/makeCandidateRegions.py \
--narrowPeak ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted \
--bam reference/dnase/ENCFF156LGK_ENCFF134DLD.bam \
--outDir ABC_output/Peaks/ \
--chrom_sizes reference/chr_sizes \
--regions_blacklist reference/blacklist/hg38-blacklist.v2.bed.gz \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000
#Expected output: params.txt, ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed, ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK_ENCFF134DLD.bam.Counts.bed
```

> ```bash
> ├── ABC_output
> │   ├── Neighborhoods
> │   └── Peaks
> │       ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak # 1,318,093
> │       ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_model.r
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak # 1,321,235
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted # 1,318,093
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed # 124,770
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK_ENCFF134DLD.bam.Counts.bed # 1,317,935
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_summits.bed
> │       └── params.txt
> ```

Later on we may go back to the beginning of that step and try taking into account only the reads from `ENCFF156LGK` instead of `ENCFF156LGK_ENCFF134DLD`, to see if this has a significant influence on the final results.

By the way:

```bash
$ python ../src/makeCandidateRegions.py --help
```

> ```bash
> --regions_whitelist REGIONS_WHITELIST
>                         Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blacklist (default: )
> ```

### Step 2: quantifying enhancer activity

**We tried to let the `genes` argument empty but this led to the following error. Initially the whitelist was given as `gene` argument.**

> ```bash
> run.neighborhoods.py: error: the following arguments are required: --genes
> ```

`--help` yields:

>
> ```bash
> --genes GENES         bed file with gene annotations. Must be in bed-6 format. Will be used to assign TSS to genes. (default: None)
> ```

<span style="color:red">**So it seems that `../reference/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed` is nothing else but gene annotation over the human genome. So why is it called a `whitelist` in step 1 ?**</span>

```bash
python ../src/run.neighborhoods.py \
--candidate_enhancer_regions ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes ../reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed \
--H3K27ac reference/H3K27ac/ENCFF301TVL.bam,reference/H3K27ac/ENCFF879BWC.bam \
--DHS reference/dnase/ENCFF156LGK.bam,reference/dnase/ENCFF134DLD.bam \
--expression_table reference/expression/K562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt \
--chrom_sizes reference/chr_sizes \
--ubiquitously_expressed_genes ../reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir ABC_output/Neighborhoods/
```

Quite long, yields "BEDTools completed successfully." at each step. Then finally:

> ```bash
> Feature DHS completed in 575.4569108486176
> Assigning classes to enhancers
> Total enhancers: 124770
>          Promoters: 3338
>          Genic: 0
>          Intergenic: 121432
> Neighborhoods Complete!
> ```

All output files in `Neightbhorhoods` seems ok.

By the way: the output `Neightbhorhoods/GeneList.bed` is very similar to the whitelist `RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed`. Same columns, same line count, but genes loci are modified so that gene length are all 500 bp. The center of the gene is taken to be either $2 or $3 (depending upstream or downstream, +/- in $6), and expanded 250 bp on both side.

### Step 3: computing the ABC score

If experimentally derived contact data is not available, one can run the ABC model using the powerlaw estimate only. In this case the ```--HiCdir``` argument should be excluded from ```predict.py``` and the ```--score_column powerlaw.Score``` argument should be included instead. In this case the ```ABC.Score``` column of the predictions file will be set to ```NaN```. The ```powerlaw.Score``` column of the output prediction files will be the relevant Score column to use.

```bash
python ../src/predict.py \
--enhancers ABC_output/Neighborhoods/EnhancerList.txt \
--genes ABC_output/Neighborhoods/GeneList.txt \
--score_column powerlaw.Score \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType K562 \
--outdir ABC_output/Predictions/ \
--make_all_putative
```

> ```bash
> ABC_output/
> ├── Neighborhoods
> │   ├── EnhancerList.bed
> │   ├── EnhancerList.txt
> │   ├── Enhancers.DHS.ENCFF134DLD.bam.CountReads.bedgraph
> │   ├── Enhancers.DHS.ENCFF156LGK.bam.CountReads.bedgraph
> │   ├── Enhancers.H3K27ac.ENCFF301TVL.bam.CountReads.bedgraph
> │   ├── Enhancers.H3K27ac.ENCFF879BWC.bam.CountReads.bedgraph
> │   ├── GeneList.bed
> │   ├── GeneList.TSS1kb.bed
> │   ├── GeneList.txt
> │   ├── Genes.DHS.ENCFF134DLD.bam.CountReads.bedgraph
> │   ├── Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph
> │   ├── Genes.H3K27ac.ENCFF301TVL.bam.CountReads.bedgraph
> │   ├── Genes.H3K27ac.ENCFF879BWC.bam.CountReads.bedgraph
> │   ├── Genes.TSS1kb.DHS.ENCFF134DLD.bam.CountReads.bedgraph
> │   ├── Genes.TSS1kb.DHS.ENCFF156LGK.bam.CountReads.bedgraph
> │   ├── Genes.TSS1kb.H3K27ac.ENCFF301TVL.bam.CountReads.bedgraph
> │   └── Genes.TSS1kb.H3K27ac.ENCFF879BWC.bam.CountReads.bedgraph
> ├── Peaks
> │   ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak
> │   ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_model.r
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK_ENCFF134DLD.bam.Counts.bed
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls
> │   ├── ENCFF156LGK_ENCFF134DLD.macs2_summits.bed
> │   └── params.txt
> └── Predictions
>     ├── EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz
>     ├── EnhancerPredictionsAllPutative.txt.gz
>     ├── EnhancerPredictions.bedpe
>     ├── EnhancerPredictionsFull.txt
>     ├── EnhancerPredictions.txt
>     ├── GenePredictionStats.txt
>     └── parameters.predict.txt
> ```

## Results

### Summary statistics

We create a `sumstats` sub-directory in `ABC_output/Predictions` then compute summary statistics using `bedpe.sumstats.sh` script.

```bash
cd ABC_output/Predictions/
gzip -c EnhancerPredictions.bedpe > EnhancerPredictions.bedpe.gz
cd sumstats/
../../../../../Scripts/bedpe.sumstats.sh ../EnhancerPredictions.bedpe.gz "0.02-1" "500-500"
```

Or replace `0.02-1` with `<choose a min>-<choose a max>` where `min` and `max` can be found with

```bash
awk 'BEGIN{FS="\t"; max=0;} {if($8>max){max=$8}} END{print(max);}' ../EnhancerPredictions.bedpe

awk 'BEGIN{FS="\t"; min=1;} {if($8<min){min=$8}} END{print(min);}' ../EnhancerPredictions.bedpe
```

![Density of Enhancer predictions wrt Distance](/home/hoellinger/Documents/INSERM/shared/notes_perso/docs/ABC_model/results1_K562/Distance.png)

*Density of Enhancer predictions wrt Distance.*

![Density of Enhancer predictions wrt Distance + score.quartile](/home/hoellinger/Documents/INSERM/ABC-Enhancer-Gene-Prediction/K562/ABC_output/Predictions/sumstats/Distance_by_score.quantile.density.png)

*Density of Enhancer predictions wrt Distance + score.quartile.*



![](/home/hoellinger/Documents/INSERM/shared/notes_perso/docs/ABC_model/results1_K562/refelt.scorequantile.nbconn.nbtimes.png)

*Number of connections of element 1 / of element 2. Here elt 1 are  enhancers and elt 2 are TSS: a TSS usually makes more connections to  enhancers, than enhancers make connections to TSS (the four screens distinguish between the score-quartiles)*



![Density of Enhancer Predictions wrt Score](/home/hoellinger/Documents/INSERM/shared/notes_perso/docs/ABC_model/results1_K562/Score.png)

*Density of Enhancer Predictions wrt Score.*

