# THIS MARKDOWN CONTAINS GUIDELINES THAT LEAD TO AN ERROR AT THE END OF STEP 1, AND THE LOG OF HOW THE ERROR HAS BEEN FOUND. FOR ARCHIVE ONLY.

# ABC model on K562

Adapted from [ABC model's README on github](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

The Activity-by-contact model predicts which enhancers regulate which genes on a cell type specific basis.

## Characteristics

* No whitelist
* Use both replicates whenever possible

## Requirements

For each cell-type, the inputs to the ABC model are:

 * Required Inputs

   * bam file for Dnase-Seq or ATAC-Seq (indexed and sorted)

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

Experiment: `ENCSR000EOT` [K562 DNase-seq whole genome on GRCh38 genome assembly](https://www.encodeproject.org/search/?type=Experiment&status=released&perturbed=false&assay_title=DNase-seq&biosample_ontology.term_name=K562&assembly=GRCh38&files.file_type=bam&lab.title=John+Stamatoyannopoulos%2C+UW).

Replicate 1: `ENCFF156LGK` ; Replicate 2:  `ENCFF134DLD`

```bash
cd K562/
wget https://www.encodeproject.org/files/ENCFF156LGK/@@download/ENCFF156LGK.bam -P reference/dnase/
wget https://www.encodeproject.org/files/ENCFF156LGK/@@download/ENCFF134DLD.bam -P reference/dnase/
```

### Blacklist

We use the ENCODE Blacklist for Human GRCh38 downloaded from [this link](https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz).

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

We use the same as usual:

```bash
cd K562/
cp ../reference/chr_size* reference/
```

## Data processing

At the end of the day we should have the following in `K562/reference/`:

> ```
> reference/
> ├── blacklist
> │   └── hg38-blacklist.v2.bed.gz
> ├── chr_sizes
> ├── chr_sizes.bed
> ├── dnase
> │   ├── ENCFF134DLD.bam
> │   └── ENCFF156LGK.bam
> ├── expression
> │   ├── ENCFF172GIN.tsv
> │   ├── ENCFF768TKT.tsv
> │   ├── K562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt
> │   ├── v1_K562.ENCFF172GIN_ENCFF768TKT.TPM.txt
> │   └── v2_K562.ENCFF172GIN_ENCFF768TKT.TPM.txt
> ├── genes_names
> │   └── gencode.v35.annotation.gnid.gnname.tsv
> └── H3K27ac
>     ├── ENCFF301TVL.bam
>     └── ENCFF879BWC.bam
> ```

### Create expression table

We first compared the resulting expression table using the 1st or the 2nd replicate. This resulted in significant differences in TPM values, but all 55777 genes listed where the same.

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

We write the following in `compute_mean_expression.awk`.

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
				print("well we should not be here");
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

Then :

```bash
awk -f compute_mean_expression.awk reference/genes_names/gencode.v35.annotation.gnid.gnname.tsv reference/expression/ENCFF172GIN.tsv reference/expression/ENCFF768TKT.tsv > reference/expression/K562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt
```

## Run the ABC model

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

Seems to have worked but the script had not stopped by itself (we had to ctrl+c it). Size delta between files is concordant with the one obtained on test_chr22 though, and content seems ok => it took 25 minutes. Output `foo.macs2_peaks.narrowPeak` contains peaks of various lengths.

<span style="color:red">**So it seems it worked nonetheless but maybe some error has been induced here.**</span>

> ```
> INFO  @ Tue, 01 Dec 2020 17:57:35:  29000000 
> INFO  @ Tue, 01 Dec 2020 17:57:37:  30000000 
> INFO  @ Tue, 01 Dec 2020 17:57:43: #1 tag size is determined as 36 bps 
> INFO  @ Tue, 01 Dec 2020 17:57:43: #1 tag size = 36.0 
> INFO  @ Tue, 01 Dec 2020 17:57:43: #1  total tags in treatment: 195415953 
> INFO  @ Tue, 01 Dec 2020 17:57:43: #1 user defined the maximum tags... 
> INFO  @ Tue, 01 Dec 2020 17:57:43: #1 filter out redundant tags at the same location and the same strand by allowing at most 1 tag(s) 
> INFO  @ Tue, 01 Dec 2020 17:57:46: #1  tags after filtering in treatment: 88526164 
> INFO  @ Tue, 01 Dec 2020 17:57:46: #1  Redundant rate of treatment: 0.55 
> INFO  @ Tue, 01 Dec 2020 17:57:46: #1 finished! 
> INFO  @ Tue, 01 Dec 2020 17:57:46: #2 Build Peak Model... 
> INFO  @ Tue, 01 Dec 2020 17:57:46: #2 looking for paired plus/minus strand peaks... 
> INFO  @ Tue, 01 Dec 2020 17:57:56: #2 number of paired peaks: 100801 
> INFO  @ Tue, 01 Dec 2020 17:57:56: start model_add_line... 
> INFO  @ Tue, 01 Dec 2020 17:57:57: start X-correlation... 
> INFO  @ Tue, 01 Dec 2020 17:57:57: end of X-cor 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #2 finished! 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #2 predicted fragment length is 60 bps 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #2 alternative fragment length(s) may be 60 bps 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #2.2 Generate R script for model : ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_model.r 
> WARNING @ Tue, 01 Dec 2020 17:57:57: #2 Since the d (60) calculated from paired-peaks are smaller than 2*tag length, it may be influenced by unknown sequencing problem! 
> WARNING @ Tue, 01 Dec 2020 17:57:57: #2 You may need to consider one of the other alternative d(s): 60 
> WARNING @ Tue, 01 Dec 2020 17:57:57: #2 You can restart the process with --nomodel --extsize XXX with your choice or an arbitrary number. Nontheless, MACS will continute computing. 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #3 Call peaks... 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #3 Going to call summits inside each peak ... 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #3 Call peaks with given -log10pvalue cutoff: 1.00000 ... 
> INFO  @ Tue, 01 Dec 2020 17:57:57: #3 Pre-compute pvalue-qvalue table... 
> INFO  @ Tue, 01 Dec 2020 18:00:50: #3 Call peaks for each chromosome... 
> INFO  @ Tue, 01 Dec 2020 18:06:39: #4 Write output xls file... ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls 
> INFO  @ Tue, 01 Dec 2020 18:06:47: #4 Write peak in narrowPeak format file... ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak 
> ^CUser interrupted me! ;-) Bye!
> ```

```bash
conda deactivate
conda activate base && conda activate abcmodel
```

May need to install `juicer` too (in `abcmodel`). `samtools` and `bedtools` are required but already installed in `abcmodel` environment.

#### Use ABC `makeCandidateRegions.py` to define candidate regions

> Now `makeCandidateRegions.py` will take as input the narrowPeak file produced by MACS2 and then perform the following processing steps:
> 
> 1. Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
> 2. Resize each of these N peaks to be a fixed number of base pairs centered on the peak summit
> 3. Remove any blacklisted regions and include any whitelisted regions
> 4. Merge any overlapping regions

**This time we do not use any whitelist.** Although initially there were ~ 25 000 whitelisted regions. Later we shall verify whether this has much influence or not on selected candidate regions.

```bash
#First sort narrowPeak file
bedtools sort -faidx reference/chr_sizes -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
```

Did not work! `.sorted` file is empty. But others files seem ok. It yielded the following error:

> ```bash
> Chromosome "chr11_KI270721v1_random" undefined in reference/chr_sizes
> ```

Indeed, we find that 112 scaffolds are referenced in `ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak` - not only chromosomes (`awk '{print $1}' ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak |uniq |wc -l`).

We decide to keep chromosomes only by modifying `ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak` (1321235 lines initially) and `ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls` (1321263 lines). But we have to keep in mind that neither that `ENCFF156LGK_ENCFF134DLD.macs2_model.r` won't be affected and will not mean anything anymore.

Maybe we better have kept whole `ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak` and filtered using another more complete "chrom_sizes" list, because later we see that `makeCandidateRegions.py` performs such a filter at first.

```bash
awk '$1 ~ /^(chr)([1-9]{1,2}$)|(M$)/ {print $1}' ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak |uniq
```

```bash
awk '$1 ~ /^(chr)([1-9]{1,2}$)|(M$)/ {print $0}' ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak
```

> ```bash
> $ wc -l ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak 
> 1200345 ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak
> ```

```bash
awk 'NR <= 28 || $1 ~ /^(chr)([1-9]{1,2}$)|(M$)/ {print $0}' ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
```

> ```bash
> $ wc -l ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls 
> 1200373 ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
> ```

Finally there were 1 321 235 - 1 200 345 = 120 890 DNase peaks on scaffolds others than chromosomes.

No we re-try the following.

```bash
#First sort narrowPeak file
bedtools sort -faidx reference/chr_sizes -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
```

This time it yielded no error and the output is not empty.

> ```bash
> $ wc -l ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
> 1200345 ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
> ```

> ```
> K562
> ├── ABC_output
> │   ├── Neighborhoods
> │   └── Peaks
> │       ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.narrowPeak
> │       ├── ENCFF156LGK_ENCFF134DLD.chromosomes_only.macs2_peaks.xls
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_model.r
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak
> │       ├── ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted
> │       └── ENCFF156LGK_ENCFF134DLD.macs2_peaks.xls
> ├── compute_mean_expression.awk
> └── reference
>     ├── blacklist
>     │   └── hg38-blacklist.v2.bed.gz
>     ├── chr_sizes
>     ├── chr_sizes.bed
>     ├── dnase
>     │   ├── ENCFF134DLD.bam
>     │   └── ENCFF156LGK.bam
>     ├── expression
>     │   ├── ENCFF172GIN.tsv
>     │   ├── ENCFF768TKT.tsv
>     │   ├── K562.ENCFF172GIN_ENCFF768TKT.mean.TPM.txt
>     │   ├── v1_K562.ENCFF172GIN_ENCFF768TKT.TPM.txt
>     │   └── v2_K562.ENCFF172GIN_ENCFF768TKT.TPM.txt
>     ├── genes_names
>     │   └── gencode.v35.annotation.gnid.gnname.tsv
>     └── H3K27ac
>         ├── ENCFF301TVL.bam
>         └── ENCFF879BWC.bam
> ```

```bash
#Then
python ../src/makeCandidateRegions.py \
--narrowPeak ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted \
--bam reference/dnase/ENCFF156LGK.bam,reference/dnase/ENCFF134DLD.bam \
--outDir ABC_output/Peaks/ \
--chrom_sizes reference/chr_sizes \
--regions_blacklist reference/blacklist/hg38-blacklist.v2.bed.gz \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000
#Expected output: params.txt, ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed, ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.ENCFF134DLD.bam.Counts.bed
```

Does not work with 2 `bam` inputs. Let's try with 1.

By the way:

```bash
$ python ../src/makeCandidateRegions.py --help
```

> ```bash
> --regions_whitelist REGIONS_WHITELIST
>                         Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blacklist (default: )
> ```

```bash
python ../src/makeCandidateRegions.py \
--narrowPeak ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted \
--bam reference/dnase/ENCFF156LGK.bam \
--outDir ABC_output/Peaks/ \
--chrom_sizes reference/chr_sizes \
--regions_blacklist reference/blacklist/hg38-blacklist.v2.bed.gz \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000
#Expected output: params.txt, ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed, ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed
```

Seemed ok but it appears that output `ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed` is empty... (the other `.bed` output is not empty at all though...) Full log:

```
awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' reference/chr_sizes <(samtools view -H reference/dnase/ENCFF156LGK.bam | grep SQ | cut -f 2 | cut -c 4- )  > temp && bedtools sort -faidx temp -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted |wc -l $$ rm temp
```



> ```
> Running: awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' reference/chr_sizes <(samtools view -H reference/dnase/ENCFF156LGK.bam | grep SQ | cut -f 2 | cut -c 4- )  > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed.temp_sort_order
> Running: bedtools sort -faidx ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed.temp_sort_order -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted | bedtools coverage -g ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed.temp_sort_order -counts -sorted -a stdin -b reference/dnase/ENCFF156LGK.bam | awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}'  | bedtools sort -faidx reference/chr_sizes -i stdin > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed; rm ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed.temp_sort_order
> BEDTools completed successfully. 
> 
> Running: bedtools sort -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed -faidx reference/chr_sizes | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n 150000 |bedtools intersect -b stdin -a ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted -wa |awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' |bedtools slop -i stdin -b 250 -g reference/chr_sizes |bedtools sort -i stdin -faidx reference/chr_sizes |bedtools merge -i stdin | bedtools intersect -v -wa -a stdin -b reference/blacklist/hg38-blacklist.v2.bed.gz | cut -f 1-3 | bedtools sort -i stdin -faidx reference/chr_sizes | bedtools merge -i stdin > ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed
> ```

So looking at the log plus at the following, it appears that the first step (computing `....bam.Counts.bed`) worked correctly:

> ```bash
> $ head ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed 
> chr1	10132	10192	9
> chr1	10293	10367	18
> chr1	16205	16300	75
> chr1	30837	30910	11
> chr1	101802	101862	12
> chr1	104941	105028	80
> chr1	106777	106857	17
> chr1	109125	109219	31
> chr1	115612	115828	1425
> chr1	116373	116445	20
> ```

Meaning the problem lies in the 2nd step.

Maybe, instead of using only 1 `bam` , we should have concatenated the two ones (as we've used both previously)?

**1st step**

* `samtools view -H reference/dnase/ENCFF156LGK.bam | grep SQ | cut -f 2 | cut -c 4-` in the first step (see log) searches the 195 scaffolds (try it with `head` / `wc -l`).
* So `awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' reference/chr_sizes <(samtools view -H reference/dnase/ENCFF156LGK.bam | grep SQ | cut -f 2 | cut -c 4- )` creates a new "chrom_sizes" with only the chromosomes found in both `bam` and `reference/chr_sizes`, in the order they appear in the `bam`, so here it results in the very same as directly taking `reference/chr_sizes`. The result is stored in `foo.Counts.bed.temp_sort_order`.
* The next step re-sort the (> 1 million, 1200345) DNase peaks from `ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted` only taking into account the ones that are within reference chromosomes (so at that steps other scaffolds are removed... so did we really need to remove them manually at first?) => it should not change anything as we've manually performed such a filter before. We tried to execute these first them manually then to meld the files to find differences and surprisingly there are quite much line inversions between the 2 files... not sure why
* Now the result is given as an argument to `bedtools coverage`. It counts the DNase reads in `reference/dnase/ENCFF156LGK.bam` (altogether there are more than 60 million reads) that fall into each of the ~ 1.2 million DNase peaks given as stdin.
* These peaks are stored in `foo.bam.Counts.bed` and the temporary file is deleted.

=> We should definitely have concatenated the BAM together... Indeed, it does not make any sense to use 2 replicates to call peaks, then to use only one to get back the counts in these peaks.

At the steps the peaks are of variable length.

**2nd step:**

* Ensure sorted, then `bedtools merge`. It takes the previous output (1 200 277 peaks with counts), and merge overlapping peaks, resulting in 1 062 773 peaks.
* The result is sorted by counts (descending order) and only the first N=150 000 peaks are kept
* Now with `bedtools intersect`, instead of these N = 150 000 top peaks among merged peaks, we consider all the **unmerged** peaks, from, `ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted`, that overlap the merged peaks => resulting in 273 588 peaks (among ~ 1.2 millions).
* Only a few columns of these unmerged peaks are kept
* `bedtools slop`: The peaks are resized to 500 bp, using again the reference genome file `chrom_sizes`. So at this time there are 273 588 peaks.
* Peaks are sorted again. But at this precise step there is an error, resulting in 0 peaks for all downstream steps...

> ```bash
> $ bedtools sort -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed -faidx reference/chr_sizes | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n 150000 |bedtools intersect -b stdin -a ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted -wa |awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' |bedtools slop -i stdin -b 250 -g reference/chr_sizes |bedtools sort -i stdin -faidx reference/chr_sizes |wc -l
> Error: malformed BED entry at line 73598. Start was greater than end. Exiting.
> ```

Indeed:

> ```bash
> $ bedtools sort -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed -faidx reference/chr_sizes | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n 150000 |bedtools intersect -b stdin -a ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted -wa |awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' |bedtools slop -i stdin -b 250 -g reference/chr_sizes |awk 'NR == 73598 {print $0}'
> chr3	198022430	198022429
> ```
>
> ```bash
> $ bedtools sort -i ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed -faidx reference/chr_sizes | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n 150000 |bedtools intersect -b stdin -a ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted -wa |awk 'NR == 73598 {print $0}'
> chr3	198027103	198027363	ENCFF156LGK_ENCFF134DLD.macs2_peak_757594	344	.	7.88043	34.48334	31.73137	128
> ```

> ```bash
> $cat reference/chr_sizes
> chr1	249250621
> chr2	243199373
> chr3	198022430
> chr4	191154276
> ```

Okay the problem is that our reference chr3 is too short !!! "198027103" cannot be a position on chr3 wrt our reference!

We are going to redo all the previous step working with reference chromosomes sizes obtained from the header of the DNase-seq `.bam` file:

```bash
samtools view -H reference/dnase/ENCFF156LGK.bam |head
```

We keep this markdown with the mistake for future debugging purpose, and we create another one!

### Step 2: quantifying enhancer activity

(One must first ensure there exists a copy of `chr_sizes` but with extension `.bed`).

**We tried to let the `genes` argument empty but this led to the following error. Initially the whitelist was given as `gene` argument.**

> ```bash
> run.neighborhoods.py: error: the following arguments are required: --genes
> ```

`--help` yields:

>
> ```bash
> --genes GENES         bed file with gene annotations. Must be in bed-6 format. Will be used to assign TSS to genes. (default: None)
> ```

<span style="color:red">**Although it seems that `../reference/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed` is nothing else but gene annotation of the TSS in the human genome. So why is it called a `whitelist` in step 1 ?**</span>

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

Failed.

> ```bash
> ...
> BEDTools completed successfully. 
> 
> samtools idxstats: fail to load index for "reference/dnase/ENCFF134DLD.bam", reverting to slow method
> Feature DHS completed in 455.5042428970337
> Traceback (most recent call last):
>   File "../src/run.neighborhoods.py", line 97, in <module>
>     main(args)
>   File "../src/run.neighborhoods.py", line 93, in main
>     processCellType(args)
>   File "../src/run.neighborhoods.py", line 77, in processCellType
>     load_enhancers(genes=genes_for_class_assignment, 
>   File "/home/hoellinger/Documents/INSERM/ABC-Enhancer-Gene-Prediction/src/neighborhoods.py", line 193, in load_enhancers
>     enhancers = read_bed(candidate_peaks)
>   File "/home/hoellinger/Documents/INSERM/ABC-Enhancer-Gene-Prediction/src/neighborhoods.py", line 479, in read_bed
>     skip = 1 if ("track" in open(filename, "r").readline()) else 0
> FileNotFoundError: [Errno 2] No such file or directory: 'ABC_output/Peaks/ENCFF156LGK_ENCFF134DLD.macs2_peaks.narrowPeak.sorted.candidateRegions.bed'
> ```



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

## Results

### Summary statistics

We create a `sumstats` sub-directory in `ABC_output/Predictions` then compute summary statistics using `bedpe.sumstats.sh` script.

```bash
cd ABC_output/Predictions/sumstats/
../../../../../scripts/bedpe.sumstats.sh ../EnhancerPredictions.bedpe.gz "0.02-1" "500-500"
```


