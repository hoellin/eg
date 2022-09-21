# Trying ABC pipeline on a small example

We just verify on our own the pipeline in README, on the very same example, but starting without the output.

## Running the ABC model

First of all we activate the `abcmodel` environment: `conda activate abcmodel`.

Running the ABC model consists of the following steps:

 1. Define candidate enhancer regions
 2. Quantify enhancer activity
 3. Compute ABC Scores

### Step 1: defining candidate elements

`macs2` requires Python 2, so we installed it on a separate conda environment `py2`.

```shell
conda activate py2
macs2 callpeak \
-t test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
-n wgEncodeUwDnaseK562AlnRep1.chr22.macs2 \
-f BAM \
-g hs \
-p .1 \
--call-summits \
--outdir test_chr22/ABC_output/Peaks/
```

This creates the following in ABC_output/Peaks (which was empty before):

> ```bash
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.xls
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_summits.bed
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_model.r
> ```
>

For information purpose, `input_data` and `reference` contain the following:

> ```shell
> $ l test_chr22/input_data
> Chromatin/  Expression/  HiC/
> $ l test_chr22/input_data/Chromatin/
> ENCFF384ZZM.chr22.bam
> ENCFF384ZZM.chr22.bam.bai
> wgEncodeUwDnaseK562AlnRep1.chr22.bam
> wgEncodeUwDnaseK562AlnRep1.chr22.bam.bai
> wgEncodeUwDnaseK562AlnRep2.chr22.bam
> wgEncodeUwDnaseK562AlnRep2.chr22.bam.bai
> wgEncodeUwDnaseK562.mergedPeaks.slop175.withTSS500bp.chr22.bed
> $ l test_chr22/input_data/Expression/
> K562.ENCFF934YBO.TPM.txt
> $ l test_chr22/input_data/HiC/raw/chr22/
> chr22.bedpe*  chr22.bedpe.gz*  chr22.KRnorm.gz  chr22.KRobserved.gz
> $ l test_chr22/input_data/HiC/raw/powerlaw/
> hic.mean_var.txt  hic.powerlaw.txt
> $ l test_chr22/reference/
> chr22      RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed
> chr22.bed  RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.chr22.bed
> ```
>

Now we sort narrowPeak file and invoke makeCandidateRegions.

```bash
conda deactivate #after that we must be in abcmodel conda environment
#Sort narrowPeak file
bedtools sort -faidx test_chr22/reference/chr22 -i test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak > test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted

python src/makeCandidateRegions.py \
--narrowPeak test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted \
--bam test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
--outDir test_chr22/ABC_output/Peaks/ \
--chrom_sizes test_chr22/reference/chr22 \
--regions_blacklist reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--regions_whitelist test_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.chr22.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 3000
```

We should insert here some information about blacklists and **whitelists** for clarification.

> ```bash
> $ l test_chr22/ABC_output/Peaks/
> params.txt
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_model.r
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.candidateRegions.bed
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.wgEncodeUwDnaseK562AlnRep1.chr22.bam.Counts.bed
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.xls
> wgEncodeUwDnaseK562AlnRep1.chr22.macs2_summits.bed
> ```
>

### Step 2: quantifying enhancer activity

```bash
python src/run.neighborhoods.py \
--candidate_enhancer_regions test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes test_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed \
--H3K27ac test_chr22/input_data/Chromatin/ENCFF384ZZM.chr22.bam \
--DHS test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam,test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep2.chr22.bam \
--expression_table test_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
--chrom_sizes test_chr22/reference/chr22 \
--ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir test_chr22/ABC_output/Neighborhoods/
```

Results in the following:

> ```bash
> [...]
> BEDTools completed successfully. 
> 
> Feature DHS completed in 1.1089489459991455
> Assigning classes to enhancers
> Total enhancers: 3330
>          Promoters: 572
>          Genic: 1553
>          Intergenic: 1205
> Neighborhoods Complete!
> ```

and in a new output directory `Neighborhoods` containing the following:

> ```bash
> $ l test_chr22/ABC_output/Neighborhoods/
> EnhancerList.bed
> EnhancerList.txt
> Enhancers.DHS.wgEncodeUwDnaseK562AlnRep1.chr22.bam.CountReads.bedgraph
> Enhancers.DHS.wgEncodeUwDnaseK562AlnRep2.chr22.bam.CountReads.bedgraph
> Enhancers.H3K27ac.ENCFF384ZZM.chr22.bam.CountReads.bedgraph
> GeneList.bed
> GeneList.TSS1kb.bed
> GeneList.txt
> Genes.DHS.wgEncodeUwDnaseK562AlnRep1.chr22.bam.CountReads.bedgraph
> Genes.DHS.wgEncodeUwDnaseK562AlnRep2.chr22.bam.CountReads.bedgraph
> Genes.H3K27ac.ENCFF384ZZM.chr22.bam.CountReads.bedgraph
> Genes.TSS1kb.DHS.wgEncodeUwDnaseK562AlnRep1.chr22.bam.CountReads.bedgraph
> Genes.TSS1kb.DHS.wgEncodeUwDnaseK562AlnRep2.chr22.bam.CountReads.bedgraph
> Genes.TSS1kb.H3K27ac.ENCFF384ZZM.chr22.bam.CountReads.bedgraph
> ```

### Step 3: computing the ABC score

```bash
python src/predict.py \
--enhancers test_chr22/ABC_output/Neighborhoods/EnhancerList.txt \
--genes test_chr22/ABC_output/Neighborhoods/GeneList.txt \
--HiCdir test_chr22/input_data/HiC/raw/ \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType K562 \
--outdir test_chr22/ABC_output/Predictions/ \
--make_all_putative
```

> ```bash
> $ l test_chr22/ABC_output/Predictions/
> EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz  EnhancerPredictions.txt
> EnhancerPredictionsAllPutative.txt.gz                   GenePredictionStats.txt
> EnhancerPredictions.bedpe                               parameters.predict.txt
> EnhancerPredictionsFull.txt
> ```

# Test ABC pipeline through SSH

## Load environments

Load the `abcmodel` environment from `base` and load all required modules.

```bash
conda activate base && conda activate abcmodel && module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
```

### Using `macs2`

Use the conda environment `py2` and load the `system/Python-2.7.2` module whenever you need to use `macs2`.

```shell
conda activate py2 && module load system/Python-2.7.2  # comes with macs2
```

Remember to deactivate the `py2` environment **and to unload the `Python2` module** when done with `macs2`.

```bash
conda deactivate && module unload system/Python-2.7.2
```

## Create `test_chr22`

```bash
cp -r example_chr22/ test_chr22/ && rm -r test_chr22/ABC_output/
```

Then follow the above pipeline (later this markdown will be included at the end of `test_abc_pipeline.md`).



## Run the ABC model

Don't forget to use `srun --pty bash` to connect to a node before processing.

### Step 1: define candidate elements

```bash
conda activate base && conda activate abcmodel && module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6 && conda activate py2 && module load system/Python-2.7.2
```

```shell
srun --pty bash
macs2 callpeak \
-t test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
-n wgEncodeUwDnaseK562AlnRep1.chr22.macs2 \
-f BAM \
-g hs \
-p .1 \
--call-summits \
--outdir test_chr22/ABC_output/Peaks/
```

```bash
conda deactivate && module unload system/Python-2.7.2
```

```bash
#First sort narrowPeak file
srun --pty bash
bedtools sort -faidx test_chr22/reference/chr22 -i test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak > test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted

#Then
srun --pty bash
python src/makeCandidateRegions.py \
--narrowPeak test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted \
--bam test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
--outDir test_chr22/ABC_output/Peaks/ \
--chrom_sizes test_chr22/reference/chr22 \
--regions_blacklist reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--regions_whitelist test_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.chr22.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 3000
```

### Step 2: quantifying enhancer activity

```bash
srun --pty bash
python src/run.neighborhoods.py \
--candidate_enhancer_regions test_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes test_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed \
--H3K27ac test_chr22/input_data/Chromatin/ENCFF384ZZM.chr22.bam \
--DHS test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam,test_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep2.chr22.bam \
--expression_table test_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
--chrom_sizes test_chr22/reference/chr22 \
--ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir test_chr22/ABC_output/Neighborhoods/
```

### Step 3: computing the ABC score

```bash
srun --pty bash 
python src/predict.py \
--enhancers test_chr22/ABC_output/Neighborhoods/EnhancerList.txt \
--genes test_chr22/ABC_output/Neighborhoods/GeneList.txt \
--HiCdir test_chr22/input_data/HiC/raw/ \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType K562 \
--outdir test_chr22/ABC_output/Predictions/ \
--make_all_putative
```

