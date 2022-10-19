
# ABC model over GM12878 with same accessions as for BENGI benchmarks

## How to use this notebook?
First, make a copy of [this notebook](../../../scripts/notebooks/ipynb/ABC/BENGI_GM12878.ipynb) on your computer / cluster.

Then, to use this notebook, one should only have to carefully modify the content of the [Set working directory](#set-working-directory) section, then to execute the notebook cell by cell, in the correct order. After execution of each cell, remember to check for errors before executing the next one.

One can use the following to switch notebook theme:


```python
# night theme
#!jt -t monokai -f fira -fs 10 -nf ptsans -nfs 11 -N -kl -cursw 2 -cursc r -cellw 95% -T
```


```python
# standard theme
#!jt -r
```

## Import required librairies


```python
import os.path #Or see https://stackoverflow.com/a/82852 for an object-oriented approach
from IPython.core.magic import register_line_cell_magic
```


```python
@register_line_cell_magic
def writetemplate(line, cell):
    with open(line, 'w') as f:
        f.write(cell.format(**globals()))
```

## Set working directory


```python
mail_user = "tristan.hoellinger@inserm.fr" # slurm notifications
version = "hg19" # used only in the current cell to compute some paths 
cell_type = "gm12878" # used all along the notebook
specie = "homo_sapiens/"+version+"/" # used only in the current cell to compute some paths 

# Where to store run-specific references, scripts, intermediate files, predictions, etc
work_dir = "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/BENGI_GM12878/"

scripts = work_dir+"hi_slurm/" # where to store scripts specific to this run
results_dir = "/work2/project/regenet/results/" # used to compute paths where data for this cell type is
                                                # stored / will be downloaded
dnase_dir = results_dir+"dnaseseq/"+specie+cell_type+'/'
chipseq_dir = results_dir+"chipseq/h3k27ac/"+specie+cell_type+'/'
expression_dir = results_dir+"rnaseq/"+specie+cell_type+'/'
blacklist_dir = "/work2/project/regenet/data/species/"+specie # where the ENCODE blacklist is (going to be) stored

blacklist = blacklist_dir+version+"-blacklist.v2.bed.gz" # need not exist yet

# Gene annotation in the gtf format and gene name / gene id table (1st column id, 2nd column name)
gene_annotation = "/work2/project/fragencode/data/species/homo_sapiens/hg19.gencv19/homo_sapiens.gtf"
gnid_gname = "/work2/project/fragencode/data/species/homo_sapiens/hg19.gencv19/homo_sapiens.gnid.gnname.tsv"

# Accessions for the current run. Just put accession numbers here, no need to download these accessions "by hand".
dnase_rep1 = "ENCFF664VNB" # exp ENCSR000COQ
dnase_rep2 = "ENCFF790WER" # exp ENCSR000COQ
h3k27_rep1 = "ENCFF197QHX" # exp ENCSR000AKC
h3k27_rep2 = "ENCFF882PRP" # exp ENCSR000AKC
rnaseq_rep1 = "ENCFF587ALK" # exp ENCSR000COQ
rnaseq_rep2 = "ENCFF766CPS" # exp ENCSR000COQ

# Where to download the blacklist. Will not be used unless blacklist not found yet.
blacklist_link = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"

dnase_file_rep1 = dnase_dir+dnase_rep1+".bam" # do not change
dnase_file_rep2 = dnase_dir+dnase_rep2+".bam"  # do not change
h3k27_file_rep1 = chipseq_dir+h3k27_rep1+".bam"  # do not change
h3k27_file_rep2 = chipseq_dir+h3k27_rep2+".bam"  # do not change
rnaseq_file_rep1 = expression_dir+rnaseq_rep1+".tsv"  # do not change
rnaseq_file_rep2 = expression_dir+rnaseq_rep2+".tsv"  # do not change

reference_dir = work_dir+"reference/" # do not change
annotations_dir = reference_dir+"gene_annotation/" # do not change
peaks = work_dir+"ABC_output/Peaks/" # do not change
neighborhoods = work_dir+"ABC_output/Neighborhoods/" # do not change
predictions = work_dir+"ABC_output/Predictions/" # do not change

light_annotation = annotations_dir+"gene_ids.bed"

# Ubiquitously expressed genes (gene names). If provided, `ubiquitously_expressed_genes` will be
# automatically generated with gene ids instead of gene names. If not (left empty), 
# `ubiquitously_expressed_genes` must be provided.
ubiquitous_gene_names = work_dir+"../reference/UbiquitouslyExpressedGenesHG19.txt"
# Ubiquitously expressed genes (gene ids). Please provide either `ubiquitously_expressed_genes` or
# `ubiquitous_gene_names`
ubiquitously_expressed_genes = reference_dir+"UbiquitouslyExpressedGenes_ids.txt" # need not exist

compute_mean_expression = work_dir+"../compute_mean_expression.awk"
# compute_mean_expression.awk computes an expression table (<gene id> <expression>) taking the mean
# TPM expression found in the 2 input polyA+ RNA-seq tsv.

make_candidate_regions = work_dir+"../src/makeCandidateRegions.py"
run_neighborhoods = work_dir+"../src/run.neighborhoods.py"
predict = work_dir+"../src/predict.py"

genome_file = reference_dir+"chr_sizes" # shall not exist yet

gene_expression_table = reference_dir+"expression/"+cell_type+"."+rnaseq_rep1+"_"+rnaseq_rep2+".mean.TPM.txt"
```


```python
if not os.path.isdir(reference_dir):
    !mkdir $reference_dir
if not os.path.isdir(scripts):
    !mkdir $scripts
```


```python
# One may only change what is between quotes here
slurm_1_1 = scripts+"step1.1.sh"
slurm_1_2 = scripts+"step1.2.sh"
slurm_1_3 = scripts+"step1.3.sh"
slurm_2 = scripts+"step2.sh"
slurm_3 = scripts+"step3.sh"
```


```python
ccREs_dir = results_dir+"ccREs/"+specie+"V1/"+cell_type+'/'

ccREs_accession = "ENCFF028SGJ" # ENCODE ccREs with annotations for GM12878 cell type
TSS_list = "/work2/project/regenet/workspace/thoellinger/BENGI/Benchmark/Annotations/GENCODEv19-TSSs.bed"

ccREs = ccREs_dir+ccREs_accession+".bed" # need not exist yet
ccREs_ELS = ccREs_dir+cell_type+"_ccRE_ELSs.bed" # need not exist yet
ccREs_dELS = ccREs_dir+cell_type+"_ccRE_dELSs.bed" # need not exist yet
```

## Data acquisition

### Chromatin accessibility (DNase-seq)

We use the same DNase-seq data as Moore et al. 2020 for DNase-DNase correlation method - that is, the data used by Thurman et al. 2012.

Experiment: `ENCSR000COQ` https://www.encodeproject.org/experiments/ENCSR000COQ/

Replicate 1: `ENCFF664VNB.bam` ; Replicate 2: `ENCFF790WER.bam`


```python
if not os.path.isfile(dnase_file_rep1):
    !wget https://www.encodeproject.org/files/$dnase_rep1/@@download/$dnase_rep1$".bam" -P $dnase_dir
            
if not os.path.isfile(dnase_file_rep2):
    !wget https://www.encodeproject.org/files/$dnase_rep2/@@download/$dnase_rep2$".bam" -P $dnase_dir
```

### Blacklist

We use [the ENCODE Blacklist](https://github.com/Boyle-Lab/Blacklist) for Human GRCh19 downloaded from [this link](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz).


```python
if not os.path.isfile(blacklist):
    !wget $blacklist_link -P $blacklist_dir
```

### Histone mark H3K27ac ChIP-seq

In Moore et al. 2020, histone marks are used:

* to build the ccREs
* for the training of Target Finder

For the sake of comparability between methods evaluation, if ever it appears that H3K27ac ChIP-seq accessions are not the same for ccREs and Target Finder training, we should use here the same as the ones used for Target Finder training. We shall verify later whether the accessions are the same or not.

In Moore et al. 2020:

> ***TargetFinder***
> We reimplemented TargetFinder to run on our cCRE-gene pairs with chromCV. For features, we used the identical datasets described by Whalen et al. for each cell type.

The accessions of H3k27ac ChIP-seq files used are available in supplementary Table S2 of [Whalen et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4910881/):

|     Cell line      |   K562   | GM12878  | HeLa-S3  |  HUVEC   |  IMR90   |   NHEK   |
| :----------------: | :------: | -------- | -------- | :------: | :------: | :------: |
| ChIP-seq (H3K27ac) | GSE29611 | GSE29611 | GSE29611 | GSE29611 | GSE16256 | GSE29611 |

So either we download the big file (> 330 gb) for the `GSE29611` accession on GEO, or only the info file:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29611/suppl/GSE29611_run_info.txt.gz -P /work2/project/regenet/results/chipseq/h3k27ac/homo_sapiens/hg19/gm12878/
```
in which we find the precise accession for H3K27ac ChIP-seq, which is: `GSM733771`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM733771.
(see file names in "Set working directory" section).

On the [ENCODE portal](https://www.encodeproject.org/) we find that this experiment corresponds to [experiment ENCSR000AKC](https://www.encodeproject.org/experiments/ENCSR000AKC/), and that the accessions available [on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM733771) correspond to [ENCFF197QHX](https://www.encodeproject.org/files/ENCFF197QHX/) and [ENCFF882PRP](https://www.encodeproject.org/files/ENCFF882PRP/). Note that sizes do not match (especially for replicate 1), but this is probably due to specific ENCODE processing pipeline.

Hence we download the following accessions:


```python
if not os.path.isfile(h3k27_file_rep1):
    !wget https://www.encodeproject.org/files/$h3k27_rep1/@@download/$h3k27_rep1$".bam" -P $chipseq_dir
            
if not os.path.isfile(h3k27_file_rep2):
    !wget https://www.encodeproject.org/files/$h3k27_rep2/@@download/$h3k27_rep2$".bam" -P $chipseq_dir
```

### Gene expression (polyA+ RNA-seq)

We use the two replicates provided by BENGI Moore et al. 2020 in Additional File 1 Table S1b.

Experiment: `ENCSR000COQ` Thomas Gingeras, CSHL https://www.encodeproject.org/experiments/ENCSR000COQ/

Replicate 1: `ENCFF587ALK` ; Replicate 2: `ENCFF766CPS`


```python
if not os.path.isfile(rnaseq_file_rep1):
    !wget https://www.encodeproject.org/files/$rnaseq_rep1/@@download/$rnaseq_rep1$".bam" -P $expression_dir
            
if not os.path.isfile(rnaseq_file_rep2):
    !wget https://www.encodeproject.org/files/$rnaseq_rep2/@@download/$rnaseq_rep2$".bam" -P $expression_dir
```

### Gene annotation

We use the encode version 19 of hg19 (`gene_annotation` here). We also need a files with 2 columns making the association between gene names / gene ids (`gnid_gname` here).

### Ubiquitously expressed genes

We use the `UbiquitouslyExpressedGenesHG19.txt` list of 847 ubiquitously expressed genes used by ABC authors. We replace gene symbols (names) with gene ids, which result in a list of 843 ubiquitously expressed genes for which we find an id.


```bash
%%bash -s "$gnid_gname" "$ubiquitous_gene_names" "$ubiquitously_expressed_genes"
if [[ -f $2 ]] && [[ ! -f $3 ]]; then
    awk 'BEGIN{sep="\t"} NR==FNR {id[$2]=$1; next} id[$1] {print id[$1]}' $1 $2 > $3
fi
```


```python
!wc -l $ubiquitously_expressed_genes
```

    843 /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/final_BENGI_GM12878/reference/UbiquitouslyExpressedGenes_ids.txt


### Genome file

We extract chromosomes sizes from the header of `ENCFF664VNB.bam` (DNase-seq, `dnase_rep1`).


```bash
%%bash -s "$work_dir" "$dnase_file_rep1" "$genome_file"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && module load bioinfo/samtools-1.9
samtools view -H $2 | grep SQ | cut -f 2,3 |cut -c 4- |awk 'BEGIN{FS="\t"} {split($2,locus,":"); if($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/){print($1"\t"locus[2])}}' > $3
awk '{print $1"\t"0"\t"$2}' $3 > $3".bed"
```

### ccREs
We  need to download ccREs with annotations for this cell type, in order to intersect ABC predictions with distal ccREs with Enhancer-like signature, to later intersect them with BENGI datasets. The correct accession can be found in Supplentary Figure 1 Table S1a of Moore et al. 2020.


```python
if not os.path.isfile(ccREs+".gz"):
    !wget https://www.encodeproject.org/files/$ccREs_accession/@@download/$ccREs_accession$".bed.gz" -P $ccREs_dir
        
if not os.path.isfile(ccREs):
    !gzip -cd $ccREs".gz" > $ccREs
```

## Data reprocessing

### Create expression table

We use `compute_mean_expression.awk` to take as the reference TPM value, the mean value for the 2 replicates.


```bash
%%bash -s "$compute_mean_expression" "$reference_dir" "$rnaseq_file_rep1" "$rnaseq_file_rep2" "$cell_type" "$rnaseq_rep1" "$rnaseq_rep2"
if [[ ! -d $2"expression/" ]]; then mkdir $2"expression/"; fi
awk -f $1 $3 $4 > $2"expression/"$5"."$6"_"$7".mean.TPM.txt"
```

### Filter gene annotation

We keep genes only, and select appropriate columns to fit the format of `reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed` available on [ABC repository on Github](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction), that is as follows:
> ```bash
> chrX	99883667	99894988	ENSG00000000003.10	0	-
> chrX	99839799	99854882	ENSG00000000005.5	0	+
> chr20	49551404	49575092	ENSG00000000419.8	0	-
> ```

WARNING: ABC assumes gene annotation is sorted at step 2, so we NEED to sort it, otherwise step 2 won't work.

Moreover, we must keep only genes that lie on a chromosome listed in our curated `genome_file`.


```bash
%%bash -s "$work_dir" "$genome_file" "$light_annotation" "$gene_annotation"
if [[ ! -d $3 ]]; then mkdir $3; fi
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && module load bioinfo/bedtools-2.27.1
awk '{if(NR==FNR){!genome[$1]++; next}; if(genome[$1]){if($3 ~ /(^gene$)/){print $1"\t"$4"\t"$5"\t"substr($10,2,length($10)-3)"\t"0"\t"$7}}}' $2 $4 | bedtools sort -g $2 -i stdin > $3
```

## Running the ABC model

### Step 1: define candidate elements

#### Call peaks with `macs2`


```python
%%writetemplate $slurm_1_1
#!/bin/sh                                                                                                                                                                                                                    
# dependencies: python2
macs2 callpeak \
-t {dnase_file_rep1} {dnase_file_rep2} \
-n {dnase_rep1}_{dnase_rep2}.macs2 \
-f BAM \
-g hs \
-p .1 \
--call-summits \
--outdir {peaks}
```


```bash
%%bash -s "$work_dir" "$slurm_1_1" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
conda activate py2 && module load system/Python-2.7.2
sbatch --mem=4G --cpus-per-task=1 -J step1.1 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
conda deactivate && module unload system/Python-2.7.2
```

    Submitted batch job 22128421



```python
!seff 22128421
```

    Job ID: 22128421
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:51:18
    CPU Efficiency: 98.40% of 00:52:08 core-walltime
    Job Wall-clock time: 00:52:08
    Memory Utilized: 1.40 GB
    Memory Efficiency: 35.12% of 4.00 GB


#### Use ABC `makeCandidateRegions.py` to define candidate regions

As there are peaks over many kinds of scaffolds listed in the `.narrowPeak` and `.xls` output files, we modify them keep only chromosomes only (otherwise `bedtools sort` below would not work). **But we have to keep in mind that the `.r` will still corresponds to DNase peaks over all scaffolds**.


```bash
%%bash -s "$peaks" "$dnase_rep1" "$dnase_rep2"
awk '$1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/ {print $0}' $1$2"_"$3".macs2_peaks.narrowPeak" > $1/$2"_"$3".chromosomes_only.macs2_peaks.narrowPeak"
awk 'NR <= 28 || $1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/ {print $0}' $1$2"_"$3".macs2_peaks.xls" > $1$2"_"$3".chromosomes_only.macs2_peaks.xls"
```

Now we sort the result (note that it seems quite tricky to use `srun` from a notebook as we first need to export environment variables - maybe there is another way but as for now I could not find better):


```bash
%%bash -s "$peaks" "$dnase_rep1" "$dnase_rep2" "$genome_file"
export V1=$1 V2=$2 V3=$3 V4=$4
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/bedtools-2.27.1
srun bash
bedtools sort -faidx $V4 -i $V1$V2"_"$V3".chromosomes_only.macs2_peaks.narrowPeak" > $V1$V2"_"$V3".macs2_peaks.narrowPeak.sorted"
```

    srun: job 22131289 queued and waiting for resources
    srun: job 22131289 has been allocated resources


Now as `work_dir ../src/makeCandidateRegions.py` does not support multiple `bam` inputs, we concatenate the bam of the 2 replicates.  Indeed, `makeCandidateRegions.py` will use the bam input to count DNase reads on each peak of the`.macs2_peaks.narrowPeak.sorted` output found with `macs2`.

(WARNING: depending on the size of the input bam, which is much variable, the execution on 1 thread only can take from a few minutes to several hours, so better use at least 8 threads)


```python
n_threads = 8 #we recommend using at least 8 cpu
```


```python
%%writetemplate $slurm_1_2
#!/bin/sh
# dependencies: bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1

if [[ ! -f {dnase_dir}{dnase_rep1}_{dnase_rep2}.bam ]]
then
    if [[ ! -f {dnase_dir}{dnase_rep1}_{dnase_rep2}.sam ]]
    then
        samtools view -@ {n_threads} -h {dnase_file_rep1} > {dnase_dir}{dnase_rep1}_{dnase_rep2}.sam
        samtools view -@ {n_threads} {dnase_file_rep2} >> {dnase_dir}{dnase_rep1}_{dnase_rep2}.sam
        samtools view -@ {n_threads} -S -b {dnase_dir}{dnase_rep1}_{dnase_rep2}.sam > {dnase_dir}{dnase_rep1}_{dnase_rep2}.bam
        rm {dnase_dir}{dnase_rep1}_{dnase_rep2}.sam
    else
        samtools view -@ {n_threads} -S -b {dnase_dir}{dnase_rep1}_{dnase_rep2}.sam > {dnase_dir}{dnase_rep1}_{dnase_rep2}.bam
    fi
fi
```


```bash
%%bash -s "$work_dir" "$slurm_1_2" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9
sbatch --mem=2G --cpus-per-task=8 -J step1.2 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 22134558



```python
!seff 22134558
```

    Job ID: 22134558
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Nodes: 1
    Cores per node: 8
    CPU Utilized: 02:32:39
    CPU Efficiency: 65.72% of 03:52:16 core-walltime
    Job Wall-clock time: 00:29:02
    Memory Utilized: 8.98 MB
    Memory Efficiency: 0.44% of 2.00 GB


Now we can launch `makeCandidateRegions` on slurm:


```python
%%writetemplate $slurm_1_3
#!/bin/sh
# dependencies: python3 bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
python {make_candidate_regions} \
--narrowPeak {peaks}{dnase_rep1}_{dnase_rep2}.macs2_peaks.narrowPeak.sorted \
--bam {dnase_dir}{dnase_rep1}_{dnase_rep2}.bam \
--outDir {peaks} \
--chrom_sizes {genome_file} \
--regions_blacklist {blacklist} \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000
#Expected output: params.txt, foo1_foo2.macs2_peaks.narrowPeak.sorted.candidateRegions.bed, foo1_foo2.macs2_peaks.narrowPeak.sorted.foo1_foo2.bam.Counts.bed
```

WARNING: depending on the length of the inputs, this script may be very RAM-demanding. So we recommend to allocate at least 16 GB RAM to avoid "out of memory" errors.


```bash
%%bash -s "$work_dir" "$slurm_1_3" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=16G --cpus-per-task=1 -J step1.3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 22139064



```python
!seff 22139064
```

    Job ID: 22139064
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:20:54
    CPU Efficiency: 99.68% of 00:20:58 core-walltime
    Job Wall-clock time: 00:20:58
    Memory Utilized: 8.54 GB
    Memory Efficiency: 53.36% of 16.00 GB


### Step 2: quantifying enhancer activity

#### Use ABC `run.neighborhoods.py`


```python
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

WARNING: depending on the length of the inputs, this script may be very RAM-demanding. So we recommend to allocate at least 32 GB RAM in general, to avoid "out of memory" errors.


```bash
%%bash -s "$work_dir" "$slurm_2" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=16G --cpus-per-task=1 -J step2 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 22142267



```python
!seff 22142267
```

    Job ID: 22142267
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 03:41:59
    CPU Efficiency: 99.83% of 03:42:22 core-walltime
    Job Wall-clock time: 03:42:22
    Memory Utilized: 11.83 GB
    Memory Efficiency: 73.94% of 16.00 GB


### Step 3: computing the ABC score

If experimentally derived contact data is not available, one can run the ABC model using the powerlaw estimate only. In this case the ```--HiCdir``` argument should be excluded from ```predict.py``` and the ```--score_column powerlaw.Score``` argument should be included in ```predict.py```. In this case the ```ABC.Score``` column of the predictions file will be set to ```NaN```. The ```powerlaw.Score``` column of the output prediction files will be the relevant Score column to use.


```python
%%writetemplate $slurm_3
#!/bin/sh
python {predict} \
--enhancers {neighborhoods}EnhancerList.txt \
--genes {neighborhoods}GeneList.txt \
--score_column powerlaw.Score \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType {cell_type} \
--outdir {predictions} \
--make_all_putative
```

Again, that step is very RAM-demanding so we recommend to allocate at least 64G of RAM in general.


```bash
%%bash -s "$work_dir" "$slurm_3" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=32G --cpus-per-task=1 -J step3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 22149102



```python
!seff 22149102
```

    Job ID: 22149102
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:29:51
    CPU Efficiency: 99.83% of 00:29:54 core-walltime
    Job Wall-clock time: 00:29:54
    Memory Utilized: 20.07 GB
    Memory Efficiency: 62.70% of 32.00 GB


## Evaluation of ABC predictions on BENGI datasets

### Find ccRE-dELSs


```bash
%%bash -s "$ccREs" "$ccREs_ELS" "$genome_file"
if [[ ! -f $2 || ! -f $2".sorted" ]]
then
    export V1=$1 V2=$2 V3=$3
    awk -v OFS='\t' '$9 ~ /(^255,205,0)/ {print $1, $2, $3, ".", $4, "Enhancer-like"}' $V1 > $V2
    srun bash
    bedtools sort -faidx $V3 -i $V2 > $V2".sorted"
fi
```

Then using `bedtools closest` we keep only distal ccRE-ELS (ie greater than 2 kb from an annotated TSS, GENCODE v19).


```bash
%%bash -s "$genome_file" "$ccREs_ELS" "$TSS_list" "$ccREs_dELS" "$genome_file"
if [[ ! -f $4 ]]
then
    export V1=$1 V2=$2 V3=$3 V4=$4 V5=$5
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun bash
    awk '{if(NR==FNR){chr[$1]++; next}; if(chr[$1]){print $0}}' $V5 $V3 | bedtools sort -faidx $V5 -i stdin | bedtools closest -d -t "first" -g $V1 -a $V2".sorted" -b stdin | awk -v OFS='\t' '$14 >= 2000 {print $1, $2, $3, $4, $5, "distal-ELS"}' > $V4
fi
```

### Data processing


```python
predictions_non_expressed = predictions+"EnhancerPredictionsAllPutativeNonExpressedGenes.txt"
predictions_expressed = predictions+"EnhancerPredictionsAllPutative.txt"
all_predictions = predictions+"AllPredictions.bedpe" # does not exist yet
all_predictions_sorted = all_predictions+".sorted"
all_predictions_ccre_els = all_predictions_sorted+".ccRE_intersect.txt"
all_predictions_ccre_els_bedpe = all_predictions_sorted+".ccRE_intersect.bedpe"

bengi = "/work2/project/regenet/workspace/thoellinger/BENGI/Benchmark/All-Pairs.Natural-Ratio/"
nb_benchmarks = 6
file_names = list()
# benchmark file names without extensions
file_names.append("GM12878.CHiC-Benchmark.v3")
file_names.append("GM12878.CTCF-ChIAPET-Benchmark.v3")
file_names.append("GM12878.GEUVADIS-Benchmark.v3")
file_names.append("GM12878.GTEx-Benchmark.v3")
file_names.append("GM12878.HiC-Benchmark.v3")
file_names.append("GM12878.RNAPII-ChIAPET-Benchmark.v3")

# short custom names for benchmarks, same order as above
names = ["CHiC", "CTCF", "GEUVADIS", "GTEx", "HiC", "RNAPII"]

files_dict = {}
for k in range(nb_benchmarks):
    files_dict[names[k]] = bengi+file_names[k]+".new.bedpe"
```

Unzip ABC predictions.


```python
if not os.path.isfile(predictions_non_expressed):
    !gzip -cd $predictions_non_expressed".gz" > $predictions_non_expressed
if not os.path.isfile(predictions_expressed):
    !gzip -cd $predictions_expressed".gz" > $predictions_expressed    
```

Merge ABC predictions for expressed and non expressed genes, keep only columns of interest and sort the result.


```bash
%%bash -s "$all_predictions" "$predictions_expressed" "$predictions_non_expressed"
if [[ ! -f $1 ]]
then
    awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' $2 > $1;
    tail -n+2 $3 | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' >> $1;
fi
```

The following is very RAM-demanding, we recommend using at least 64GB of ram to avoid "out of memory" errors.


```bash
%%bash -s "$all_predictions" "$all_predictions_sorted" "$genome_file"
if [[ ! -f $V2 ]]
then
    export V1=$1 V2=$2 V3=$3
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun --mem=64G bash
    tail -n+2 $V1 | bedtools sort -faidx $V3 -i stdin  > $V2
fi
```

    srun: job 22250266 queued and waiting for resources
    srun: job 22250266 has been allocated resources



```python
# !head -n 1 $all_predictions_sorted
# COLUMNS DESCRIPTION: (DO NOT DELETE)
# chr, start1 (candidate enhancer), end1 (candidate enhancer), chr, start2 (gene TSS), end2 (gene TSS),
# name (target gene id), score (ABC powerlaw.Score), empty, empty, gene expressed?, gene expression
```

Intersect ABC candidate regions with ccRE-dELS.


```bash
%%bash -s "$ccREs_dELS" "$all_predictions_sorted" "$all_predictions_ccre_els" "$genome_file"
if [[ ! -f $3 ]]
then
    export V1=$1 V2=$2 V3=$3 V4=$4
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun --mem=32G bash
    bedtools intersect -sorted -wo -g $V4 -a $V1 -b $V2 > $V3
fi
```

    srun: job 22250437 queued and waiting for resources
    srun: job 22250437 has been allocated resources


When multiple ABC candidate enhancers for a given gene overlap the same ccRE-ELS, we keep only the one that maximize the ABC score (the "powerlaw.Score" if ABC ran without 3D data).


```bash
%%bash -s "$light_annotation" "$all_predictions_ccre_els" "$all_predictions_ccre_els_bedpe"
if [[ ! -f $3 ]]
then
    export V1=$1 V2=$2 V3=$3
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun --mem=32G bash
    awk '{ \
        if(NR==FNR){ \
            start[$4]=$2; \
            end[$4]=$3; \
            next; \
        }; \
        if(!pair[$5$13]){ \
            pair[$5$13]=$1"\t"start[$13]"\t"end[$13]"\t"$1"\t"$2"\t"$3"\t"$13":"$5"\t"$14"\t"$15"\t"$16"\t"$2"\t"$19"\t"$17"\t"$18; \
            next}; \
        split(pair[$5$13],current_pair,"\t"); \
        if($14>current_pair[$14]){ \
            pair[$5$13]=$1"\t"start[$13]"\t"end[$13]"\t"$1"\t"$2"\t"$3"\t"$13":"$5"\t"$14"\t"$15"\t"$16"\t"$2"\t"$19"\t"$17"\t"$18; \
        }} \
    END{for(u in pair){ \
        print pair[u] \
        }}' $V1 $V2 | bedtools sort -i stdin > $V3
fi
```

    srun: job 22250445 queued and waiting for resources
    srun: job 22250445 has been allocated resources



```python
# !head -n 1 $all_predictions_ccre_els_bedpe
# COLUMNS DESCRIPTION: (DO NOT DELETE)
# chr, start1 (gene), end1 (gene), chr, start2 (ccRE-ELS), end2 (ccRE-ELS),
# name (<target gene id>:<ccRE id>), score (ABC powerlaw.Score), empty, empty, 
# TSS (unique yet), overlap between ABS candidate region and ccRE, gene expressed?,
# gene expression
```


```python
!wc -l $all_predictions_sorted
```

    37308539 /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/BENGI_GM12878/ABC_output/Predictions/AllPredictions.bedpe.sorted



```python
!wc -l $all_predictions_ccre_els
```

    1177503 /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/BENGI_GM12878/ABC_output/Predictions/AllPredictions.bedpe.sorted.ccRE_intersect.txt



```python
!wc -l $all_predictions_ccre_els_bedpe
```

    1084021 /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/BENGI_GM12878/ABC_output/Predictions/AllPredictions.bedpe.sorted.ccRE_intersect.bedpe


### Preliminary remarks

Among the 37,308,539 ABC predictions (positive + negative pairs, expressed + non-expressed genes), 1,084,021 unique pairs involve a ccRE-dELS.

Now we can compute the overlap between these predictions (positives + negatives) and the BENGI datasets. At this point, a small overlap would reflect a bad implementation of the ABC model wrt our objective, rather than poor performance of the ABC model over BENGI datasets.

### Intersection between ABC predictions / BENGI

#### Import files


```python
import pandas as pd
```


```python
col_names = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "is.positive", "strand1", "strand2", "TSSs"] # column names in bedpe benchmarks
col_names_abc = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2", "TSS", "overlap", "is.expressed", "expression"] # column names in abc bedpe
benchmarks = {}
for k in range(nb_benchmarks):
    benchmarks[names[k]] = pd.read_csv(files_dict[names[k]], sep='\t', header=None, usecols=[0, 1, 2, 4, 5, 6, 7, 10], names=col_names, dtype='str', engine='c')

abc = pd.read_csv(all_predictions_ccre_els_bedpe, sep='\t', header=None, usecols=[0, 1, 2, 4, 5, 6, 7, 10, 11], names=col_names_abc, dtype='str', engine='c')
```


```python
intersection = {}
for k in range(nb_benchmarks):
    intersection[names[k]] = benchmarks[names[k]].copy().merge(abc.iloc[:,[5,6,8]], how='inner', left_on=col_names[6], right_on=col_names_abc[6])
    
col_summary = ['ABC.length', 'BENGI.length', 'intersection.length', 'intersection.length (%)', 'BENGI.positives', 'intersection.positives', 'intersection.positives (%)']
summary = pd.DataFrame(columns=col_summary)
for k in range(nb_benchmarks):
    len_abc = len(abc)
    len_bengi = len(benchmarks[names[k]])
    len_inter = len(intersection[names[k]])
    ratio = len_inter/len_bengi*100
    bengi_positives = sum(benchmarks[names[k]][col_names[7]] == "1")
    inter_positives = sum(intersection[names[k]][col_names[7]] == "1")
    ratio_positives = inter_positives/bengi_positives*100
    summary = summary.append(pd.DataFrame([[len_abc, len_bengi, len_inter, ratio, bengi_positives, inter_positives, ratio_positives]], columns=col_summary), ignore_index=True)
summary.index = names
```


```python
summary
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
    
    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ABC.length</th>
      <th>BENGI.length</th>
      <th>intersection.length</th>
      <th>intersection.length (%)</th>
      <th>BENGI.positives</th>
      <th>intersection.positives</th>
      <th>intersection.positives (%)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CHiC</th>
      <td>1084021</td>
      <td>375728</td>
      <td>70543</td>
      <td>18.775018</td>
      <td>88245</td>
      <td>16783</td>
      <td>19.018641</td>
    </tr>
    <tr>
      <th>CTCF</th>
      <td>1084021</td>
      <td>105016</td>
      <td>18348</td>
      <td>17.471623</td>
      <td>7591</td>
      <td>1311</td>
      <td>17.270452</td>
    </tr>
    <tr>
      <th>GEUVADIS</th>
      <td>1084021</td>
      <td>50999</td>
      <td>11471</td>
      <td>22.492598</td>
      <td>2073</td>
      <td>534</td>
      <td>25.759768</td>
    </tr>
    <tr>
      <th>GTEx</th>
      <td>1084021</td>
      <td>38200</td>
      <td>9770</td>
      <td>25.575916</td>
      <td>1301</td>
      <td>377</td>
      <td>28.977709</td>
    </tr>
    <tr>
      <th>HiC</th>
      <td>1084021</td>
      <td>153739</td>
      <td>31952</td>
      <td>20.783276</td>
      <td>3404</td>
      <td>707</td>
      <td>20.769683</td>
    </tr>
    <tr>
      <th>RNAPII</th>
      <td>1084021</td>
      <td>157235</td>
      <td>31289</td>
      <td>19.899513</td>
      <td>23699</td>
      <td>4981</td>
      <td>21.017764</td>
    </tr>
  </tbody>
</table>
</div>



#### Export files


```python
if not os.path.isdir(predictions+"BENGI_intersect/"):
    !mkdir $predictions"BENGI_intersect/"
for k in range(nb_benchmarks):
    intersection[names[k]].drop(intersection[names[1]].columns[-1], axis=1).to_csv(predictions+"BENGI_intersect/AllPredictions.intersect."+file_names[k]+".bedpe",sep='\t',header=False,index=False)
```

