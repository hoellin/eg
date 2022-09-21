
# ABC over K562

## How to use this notebook?
First, make a copy of [this notebook](http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/ABC/ABC_generic.ipynb) on your computer / cluster.

Then, to use this notebook, one should only have to carefully modify the content of the "Set working directory" section, then to execute the notebook cell by cell, in the correct order. After execution of each cell, remember to check for errors before executing the next one.

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
cell_type = "k562" # used all along the notebook
specie = "homo_sapiens/"+version+"/" # used only in the current cell to compute some paths 

# Where to store run-specific references, scripts, intermediate files, predictions, etc
work_dir = "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/new_K562/"

scripts = work_dir+"hi_slurm/" # where to store scripts specific to this run
results_dir = "/work2/project/regenet/results/" # used to compute paths where data for this cell type is
                                                # stored / will be downloaded
dnase_dir = results_dir+"dnaseseq/"+specie+cell_type+'/'
chipseq_dir = results_dir+"chipseq/h3k27ac/"+specie+cell_type+'/'
expression_dir = results_dir+"rnaseq/"+specie+cell_type+'/'
blacklist_dir = results_dir+"multi/"+specie # where the ENCODE blacklist is (going to be) stored

blacklist = blacklist_dir+version+"-blacklist.v2.bed" # need not exist yet

# Gene annotation in the gtf format and gene name / gene id table (1st column id, 2nd column name)
gene_annotation = "/work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf"
gnid_gname = "/work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gnid.gnname.tsv"

# Accessions for the current run. Just put accession numbers here, no need to download these accessions "by hand".
dnase_rep1 = "ENCFF001DOX" # see https://www.encodeproject.org/files/ENCFF001DOX/
                           # => "Original file name hg19/wgEncodeUwDnase/wgEncodeUwDnaseK562AlnRep1.bam"
                           # (actually wgEncodeUwDnaseK562AlnRep1.bam is the name fiven in Supplementary Table 4)
dnase_rep2 = "wgEncodeUwDnaseK562AlnRep2" # well we found new accession corresponding to wgEncodeUwDnaseK562AlnRep1.bam,
                # ENCFF001DOX (archived), but not to wgEncodeUwDnaseK562AlnRep2.bam.
                # Fortunately we found the file here: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwDnase/
h3k27_rep1 = "ENCFF384ZZM" # exp ENCSR000AKP
h3k27_rep2 = "ENCFF070PWH" # exp ENCSR000AKP
rnaseq_rep1 = "ENCFF172GIN" # exp ENCSR000CPH
rnaseq_rep2 = "ENCFF768TKT" # exp ENCSR000CPH

# Where to download the blacklist. Will not be used unless blacklist not found yet.
blacklist_link = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"

dnase_file_rep1 = dnase_dir+dnase_rep1+".bam" # you shall not change this
dnase_file_rep2 = dnase_dir+dnase_rep2+".bam"  # you shall not change this
h3k27_file_rep1 = chipseq_dir+h3k27_rep1+".bam"  # you shall not change this
h3k27_file_rep2 = chipseq_dir+h3k27_rep2+".bam"  # you shall not change this
rnaseq_file_rep1 = expression_dir+rnaseq_rep1+".tsv"  # you shall not change this
rnaseq_file_rep2 = expression_dir+rnaseq_rep2+".tsv"  # you shall not change this

reference_dir = work_dir+"reference/" # you shall not change this
annotations_dir = reference_dir+"gene_annotation/" # you shall not change this
peaks = work_dir+"ABC_output/Peaks/" # you shall not change this
neighborhoods = work_dir+"ABC_output/Neighborhoods/" # you shall not change this
predictions = work_dir+"ABC_output/Predictions/" # you shall not change this

light_annotation = annotations_dir+"gene_ids.bed" # you shall not change this

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
CRiFF_dir = "/work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/"
all_criff = CRiFF_dir+"3863.fulco.bedpe.sorted"

candidateRegions = reference_dir+"candidateRegions.bed"
```

## Data acquisition

### Chromatin accessibility (DNase-seq)


```python
if not os.path.isfile(dnase_file_rep1):
    !wget https://www.encodeproject.org/files/$dnase_rep1/@@download/$dnase_rep1$".bam" -P $dnase_dir
            
if not os.path.isfile(dnase_file_rep2):
    !wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwDnase/$dnase_rep2$".bam" -P $dnase_dir
```

### Blacklist


```python
if not os.path.isfile(blacklist):
    !wget $blacklist_link -P $blacklist_dir
```

### Histone mark H3K27ac ChIP-seq


```python
if not os.path.isfile(h3k27_file_rep1):
    !wget https://www.encodeproject.org/files/$h3k27_rep1/@@download/$h3k27_rep1$".bam" -P $chipseq_dir
            
if not os.path.isfile(h3k27_file_rep2):
    !wget https://www.encodeproject.org/files/$h3k27_rep2/@@download/$h3k27_rep2$".bam" -P $chipseq_dir
```

### Gene expression (polyA+ RNA-seq)


```python
if not os.path.isfile(rnaseq_file_rep1):
    !wget https://www.encodeproject.org/files/$rnaseq_rep1/@@download/$rnaseq_file_rep1 -P $expression_dir
            
if not os.path.isfile(rnaseq_file_rep2):
    !wget https://www.encodeproject.org/files/$rnaseq_rep2/@@download/$rnaseq_file_rep2 -P $expression_dir
```

### Ubiquitously expressed genes


```bash
%%bash -s "$gnid_gname" "$ubiquitous_gene_names" "$ubiquitously_expressed_genes"
if [[ -f $2 ]] && [[ ! -f $3 ]]; then
    awk 'BEGIN{sep="\t"} NR==FNR {id[$2]=$1; next} id[$1] {print id[$1]}' $1 $2 > $3
fi
```

### Genome file


```bash
%%bash -s "$work_dir" "$dnase_file_rep1" "$genome_file"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && module load bioinfo/samtools-1.9
samtools view -H $2 | grep SQ | cut -f 2,3 |cut -c 4- |awk 'BEGIN{FS="\t"} {split($2,locus,":"); if($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/){print($1"\t"locus[2])}}' > $3
awk '{print $1"\t"0"\t"$2}' $3 > $3".bed"
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

We keep genes only, and `gene_ids.bed` must fit following format:
> ```bash
> chrX	99883667	99894988	ENSG00000000003.10	0	-
> chrX	99839799	99854882	ENSG00000000005.5	0	+
> ```


```bash
%%bash -s "$work_dir" "$genome_file" "$annotations_dir" "$gene_annotation" "$light_annotation"
if [[ ! -d $3 ]]; then mkdir $3; fi
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && module load bioinfo/bedtools-2.27.1
awk '{if(NR==FNR){!genome[$1]++; next}; if(genome[$1]){if($3 ~ /(^gene$)/){print $1"\t"$4"\t"$5"\t"substr($10,2,length($10)-3)"\t"0"\t"$7}}}' $2 $4 | bedtools sort -g $2 -i stdin > $5
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

    Submitted batch job 23346597



```python
!seff 23346597
```

    Job ID: 23346597
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:16:55
    CPU Efficiency: 98.16% of 00:17:14 core-walltime
    Job Wall-clock time: 00:17:14
    Memory Utilized: 759.34 MB
    Memory Efficiency: 18.54% of 4.00 GB


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

    srun: job 23347039 queued and waiting for resources
    srun: job 23347039 has been allocated resources


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

    Submitted batch job 23347154



```python
!seff 23347154
```

    Job ID: 23347154
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Nodes: 1
    Cores per node: 8
    CPU Utilized: 00:11:08
    CPU Efficiency: 62.31% of 00:17:52 core-walltime
    Job Wall-clock time: 00:02:14
    Memory Utilized: 8.80 MB
    Memory Efficiency: 0.43% of 2.00 GB


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
sbatch --mem=64G --cpus-per-task=1 -J step1.3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 23358051



```python
!seff 23358051
```

    Job ID: 23358051
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:13:29
    CPU Efficiency: 99.63% of 00:13:32 core-walltime
    Job Wall-clock time: 00:13:32
    Memory Utilized: 18.79 GB
    Memory Efficiency: 29.36% of 64.00 GB


Not that it did not result in the same number of candidate regions as in Fulco et al 's paper (not even close - 136,205 instead of 162,181)... Actually we expected this as we did not use the same RNA-seq experiment.

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

WARNING: depending on the length of the inputs, this script may be very RAM-demanding. So we recommend to allocate at least 32 GB RAM to avoid "out of memory" errors.


```bash
%%bash -s "$work_dir" "$slurm_2" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=32G --cpus-per-task=1 -J step2 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 23359930



```python
!seff 23359930
```

    Job ID: 23359930
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:19:07
    CPU Efficiency: 99.83% of 00:19:09 core-walltime
    Job Wall-clock time: 00:19:09
    Memory Utilized: 2.76 GB
    Memory Efficiency: 8.62% of 32.00 GB


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

Again, that step is very RAM-demanding so we recommend to allocate at least 32G of RAM.


```bash
%%bash -s "$work_dir" "$slurm_3" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=64G --cpus-per-task=1 -J step3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 23361936



```python
!seff 23361936
```

    Job ID: 23361936
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:24:13
    CPU Efficiency: 99.73% of 00:24:17 core-walltime
    Job Wall-clock time: 00:24:17
    Memory Utilized: 16.53 GB
    Memory Efficiency: 25.82% of 64.00 GB


## Results

### Preprocessing


```python
predictions_non_expressed = predictions+"EnhancerPredictionsAllPutativeNonExpressedGenes.txt"
predictions_expressed = predictions+"EnhancerPredictionsAllPutative.txt"
all_predictions = predictions+"AllPredictions.bedpe" # does not exist yet
all_predictions_sorted = all_predictions+".sorted"
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


```python
predictions_expressed
```




    '/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/new_K562/ABC_output/Predictions/EnhancerPredictionsAllPutative.txt'




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

    srun: job 23362791 queued and waiting for resources
    srun: job 23362791 has been allocated resources


### Intersection with Fulco et al 's CRISPRi-FlowFISH validation dataset

Remark: this section is not strictly speaking included in the notebook ; I performed the computations "by hand".

```bash
conda activate base && module load bioinfo/bedtools-2.27.1 && srun --mem=64G --pty bash
```

> ```bash
> wc -l /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted
> ```
>
> 3863
>
> ```bash
> bedtools intersect -sorted -u -a /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted -b ABC_output/Predictions/AllPredictions.bedpe.sorted -g reference/chr_sizes |wc -l
> ```
>
> 3640

Not perfect, but enough to continue.

```bash
bedtools intersect -sorted -wo -a ABC_output/Predictions/AllPredictions.bedpe.sorted -b /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted -g reference/chr_sizes > ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe
```

> ```bash
> wc -l ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe
> ```
>
> 1597054

When we keep only the lines for which the gene ids are the same, there are still > 3863 predictions. Actually, multiple candidates sometimes overlap the same validation element, resulting in multiple ABC scores for a single validation data.

How to choose which ones to keep?

#### First strategy: we keep the ones that maximize the ABC score

> ```bash
> awk 'BEGIN{FS="\t"} {if($7==$24){print $7, $24}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe |wc -l
> ```
>
> 3768

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $23, $13, $14, $15, $24, $8, $25}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {uniq=$1"\t"$2"\t"$3"\t"$4"\t"$5; if($6>val[uniq]){val[uniq]=$6}} END{for(u in val){print u, val[u]}}' > ABC_output/Predictions/Predictions_over_CRISPRi_FlowFISH.bedpe
```

> ```bash
> wc -l ABC_output/Predictions/Predictions_over_CRISPRi_FlowFISH.bedpe
> ```
>
> 3640

#### Second strategy: we keep the ones that maximize the overlap

#### Third strategy: weighted average of the ABC scores, where weights are the sizes of the overlaps

#### 4th strategy: weighted average of the ABC scores, where weights are the sizes of the overlaps, divided by the sizes of the predicted regions

This one might be a little more complicated to implement.

### Analysis with R

Over the 3640 predictions remaining with 1st strategy (out of 3863 contained in the validation dataset), 105 positives (among 109) remain. There are 128 `NA` values in the ABC scores, but none of them pertain to the 105 predictions for ground positives.

### Analysis with R (results)

![Image: Precision-Recall curve for the first strategy](max_abc_precision_recall_positives_105_negatives_3535_with_128_NA_negatives.png)

