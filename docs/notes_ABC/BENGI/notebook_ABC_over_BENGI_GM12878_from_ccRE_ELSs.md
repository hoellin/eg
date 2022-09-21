
# ABC model over GM12878 with same accessions as for BENGI benchmarks, but starting from ENCODE ccRE-dELS as candidate enhancers.

## How to use this notebook?
First, make a copy of [this notebook](http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/ABC/BENGI_GM12878_from_ccRE_ELSs.ipynb) on your computer / cluster.

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
cell_type = "gm12878" # used all along the notebook
specie = "homo_sapiens/"+version+"/" # used only in the current cell to compute some paths 

# Where to store run-specific references, scripts, intermediate files, predictions, etc
work_dir = "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/BENGI_GM12878_from_ccRE_ELSs/"

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
slurm_2 = scripts+"step2.sh"
slurm_3 = scripts+"step3.sh"
```


```python
ccREs_dir = results_dir+"ccREs/"+specie+"V1/"+cell_type+'/'

ccREs_accession = "ENCFF028SGJ"
TSS_list = "/work2/project/regenet/workspace/thoellinger/BENGI/Benchmark/Annotations/GENCODEv19-TSSs.bed"

ccREs = ccREs_dir+ccREs_accession+".bed" # need not exist yet
ccREs_ELS = ccREs_dir+cell_type+"_ccRE_ELSs.bed" # need not exist yet
ccREs_dELS = ccREs_dir+cell_type+"_ccRE_dELSs.bed" # need not exist yet

candidateRegions = reference_dir+"candidateRegions.bed"
```

## Data acquisition

### Chromatin accessibility (DNase-seq)


```python
if not os.path.isfile(dnase_file_rep1):
    !wget https://www.encodeproject.org/files/$dnase_rep1/@@download/$dnase_rep1$".bam" -P $dnase_dir
            
if not os.path.isfile(dnase_file_rep2):
    !wget https://www.encodeproject.org/files/$dnase_rep2/@@download/$dnase_rep2$".bam" -P $dnase_dir
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
    !wget https://www.encodeproject.org/files/$rnaseq_rep1/@@download/$rnaseq_rep1$".bam" -P $expression_dir
            
if not os.path.isfile(rnaseq_file_rep2):
    !wget https://www.encodeproject.org/files/$rnaseq_rep2/@@download/$rnaseq_rep2$".bam" -P $expression_dir
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

### ccREs
As candidate elements, we use ENCODE ccREs with distal enhancer-like signature (ccRE-dELS) for GM12878. Hence we  need to download ccREs with annotations for this cell type.


```python
if not os.path.isfile(ccREs+".gz"):
    !wget https://www.encodeproject.org/files/$ccREs_accession/@@download/$ccREs_accession$".bed.gz" -P $ccREs_dir
        
if not os.path.isfile(ccREs):
    !gzip -cd $ccREs".gz" > $ccREs
```

    --2021-01-07 12:19:58--  https://www.encodeproject.org/files/ENCFF028SGJ/@@download/ENCFF028SGJ.bed.gz
    Resolving www.encodeproject.org (www.encodeproject.org)... 34.211.244.144
    Connecting to www.encodeproject.org (www.encodeproject.org)|34.211.244.144|:443... connected.
    HTTP request sent, awaiting response... 307 Temporary Redirect
    Location: https://encode-public.s3.amazonaws.com/2017/08/21/f281e5ef-2e5d-4552-9716-a41bccbeea6a/ENCFF028SGJ.bed.gz?response-content-disposition=attachment%3B%20filename%3DENCFF028SGJ.bed.gz&AWSAccessKeyId=ASIATGZNGCNX23M77F6O&Signature=bBXj4VILUYZvYSml8FbwUAfdGVA%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEHsaCXVzLXdlc3QtMiJGMEQCIHPL2nO9fZlB01OMTFyBPq%2BAWvq9wuG%2B%2BWatPxVfheciAiBbFLXoWBz3E84MZvcSuuRNR2C5L2PbrSNg9n002iknOiq0AwhEEAAaDDIyMDc0ODcxNDg2MyIM0uGls%2BBaw6Z%2Fs4w9KpEDpEB18RAanxxcpOfloqMsXPBOAamsg32N6%2Borzy2nmXI9a3OlLccK35a4dzftO%2B4u7qQXw%2BmCOZgYBBnrjpZseQ6tX5XP3duXjuEmhOqY2VgMMR1M%2BmSJ5LILHNNMwiqj2i7Ld9oleRg9Mg9xPNzHGyJOWySoh8HJIVAB4nV1ZQ6V6Ky0uwu8rD185ahx%2BU%2BEUGD8mE8%2Bab8LD3KKfrmRh0InbFRnaNPv4ZhkF3dFD0LpMqqK9runudgZYT1JAHMZVKJ7x8PeoYpc6gGUFucaxcTUI0psdF9MP92%2FSAQob9VDnIV9mbOqk0wgqTOQIHTKs7b5VtLz%2FppBwBCfer1L1JvIbH4TA30ayx1J3L5jIRIbAiyX7p8X%2B1ekCnGYCi8ka57eljutwxv49ng3MgYuYU7gN%2B63dXUrZsmgxo7E1O9g9XO2Ha6vhBxKm5uwARO93kWKtbfdbpY6MtfmuhXHygE2sujmzPO5pg2MXfaAggW2gWjjOu5bEoj78YEiGbix%2FO1S6kimEpUY7AL4HrdqqM4wk9bb%2FwU67AH1DT3H7Kyftp1sKFWXQrlZj%2BtfeyKdtYSguln6VCf8jaUIyLq%2BNUmnqV7RstrbKDg7b3KWGTzyMWxO%2BNdv1CQxo8D7VmimQ58xganSXQXWT8Z7kLnPxzvnxhmyNCx3YdwFo2VPPiPINAZ%2Bh1EMqB6cAfGxsh03q3abmLGrWETQsdKul3NpFN6YhxH7Lp37Y0VsLHiwMc5%2BLzor83OmR5hfOG4NdbVElmyscq1mad5CuAcCFa%2F2eScFnMswZjxLW1DIFRGK219x33Brx2Ex8YQay8gnBTsvBXb%2BE04akm1txRiM5INk71wxSCgV0A%3D%3D&Expires=1610147999 [following]
    --2021-01-07 12:19:59--  https://encode-public.s3.amazonaws.com/2017/08/21/f281e5ef-2e5d-4552-9716-a41bccbeea6a/ENCFF028SGJ.bed.gz?response-content-disposition=attachment%3B%20filename%3DENCFF028SGJ.bed.gz&AWSAccessKeyId=ASIATGZNGCNX23M77F6O&Signature=bBXj4VILUYZvYSml8FbwUAfdGVA%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEHsaCXVzLXdlc3QtMiJGMEQCIHPL2nO9fZlB01OMTFyBPq%2BAWvq9wuG%2B%2BWatPxVfheciAiBbFLXoWBz3E84MZvcSuuRNR2C5L2PbrSNg9n002iknOiq0AwhEEAAaDDIyMDc0ODcxNDg2MyIM0uGls%2BBaw6Z%2Fs4w9KpEDpEB18RAanxxcpOfloqMsXPBOAamsg32N6%2Borzy2nmXI9a3OlLccK35a4dzftO%2B4u7qQXw%2BmCOZgYBBnrjpZseQ6tX5XP3duXjuEmhOqY2VgMMR1M%2BmSJ5LILHNNMwiqj2i7Ld9oleRg9Mg9xPNzHGyJOWySoh8HJIVAB4nV1ZQ6V6Ky0uwu8rD185ahx%2BU%2BEUGD8mE8%2Bab8LD3KKfrmRh0InbFRnaNPv4ZhkF3dFD0LpMqqK9runudgZYT1JAHMZVKJ7x8PeoYpc6gGUFucaxcTUI0psdF9MP92%2FSAQob9VDnIV9mbOqk0wgqTOQIHTKs7b5VtLz%2FppBwBCfer1L1JvIbH4TA30ayx1J3L5jIRIbAiyX7p8X%2B1ekCnGYCi8ka57eljutwxv49ng3MgYuYU7gN%2B63dXUrZsmgxo7E1O9g9XO2Ha6vhBxKm5uwARO93kWKtbfdbpY6MtfmuhXHygE2sujmzPO5pg2MXfaAggW2gWjjOu5bEoj78YEiGbix%2FO1S6kimEpUY7AL4HrdqqM4wk9bb%2FwU67AH1DT3H7Kyftp1sKFWXQrlZj%2BtfeyKdtYSguln6VCf8jaUIyLq%2BNUmnqV7RstrbKDg7b3KWGTzyMWxO%2BNdv1CQxo8D7VmimQ58xganSXQXWT8Z7kLnPxzvnxhmyNCx3YdwFo2VPPiPINAZ%2Bh1EMqB6cAfGxsh03q3abmLGrWETQsdKul3NpFN6YhxH7Lp37Y0VsLHiwMc5%2BLzor83OmR5hfOG4NdbVElmyscq1mad5CuAcCFa%2F2eScFnMswZjxLW1DIFRGK219x33Brx2Ex8YQay8gnBTsvBXb%2BE04akm1txRiM5INk71wxSCgV0A%3D%3D&Expires=1610147999
    Resolving encode-public.s3.amazonaws.com (encode-public.s3.amazonaws.com)... 52.218.241.178
    Connecting to encode-public.s3.amazonaws.com (encode-public.s3.amazonaws.com)|52.218.241.178|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 17995639 (17M) [binary/octet-stream]
    Saving to: ‘/work2/project/regenet/results/ccREs/homo_sapiens/hg19/V1/gm12878/ENCFF028SGJ.bed.gz’
    
    100%[======================================>] 17 995 639  3,35MB/s   in 15s    
    
    2021-01-07 12:20:14 (1,17 MB/s) - ‘/work2/project/regenet/results/ccREs/homo_sapiens/hg19/V1/gm12878/ENCFF028SGJ.bed.gz’ saved [17995639/17995639]



## Data reprocessing

### Create expression table

We use `compute_mean_expression.awk` to take as the reference TPM value, the mean value for the 2 replicates.


```bash
%%bash -s "$compute_mean_expression" "$reference_dir" "$rnaseq_file_rep1" "$rnaseq_file_rep2" "$cell_type" "$rnaseq_rep1" "$rnaseq_rep2"
if [[ ! -d $2"expression/" ]]; then mkdir $2"expression/"; fi
awk -f $1 $3 $4 > $2"expression/"$5"."$6"_"$7".mean.TPM.txt"
```

### Filter gene annotation


```bash
%%bash -s "$work_dir" "$genome_file" "$annotations_dir" "$gene_annotation" "$light_annotation"
if [[ ! -d $3 ]]; then mkdir $3; fi
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && module load bioinfo/bedtools-2.27.1
awk '{if(NR==FNR){!genome[$1]++; next}; if(genome[$1]){if($3 ~ /(^gene$)/){print $1"\t"$4"\t"$5"\t"substr($10,2,length($10)-3)"\t"0"\t"$7}}}' $2 $4 | bedtools sort -g $2 -i stdin > $5
```

## Running the ABC model

### Step 1: define candidate elements

As candidate elements, we use ENCODE ccREs with distal enhancer-like signature (ccRE-dELS) for GM12878. First we extract only ccRE-ELS.


```bash
%%bash -s "$ccREs" "$ccREs_ELS"
if [[ ! -f $2 ]]; then
    awk -v OFS='\t' '$9 ~ /(^255,205,0)/ {print $1, $2, $3, ".", $4, "Enhancer-like"}' $1 > $2
fi
```


```bash
%%bash -s "$ccREs_ELS" "$genome_file"
if [[ ! -f $1".sorted" ]]
then
    export V1=$1 V2=$2
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun bash
    bedtools sort -faidx $V2 -i $V1 > $V1".sorted"
fi
```

    srun: job 22240197 queued and waiting for resources
    srun: job 22240197 has been allocated resources


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

    srun: job 22240200 queued and waiting for resources
    srun: job 22240200 has been allocated resources



```bash
%%bash -s "$ccREs_dELS" "$candidateRegions"
if [[ ! -f $2 ]]; then
    awk -v OFS='\t' '{print $1, $2, $3}' $1 > $2
fi
```

### Step 2: quantifying enhancer activity

#### Use ABC `run.neighborhoods.py`


```python
%%writetemplate $slurm_2
#!/bin/sh
python {run_neighborhoods} \
--candidate_enhancer_regions {candidateRegions} \
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
sbatch --mem=16G --cpus-per-task=1 -J step2 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 22240256



```python
!seff 22240256 #should take 2 hours => at 14h30 it will be ok
```

    Job ID: 22240256
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 02:48:52
    CPU Efficiency: 99.85% of 02:49:07 core-walltime
    Job Wall-clock time: 02:49:07
    Memory Utilized: 11.83 GB
    Memory Efficiency: 36.96% of 32.00 GB


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


```bash
%%bash -s "$work_dir" "$slurm_3" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=8G --cpus-per-task=1 -J step3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 22242814



```python
!seff 22242814
```

    Job ID: 22242814
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:07:49
    CPU Efficiency: 100.00% of 00:07:49 core-walltime
    Job Wall-clock time: 00:07:49
    Memory Utilized: 3.47 GB
    Memory Efficiency: 43.36% of 8.00 GB


## Evaluation of ABC predictions on BENGI datasets

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
if [[ ! -f $1 ]]; then
    awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' $2 > $1;
    tail -n+2 $3 | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' >> $1;
fi
```

The following is very RAM-demanding, we recommend using at least 64GB of ram to avoid "out of memory" errors.


```bash
%%bash -s "$all_predictions" "$all_predictions_sorted" "$genome_file"
export V1=$1 V2=$2 V3=$3
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/bedtools-2.27.1
srun --mem=64G bash
if [[ ! -f $V2 ]]; then
    tail -n+2 $V1 | bedtools sort -faidx $V3 -i stdin  > $V2
fi
```

    srun: job 22242943 queued and waiting for resources
    srun: job 22242943 has been allocated resources


Intersect ABC candidate regions with ccRE-dELS.


```bash
%%bash -s "$ccREs_dELS" "$all_predictions_sorted" "$all_predictions_ccre_els"
export V1=$1 V2=$2 V3=$3
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/bedtools-2.27.1
srun --mem=32G bash
if [[ ! -f $V3 ]]; then
    bedtools intersect -sorted -wo -a $V1 -b $V2 > $V3
fi
```

    srun: job 22242951 queued and waiting for resources
    srun: job 22242951 has been allocated resources


When multiple ABC candidate enhancers for a given gene overlap the same ccRE-ELS, we keep only the one that maximize the ABC score (the "powerlaw.Score" if ABC ran without 3D data). Here it did not matter as there were no duplicates ; although in general we may not want to keep only the max, if at the end of the day it results in a smaller intersection between ABC predictions and BENGI datasets.


```bash
%%bash -s "$light_annotation" "$all_predictions_ccre_els" "$all_predictions_ccre_els_bedpe"
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
```

    srun: job 22242953 queued and waiting for resources
    srun: job 22242953 has been allocated resources


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
      <td>5879612</td>
      <td>375728</td>
      <td>375728</td>
      <td>100.0</td>
      <td>88245</td>
      <td>88245</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>CTCF</th>
      <td>5879612</td>
      <td>105016</td>
      <td>105016</td>
      <td>100.0</td>
      <td>7591</td>
      <td>7591</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>GEUVADIS</th>
      <td>5879612</td>
      <td>50999</td>
      <td>50999</td>
      <td>100.0</td>
      <td>2073</td>
      <td>2073</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>GTEx</th>
      <td>5879612</td>
      <td>38200</td>
      <td>38200</td>
      <td>100.0</td>
      <td>1301</td>
      <td>1301</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>HiC</th>
      <td>5879612</td>
      <td>153739</td>
      <td>153739</td>
      <td>100.0</td>
      <td>3404</td>
      <td>3404</td>
      <td>100.0</td>
    </tr>
    <tr>
      <th>RNAPII</th>
      <td>5879612</td>
      <td>157235</td>
      <td>157235</td>
      <td>100.0</td>
      <td>23699</td>
      <td>23699</td>
      <td>100.0</td>
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
