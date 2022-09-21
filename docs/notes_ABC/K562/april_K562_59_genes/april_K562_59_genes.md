
# ABC over K562

## How to use this notebook?
First, make a copy of [this notebook](http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/ABC/april_K562_59_genes.ipynb) on your computer / cluster.

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
work_dir = "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/april_K562_candidates_in_whitelist_59_genes/"
gene_annotation_dir = "/work2/project/regenet/workspace/thoellinger/RefSeq/"

scripts = work_dir+"hi_slurm/" # where to store scripts specific to this run
results_dir = "/work2/project/regenet/results/" # used to compute paths where data for this cell type is
                                                # stored / will be downloaded
dnase_dir = results_dir+"dnaseseq/"+specie+cell_type+'/'
chipseq_dir = results_dir+"chipseq/h3k27ac/"+specie+cell_type+'/'
expression_dir = results_dir+"rnaseq/"+specie+cell_type+'/'
blacklist_dir = results_dir+"multi/black.lists/"+specie # where the ENCODE blacklist is (going to be) stored

blacklist = blacklist_dir+version+"-blacklist.v2.bed" # need not exist yet

# Gene annotation in the gtf format and gene name / gene id table (1st column id, 2nd column name)
gene_annotation_link = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz"

gene_annotation = gene_annotation_dir+"RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed"
TSS500bp_annotation = gene_annotation_dir+"RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_TSS_that_maximize_nb_of_distinct_coding_transcripts_with_3_more_genes_added_by_hand.TSS500bp.bed"
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
rnaseq = "ENCFF934YBO" # exp ENCSR000AEM, indirectly

# Where to download the blacklist. Will not be used unless blacklist not found yet.
blacklist_link = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"

dnase_file_rep1 = dnase_dir+dnase_rep1+".bam" # you shall not change this
dnase_file_rep2 = dnase_dir+dnase_rep2+".bam"  # you shall not change this
h3k27_file_rep1 = chipseq_dir+h3k27_rep1+".bam"  # you shall not change this
h3k27_file_rep2 = chipseq_dir+h3k27_rep2+".bam"  # you shall not change this
rnaseq_file = expression_dir+rnaseq+".tsv"  # you shall not change thisUbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.txt

reference_dir = work_dir+"reference/" # you shall not change this
annotations_dir = reference_dir+"gene_annotation/" # you shall not change this
peaks = work_dir+"ABC_output/Peaks/" # you shall not change this
neighborhoods = work_dir+"ABC_output/Neighborhoods/" # you shall not change this
predictions = work_dir+"ABC_output/Predictions/" # you shall not change this

ubiquitously_expressed_genes = work_dir+"../reference/UbiquitouslyExpressedGenesHG19_with_gene_details_from_table_S5b_intersected_with_our_RefSeq_filtered_annotation_with_3_more_genes.txt"
whitelist = reference_dir+"whitelist.bed"

make_candidate_regions = work_dir+"../src/makeCandidateRegions.py"
run_neighborhoods = work_dir+"../src/run.neighborhoods.py"
predict = work_dir+"../src/predict.py"

genome_file = reference_dir+"chr_sizes" # although we generate this file later on in the notebook,
# one can directly obtain it from
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
# after removing scaffolds different from chromosomes and sorting it

gene_expression_table = reference_dir+"expression/"+cell_type+"."+rnaseq+".TPM.txt"

CRiFF_dir = "/work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/"
all_criff = CRiFF_dir+"3863.fulco.bedpe.sorted.geneNames"
candidateRegions = reference_dir+"candidateRegions.bed"
```


```python
if not os.path.isdir(reference_dir):
    !mkdir -p $reference_dir
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

## Data acquisition

### Gene annotation


```python
if not (os.path.isfile(gene_annotation_dir+"hg19.refGene.gtf.gz") and os.path.isfile(gene_annotation_dir+"hg19.refGene.gtf")):
    !wget $gene_annotation_link -P $gene_annotation_dir
if not os.path.isfile(gene_annotation_dir+"hg19.refGene.gtf"):
    !gzip -cd $gene_annotation_dir$"hg19.refGene.gtf.gz" > $gene_annotation_dir$"hg19.refGene.gtf"
```

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
if not os.path.isfile(rnaseq_file):
    !wget https://www.encodeproject.org/files/$rnaseq/@@download/$rnaseq".tsv" -P $expression_dir
```

### Genome file


```bash
%%bash -s "$work_dir" "$dnase_file_rep1" "$genome_file"
if [[ ! -f $3 ]]; then
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && module load bioinfo/samtools-1.9
    samtools view -H $2 |grep SQ |cut -f 2,3 |cut -c 4- |awk 'BEGIN{FS="\t"} {split($2,locus,":"); if($1 ~ /^(chr)([0-9]{1,2}$)|(M$)|(X$)|(Y$)/){print($1"\t"locus[2])}}' > $3
    if [[ ! -f $3".bed" ]]; then
        awk '{print $1"\t"0"\t"$2}' $3 > $3".bed"
    fi
fi
```

## Data reprocessing

### Create expression table


```bash
%%bash -s "$gnid_gname" "$reference_dir" "$rnaseq_file" "$cell_type" "$rnaseq"
if [[ ! -d $2"expression/" ]]; then mkdir $2"expression/"; fi
awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){name[$1]=$2; next} if(($1 ~ /(^ENS)/) && name[$1]){TPM[name[$1]] += $6; count[name[$1]]++}} END{for(u in TPM){print u, TPM[u]/count[u]}}' $1 $3 > $2"expression/"$4"."$5".TPM.txt"
```

### Filter gene annotation

We have done the following step "by hand".


```python
!wc -l $gene_annotation
```

    19326 /work2/project/regenet/workspace/thoellinger/RefSeq/RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed



```python
!wc -l $TSS500bp_annotation
```

    19326 /work2/project/regenet/workspace/thoellinger/RefSeq/RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_TSS_that_maximize_nb_of_distinct_coding_transcripts_with_3_more_genes_added_by_hand.TSS500bp.bed



```python
whitelist # shall not exist yet, we create during step 1
```




    '/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/april_K562_candidates_in_whitelist_59_genes/reference/whitelist.bed'



## Running the ABC model

### Step 1: define candidate elements

#### Call peaks with `macs2`

Found in Fulco et al:
> For K562, we concatenated all peaks called by ENCODE in both replicate DNase-seq experiments (Supplementary Table 4). Given that the ENCODE peaks were initially 150 bp in length, we extended each of these peaks by 175 bp to arrive at candidate elements that were 500 bp in length. We then removed any peaks overlapping regions of the genome that have been observed to accumulate anomalous number of reads in epigenetic sequencing experiments (blacklisted regions42,43, downloaded from https://sites.google.com/site/anshulkundaje/projects/blacklists). To this peak list we added 500-bp regions centered on the transcription start site of all genes. Any overlapping regions resulting from these additions or extensions were merged. In total, this procedure resulted in 162,181 candidate regions in K562, whose average length was 576 bp (Extended Data Fig. 2b).


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

    Submitted batch job 24960364



```python
!seff 24960364
```

    Job ID: 24960364
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:17:25
    CPU Efficiency: 98.58% of 00:17:40 core-walltime
    Job Wall-clock time: 00:17:40
    Memory Utilized: 746.69 MB
    Memory Efficiency: 18.23% of 4.00 GB


#### Create the whitelist


We remind that for this run, we want to investigate what happens when we forcibly include as candidate regions not only the TSS-centered 500bp regions curated above, but also all the candidate regions found in the CRISPRi-FlowFISH validation dataset.


```bash
%%bash -s "$TSS500bp_annotation" "$all_criff" "$genome_file" "$whitelist"
if [[ ! -f $V4 ]]
then
    export V1=$1 V2=$2 V3=$3 V4=$4
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun --mem=16G bash
    awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, "0"}' $V1 $V2 |bedtools sort -faidx $V3 -i > $V4
fi
```

    srun: job 25045474 queued and waiting for resources
    srun: job 25045474 has been allocated resources



```python
!wc -l $whitelist
```

    23189 /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/april_K562_candidates_in_whitelist_59_genes/reference/whitelist.bed


#### Use ABC `makeCandidateRegions.py` to define candidate regions

As there might be peaks over many kinds of scaffolds listed in the `.narrowPeak` and `.xls` output files, we modify them keep only chromosomes only (otherwise `bedtools sort` below would not work). **But we have to keep in mind that the `.r` will still corresponds to DNase peaks over all scaffolds**.

Note that in the current run, all peaks were already over chromosomes only, hence the following step makes no differences in the current run.


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

    srun: job 24961003 queued and waiting for resources
    srun: job 24961003 has been allocated resources


Now as `work_dir ../src/makeCandidateRegions.py` does not support multiple `bam` inputs, we concatenate the bam of the 2 replicates.  Indeed, `makeCandidateRegions.py` will use the bam input to count DNase reads on each peak of the`.macs2_peaks.narrowPeak.sorted` output found with `macs2`.

(WARNING: depending on the size of the input bam, which is much variable, the execution on 1 thread only can take from a few seconds to several hours, so better use at least 8 threads)


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

    Submitted batch job 24961055



```python
!seff 24961055
```

    Job ID: 24961055
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Nodes: 1
    Cores per node: 8
    CPU Utilized: 00:00:00
    CPU Efficiency: 0.00% of 00:00:08 core-walltime
    Job Wall-clock time: 00:00:01
    Memory Utilized: 1.39 MB
    Memory Efficiency: 0.07% of 2.00 GB


Now we can launch `makeCandidateRegions` on slurm (note that as we work with K562 cell line, we take `nStrongestPeaks` to be 175,000 - which is the default argument, since Fulco et al restricted to 150,000 only for other cell lines to approximately match the whole number of DNase peaks of K562).

> To approximately match the number of candidate elements considered in K562, we then counted DNase-seq (or ATAC-seq) reads overlapping these peaks and kept the 150,000 with the highest number of read counts. 


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
--regions_whitelist {whitelist} \
--peakExtendFromSummit 250 \
--nStrongestPeaks 175000
#Expected output: params.txt, foo1_foo2.macs2_peaks.narrowPeak.sorted.candidateRegions.bed, foo1_foo2.macs2_peaks.narrowPeak.sorted.foo1_foo2.bam.Counts.bed
```

WARNING: depending on the length of the inputs, this script may be very RAM-demanding. So we recommend to allocate at least 16 GB RAM to avoid "out of memory" errors.


```bash
%%bash -s "$work_dir" "$slurm_1_3" "$mail_user"
CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
conda activate base && conda activate abcmodel
module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
sbatch --mem=32G --cpus-per-task=1 -J step1.3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 25045626



```python
!seff 25045626
```

    Job ID: 25045626
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:22:13
    CPU Efficiency: 99.55% of 00:22:19 core-walltime
    Job Wall-clock time: 00:22:19
    Memory Utilized: 18.79 GB
    Memory Efficiency: 58.71% of 32.00 GB


Not that it did not result in the very same number of candidate regions as in Fulco et al 's paper (163,654 instead of 162,181).

### Step 2: quantifying enhancer activity

#### Use ABC `run.neighborhoods.py`


```python
%%writetemplate $slurm_2
#!/bin/sh
python {run_neighborhoods} \
--candidate_enhancer_regions {peaks}{dnase_rep1}_{dnase_rep2}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes {gene_annotation} \
--H3K27ac {h3k27_file_rep1} \
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

    Submitted batch job 25050270



```python
!seff 25050270
```

    Job ID: 25050270
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:16:51
    CPU Efficiency: 99.51% of 00:16:56 core-walltime
    Job Wall-clock time: 00:16:56
    Memory Utilized: 2.72 GB
    Memory Efficiency: 8.51% of 32.00 GB


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
sbatch --mem=16G --cpus-per-task=1 -J step3 --mail-user=$3 --mail-type=END,FAIL --workdir=$1 --export=ALL -p workq $2
```

    Submitted batch job 25055299



```python
!seff 25055299
```

    Job ID: 25055299
    Cluster: genobull
    User/Group: thoellinger/U1220
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:12:56
    CPU Efficiency: 97.98% of 00:13:12 core-walltime
    Job Wall-clock time: 00:13:12
    Memory Utilized: 7.36 GB
    Memory Efficiency: 46.00% of 16.00 GB


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

Keep only the columns of interests.


```bash
%%bash -s "$all_predictions" "$predictions_expressed"
if [[ ! -f $1 ]]
then
    awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' $2 > $1;
fi
# Or merge ABC predictions for expressed and non expressed genes, keep only columns of interest and sort the result:
#%%bash -s "$all_predictions" "$predictions_expressed" "$predictions_non_expressed"
#if [[ ! -f $1 ]]
#then
#    awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' $2 > $1;
#    tail -n+2 $3 | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$10"\t"$10"\t"$9"\t"$20"\t.\t.\t"$13"\t"$11}' >> $1;
#fi
```

Now sort the result. The following is very RAM-demanding, we recommend using at least 64GB of ram to avoid "out of memory" errors.


```bash
%%bash -s "$all_predictions" "$all_predictions_sorted" "$genome_file"
if [[ ! -f $V2 ]]
then
    export V1=$1 V2=$2 V3=$3
    CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base && conda activate abcmodel
    module load bioinfo/bedtools-2.27.1
    srun --mem=32G bash
    tail -n+2 $V1 | bedtools sort -faidx $V3 -i stdin  > $V2
fi
```

    srun: job 25055720 queued and waiting for resources
    srun: job 25055720 has been allocated resources


### Intersection with Fulco et al 's CRISPRi-FlowFISH validation dataset

Remark: this section is not strictly speaking included in the notebook ; I performed the computations "by hand".


```python
gene_annotation
```




    '/work2/project/regenet/workspace/thoellinger/RefSeq/RefSeq_GRCh37_p13_coding_genes_g300bp_s2Mbp_with_3_more_genes_added_by_hand.sorted.bed'



First we intersect (regulatory regions of) all predictions (for expressed genes) with the CRISPRi-FlowFISH validation dataset using `bedtools intersect`.

> ```bash
> conda activate base && module load bioinfo/bedtools-2.27.1 && srun --mem=64G --pty bash
> bedtools intersect -sorted -wo -a ABC_output/Predictions/AllPredictions.bedpe.sorted -b /work2/project/regenet/workspace/thoellinger/CRISPRi_FlowFISH/k562/3863.fulco.bedpe.sorted.geneNames.byhand -g reference/chr_sizes > ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe
> ```

Now we keep only the lines for which the gene names match:
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $0}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe |wc -l
> ```
> 3860

The ABC model did not perform predictions for at least 3 pairs over the 3863 unique pairs of the reference validation dataset... Well that's not much this time. We keep the columns of interest only, and for identical enhancer-gene pairs, we choose to keep the one that maximizes the ABC score:

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $23, $13, $14, $15, $24, $8, $25}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {uniq=$1"\t"$2"\t"$3"\t"$4"\t"$5; if($6>val[uniq]){val[uniq]=$6}} END{for(u in val){print u, val[u]}}' |wc -l
> ```
> 3860

Luckily there was no degenerescence: we don't have predictions for exactly 3 pairs among the 3863 pairs of the validation dataset. Let's continue to further analyse the results with R.

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if($7==$24){print $23, $13, $14, $15, $24, $8, $25}}' ABC_output/Predictions/enhancer_predictions_intersected_with_CRISPRi_FlowFISH.bedpe |awk 'BEGIN{FS="\t"; OFS="\t"} {uniq=$1"\t"$2"\t"$3"\t"$4"\t"$5; if($6>val[uniq]){val[uniq]=$6}} END{for(u in val){print u, val[u]}}' > ABC_output/Predictions/predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.bedpe
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {print $2, $3, $4, $5, $6, $1}' ABC_output/Predictions/predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.bedpe |bedtools sort -faidx reference/chr_sizes.bed -i |awk 'BEGIN{FS="\t"; OFS="\t"} {print $6, $1, $2, $3, $4, $5}' > ABC_output/Predictions/predictions_intersected_with_CRISPRi_FlowFISH_regions_that_maximize_the_ABC_score.bedpe.sorted
```

### Analysis with R

Over the 3860 predictions remaining (out of 3863 contained in the validation dataset), 106 positives (among 109) remain.
