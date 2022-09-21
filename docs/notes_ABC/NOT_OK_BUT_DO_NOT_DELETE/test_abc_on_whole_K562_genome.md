# Test ABC model on the whole K562 genome

Adapted from [ABC model's README on github](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

The Activity-by-contact model predicts which enhancers regulate which genes on a cell type specific basis.

## Requirements

For each cell-type, the inputs to the ABC model are:

 * Required Inputs
     * bam file for Dnase-Seq or ATAC-Seq (indexed and sorted)

       The DNase-Seq file will contain the alignment of the genome-wide sequencing of regions sensitive to cleavage by DNase I of K562 leukemia cell, over the GRCh38 human's genome assembly.

     * bam file for H3K27ac ChIP-Seq (indexed and sorted)

       (Wiki) H3K27ac is an epigenetic modification to the DNA packaging protein Histone H3 (one of the five main histones involved in the structure of chromatin in eukaryotic cells). It is a mark that indicates the acetylation at the 27th lysine residue of the histone H3 protein. H3K27ac is associated with the higher activation of transcription and therefore defined as an *active* enhancer mark. H3K27ac is found at both proximal and distal regions of transcription start site (TSS).
 * Optional Inputs
     * Hi-C data

       These Hi-C data serve to compute the contact frequency between regions. If not available, they can be estimated as it has been proven that the contact frequency approximately follows a power law of the distance.

     * A measure of gene expression

In addition the following (non-cell-type specific) genome annotation files are required

 * bed file containing gene annotations (may change across cell types if using cell-type specific TSS's)
 * bed file containing chromosome annotations



## Download data

1. We downloaded the DNase-seq BAMs `ENCFF156LGK` (replicate 1) and `ENCFF134DLD` (replicate 2) from here: [K562 DNase-seq whole genome on GRCh38 genome assembly](https://www.encodeproject.org/search/?type=Experiment&status=released&perturbed=false&assay_title=DNase-seq&biosample_ontology.term_name=K562&assembly=GRCh38&files.file_type=bam&lab.title=John+Stamatoyannopoulos%2C+UW).
2. We use the ENCODE Blacklist for Human GRCh38 downloaded from [this link](https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz).
3. Histone ChIP-seq of K562 (`h3k27ac` mark, ENCODE experiment `ENCSR000AKP`) can be downloaded here: [K562 H3K27ac whole genome, ENCFF301TVL.bam](https://www.encodeproject.org/files/ENCFF301TVL/@@download/ENCFF301TVL).

As ENCODE annotations involve not only chromosomes but also other structures such as scaffolds, we need to report all these (195) structures in the chromosome annotations. We can find the loci of these scaffolds in header of the `ENCFF156LGK.bam` BAM file. We use `samtools` to convert the BAM to a SAM. `-H` indicates to keep the header only.

```bash
samtools view /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam -H -o /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK_header.sam

cp /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK_header.sam /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/ENCFF156LGK_header.sam

awk 'BEGIN{FS="\t"} {if($1=="@SQ"){split($2,id,":"); split($3,locus,":"); print(id[2]"\t"locus[2])}}' K562_whole_cell/references/ENCFF156LGK_header.sam > K562_whole_cell/references/chr_and_scaffolds_sizes
```

> ```
> work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/
> ├── input_data
> │   ├── Chromatin
> │   ├── Expression
> │   └── HiC
> └── references
> ```
>
> ```
> work2/project/regenet/data/
> ├── reads
> │   └── chipseq
> │       ├── h3k27ac
> │       │   └── homo_sapiens
> │       │       └── k562
> │       └── h3k4me3
> └── species
>     └── homo_sapiens
>         └── hg38
>             ├── gencv35
>             └── hg38-blacklist.v2.bed.gz
> ```
>
> ```
> /work2/project/regenet/results/
> ├── chipseq
> │   └── h3k27ac
> │       └── homo_sapiens
> │           └── hg38
> │               └── k562
> │                   └── ENCFF301TVL.bam
> └── dnaseseq
>     └── homo_sapiens
>         └── hg38
>             └── k562
>                 ├── ENCFF134DLD.bam
>                 ├── ENCFF156LGK.bam
>                 └── ENCFF156LGK_header.sam
> ```



## Run the ABC model

Don't forget to use `srun --pty bash` to connect to a node before processing small computations, and to use `slurm` for very resources-demanding computations. To get a job wall clock time: `seff <job id>`. This also yields resources info once the job is completed.

> ```bash
> $ seff 20187363
> Job ID: 20187363
> Cluster: genobull
> User/Group: thoellinger/U1220
> State: COMPLETED (exit code 0)
> Cores: 1
> CPU Utilized: 00:37:51
> CPU Efficiency: 96.56% of 00:39:12 core-walltime
> Job Wall-clock time: 00:39:12
> Memory Utilized: 1.31 GB
> Memory Efficiency: 2.05% of 64.00 GB
> ```

### Step 1: define candidate elements

> 'Candidate elements' are the set of putative enhancer elements for which ABC Scores will be computed. A typical way to define candidate elements is by calling peaks on a DNase-Seq or ATAC-Seq bam file. In this implementation we first call peaks using MACS2 and then process these peaks using `makeCandidateRegions.py`.

Here we call peaks on a DNase-Seq file, `ENCFF156LGK.bam`, which contains the alignment of the genome-wide sequencing of regions sensitive to cleavage by DNase I of K562 leukemia cells, over the GRCh38 genome assembly. Note that we use only one of the two bio-replicates for that purpose, as it was done by authors of the ABC model on their dataset (they did not explained why). **Note: next time I should verify if it is possible to add 2 "`-t`" args in `macs2 ` call.**

> ```bash
> /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/
> K562_whole_cell
> ├── hi_slurm
> │   └── step1.1.sh #see next cell
> ├── input_data
> │   ├── Chromatin
> │   ├── Expression
> │   └── HiC
> └── references
> ```

Write the following in `K562_whole_cell/hi_slurm/step1.1.sh` (**note: next time I reproduce this, I should take some time searching how one can parallelise `macs2` execution**):

> ```shell
> #!/bin/sh
> macs2 callpeak \
> -t /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam \
> -n ENCFF156LGK.macs2 \
> -f BAM \
> -g hs \
> -p .1 \
> --call-summits \
> --outdir /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/ABC_output/Peaks/
> ```
>

Now:

```bash
conda activate py2 && module load system/Python-2.7.2
```

```bash
sbatch --mem=4G --cpus-per-task=1 -J step1.1 --mail-user=tristan.hoellinger@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq K562_whole_cell/hi_slurm/step1.1.sh
```

```shell
conda deactivate && module unload system/Python-2.7.2
```

```bash
conda activate base && conda activate abcmodel && module load bioinfo/samtools-1.9 bioinfo/bedtools-2.27.1 bioinfo/tabix-0.2.5 bioinfo/juicer-1.5.6
```

Execution should take ~40 minutes.

> ```bash
> /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/
> K562_whole_cell
> ├── ABC_output
> │   └── Peaks #downstream files will be produced after running macs2 (step 1.1)
> │       ├── ENCFF156LGK.macs2_model.r
> │       ├── ENCFF156LGK.macs2_peaks.narrowPeak
> │       ├── ENCFF156LGK.macs2_peaks.xls
> │       └── ENCFF156LGK.macs2_summits.bed
> ├── hi_slurm
> │   └── step1.1.sh
> ├── input_data
> │   ├── Chromatin
> │   ├── Expression
> │   └── HiC
> └── references
> ```
>
> Now `makeCandidateRegions.py` will take as input the narrowPeak file produced by MACS2 and then perform the following processing steps:
>
> 1. Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
> 2. Resize each of these N peaks to be a fixed number of base pairs centered on the peak summit
> 3. Remove any blacklisted regions and include any whitelisted regions
> 4. Merge any overlapping regions

Problem: ENCODE does not only annotate chromosomes but also other structures such as scaffolds, resulting in `ENCFF156LGK.macs2_peaks.narrowPeak` containing names such as `chr11_KI270721v1_random` which are not associated to any locus in our reference genome file `chr_sizes`. We should either pre-filter the data to remove these occurrences or use a more complete genome file.

> ```bash
> $ awk '{print $1}' K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak |sort |uniq
> > chr1
> > chr10
> > chr11
> > chr11_KI270721v1_random/work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK_header.sam
> > chr12
> > chr13
> > chr14
> > chr14_GL000009v2_random
> > chr14_GL000194v1_random
> [...]
> ```

We should be able to find the loci of these scaffolds in header of the `ENCFF156LGK.bam` BAM file. We use `samtools` to convert the BAM to a SAM. `-H` indicates to keep the header only.

```bash
samtools view /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam -H -o /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK_header.sam

cp /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK_header.sam /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/chr_and_scaffolds_sizes

awk 'BEGIN{FS="\t"} {if($1=="@SQ"){split($2,id,":"); split($3,locus,":"); print(id[2]"\t"locus[2])}}' K562_whole_cell/references/chr_and_scaffolds_sizes > K562_whole_cell/references/chr_and_scaffolds_sizes_temp

mv K562_whole_cell/references/chr_and_scaffolds_sizes_temp K562_whole_cell/references/chr_and_scaffolds_sizes
```

> ```bash
> $ wc -l K562_whole_cell/references/chr_and_scaffolds_sizes
> 195 K562_whole_cell/references/chr_and_scaffolds_sizes
> ```

```bash
#First sort narrowPeak file
srun --pty bash
bedtools sort -faidx /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/chr_and_scaffolds_sizes -i /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak > /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak.sorted

#Then (better do on slurm, quite slow)
srun --pty bash
python src/makeCandidateRegions.py \
--narrowPeak /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak.sorted \
--bam /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam \
--outDir K562_whole_cell/ABC_output/Peaks/ \
--chrom_sizes /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/chr_and_scaffolds_sizes \
--regions_blacklist /work2/project/regenet/data/species/homo_sapiens/hg38/hg38-blacklist.v2.bed.gz \
--regions_whitelist /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 3000
```

> ```bash
> l K562_whole_cell/ABC_output/Peaks/
> ENCFF156LGK.macs2_model.r                                       ENCFF156LGK.macs2_peaks.narrowPeak.sorted.ENCFF156LGK.bam.Counts.bed
> ENCFF156LGK.macs2_peaks.narrowPeak                              ENCFF156LGK.macs2_peaks.xls
> ENCFF156LGK.macs2_peaks.narrowPeak.sorted                       ENCFF156LGK.macs2_summits.bed
> ENCFF156LGK.macs2_peaks.narrowPeak.sorted.candidateRegions.bed
> params.txt
> ```

### Step 2: quantifying enhancer activity

Now we need the [replicate 2](https://www.encodeproject.org/files/ENCFF134DLD/@@download/ENCFF134DLD.bam), the name of which is `ENCFF134DLD`.

```bash
wget https://www.encodeproject.org/files/ENCFF134DLD/@@download/ENCFF134DLD.bam -P /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/
```

Then write the following in `K562_whole_cell/hi_slurm/step2.1.sh`.

(One must first create a copy of chromosomes and scaffold sizes, with the same name but a BED extension: `chr_and_scaffolds_sizes` -> `chr_and_scaffolds_sizes.bed`)

```bash
#!/bin/sh
python src/run.neighborhoods.py \
--candidate_enhancer_regions K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed \
--H3K27ac /work2/project/regenet/results/chipseq/h3k27ac/homo_sapiens/hg38/k562/ENCFF301TVL.bam \
--DHS /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam,/work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF134DLD.bam \
--expression_table /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/example_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
--chrom_sizes /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/chr_and_scaffolds_sizes \
--ubiquitously_expressed_genes /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir K562_whole_cell/ABC_output/Neighborhoods/
```

```bash
sbatch --mem=16G --cpus-per-task=1 -J step2.1 --mail-user=tristan.hoellinger@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq K562_whole_cell/hi_slurm/step2.1.sh
```



=> Does not work. Raise the following error when launching on a node:

> ```bash
> $ python src/run.neighborhoods.py \
> > --candidate_enhancer_regions K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
> > --genes reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed \
> > --H3K27ac /work2/project/regenet/results/chipseq/h3k27ac/homo_sapiens/hg38/k562/ENCFF301TVL.bam \
> > --DHS /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam \
> > --expression_table /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/example_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
> > --chrom_sizes /work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/chr_and_scaffolds_sizes \
> > --ubiquitously_expressed_genes /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/UbiquitouslyExpressedGenesHG19.txt \
> > --cellType K562 \
> > --outdir K562_whole_cell/ABC_output/Neighborhoods/
> Namespace(ATAC='', DHS='/work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam', H3K27ac='/work2/project/regenet/results/chipseq/h3k27ac/homo_sapiens/hg38/k562/ENCFF301TVL.bam', candidate_enhancer_regions='K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak.sorted.candidateRegions.bed', cellType='K562', chrom_sizes='/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/K562_whole_cell/references/chr_and_scaffolds_sizes', default_accessibility_feature=None, enhancer_class_override=None, expression_table='/work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/example_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt', gene_name_annotations='symbol', genes='reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed', genes_for_class_assignment=None, outdir='K562_whole_cell/ABC_output/Neighborhoods/', primary_gene_identifier='symbol', qnorm=None, skip_gene_counts=False, skip_rpkm_quantile=False, supplementary_features=None, tss_slop_for_class_assignment=500, ubiquitously_expressed_genes='/work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/UbiquitouslyExpressedGenesHG19.txt', use_secondary_counting_method=False)
> Using gene expression from files: ['/work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/example_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt'] 
> 
> Traceback (most recent call last):
>   File "/home/thoellinger/.conda/envs/abcmodel/lib/python3.7/site-packages/pandas/core/indexes/base.py", line 2895, in get_loc
>     return self._engine.get_loc(casted_key)
>   File "pandas/_libs/index.pyx", line 70, in pandas._libs.index.IndexEngine.get_loc
>   File "pandas/_libs/index.pyx", line 101, in pandas._libs.index.IndexEngine.get_loc
>   File "pandas/_libs/hashtable_class_helper.pxi", line 1675, in pandas._libs.hashtable.PyObjectHashTable.get_item
>   File "pandas/_libs/hashtable_class_helper.pxi", line 1683, in pandas._libs.hashtable.PyObjectHashTable.get_item
> KeyError: 'end'
> 
> The above exception was the direct cause of the following exception:
> 
> Traceback (most recent call last):
>   File "src/run.neighborhoods.py", line 97, in <module>
>     main(args)
>   File "src/run.neighborhoods.py", line 93, in main
>     processCellType(args)
>   File "src/run.neighborhoods.py", line 74, in processCellType
>     outdir = args.outdir)
>   File "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/src/neighborhoods.py", line 85, in annotate_genes_with_features
>     tss1kb = make_tss_region_file(genes, outdir, genome_sizes)
>   File "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/src/neighborhoods.py", line 111, in make_tss_region_file
>     sizes_pr = df_to_pyranges(read_bed(sizes + '.bed'))
>   File "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/src/tools.py", line 54, in df_to_pyranges
>     df['End'] = df[end_col] + end_slop
>   File "/home/thoellinger/.conda/envs/abcmodel/lib/python3.7/site-packages/pandas/core/frame.py", line 2902, in __getitem__
>     indexer = self.columns.get_loc(key)
>   File "/home/thoellinger/.conda/envs/abcmodel/lib/python3.7/site-packages/pandas/core/indexes/base.py", line 2897, in get_loc
>     raise KeyError(key) from err
> KeyError: 'end'
> ```
>
> 
>
> ```bash
> File "/work2/project/regenet/workspace/thoellinger/ABC-Enhancer-Gene-Prediction/src/neighborhoods.py", line 85, in annotate_genes_with_features
>     tss1kb = make_tss_region_file(genes, outdir, genome_sizes)
> ```







Try with another genome file:

Maybe the problem is that we defined candidate elements over all scaffolds, not only chromosomes. If some have been found somewhere else than on a chromosome, there might be a problem looking for genes *in cis*, as these scaffolds are not mentioned in *genes* list (in `... .TSS500bp.bed`)

```bash
python src/run.neighborhoods.py \
--candidate_enhancer_regions K562_whole_cell/ABC_output/Peaks/ENCFF156LGK.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed \
--H3K27ac /work2/project/regenet/results/chipseq/h3k27ac/homo_sapiens/hg38/k562/ENCFF301TVL.bam \
--DHS /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam,/work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF134DLD.bam \
--expression_table /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/example_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
--chrom_sizes /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/chr_sizes \
--ubiquitously_expressed_genes /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir K562_whole_cell/ABC_output/Neighborhoods/
```

Seems to work when using a genome file that contains only chromosomes! Maybe the python code checks if a candidate element is in a well defined chromosome before searching for target genes, but does not verify that target genes are in a well defined chromosomes and yields an error if they're not.

Seems to be new errors though. After some time:

> ```bash
> samtools idxstats: fail to load index for "/work2/project/regenet/results/chipseq/h3k27ac/homo_sapiens/hg38/k562/ENCFF301TVL.bam", reverting to slow method
> Feature H3K27ac completed in 49.11177635192871
> Regenerating K562_whole_cell/ABC_output/Neighborhoods/Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph
> Counting coverage for Genes.DHS.ENCFF156LGK.bam
> Running: awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/chr_sizes <(samtools view -H /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam | grep SQ | cut -f 2 | cut -c 4- )  > K562_whole_cell/ABC_output/Neighborhoods/Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph.temp_sort_order
> Running: bedtools sort -faidx K562_whole_cell/ABC_output/Neighborhoods/Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph.temp_sort_order -i K562_whole_cell/ABC_output/Neighborhoods/GeneList.bed | bedtools coverage -g K562_whole_cell/ABC_output/Neighborhoods/Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph.temp_sort_order -counts -sorted -a stdin -b /work2/project/regenet/results/dnaseseq/homo_sapiens/hg38/k562/ENCFF156LGK.bam | awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}'  | bedtools sort -faidx /work2/project/regenet/workspace/sdjebali/egprediction/from1D/methods/abc.model/ABC-Enhancer-Gene-Prediction-0.2.2/reference/chr_sizes -i stdin > K562_whole_cell/ABC_output/Neighborhoods/Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph; rm K562_whole_cell/ABC_output/Neighborhoods/Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph.temp_sort_order
> 
> ```

Anyways, we tried on slurm and it seemed to work:

```bash
sbatch --mem=16G --cpus-per-task=1 -J step2.1bis --mail-user=tristan.hoellinger@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq K562_whole_cell/hi_slurm/step2.1bis.sh
```



> ```bash
> $ seff 20192684
> Job ID: 20192684
> Cluster: genobull
> User/Group: thoellinger/U1220
> State: COMPLETED (exit code 0)
> Cores: 1
> CPU Utilized: 00:48:20
> CPU Efficiency: 99.86% of 00:48:24 core-walltime
> Job Wall-clock time: 00:48:24
> Memory Utilized: 174.82 MB
> Memory Efficiency: 1.07% of 16.00 GB
> ```
>
> ```bash
> l K562_whole_cell/ABC_output/Neighborhoods/
> EnhancerList.bed                                       GeneList.txt
> EnhancerList.txt                                       Genes.DHS.ENCFF134DLD.bam.CountReads.bedgraph
> Enhancers.DHS.ENCFF134DLD.bam.CountReads.bedgraph      Genes.DHS.ENCFF156LGK.bam.CountReads.bedgraph
> Enhancers.DHS.ENCFF156LGK.bam.CountReads.bedgraph      Genes.H3K27ac.ENCFF301TVL.bam.CountReads.bedgraph
> Enhancers.H3K27ac.ENCFF301TVL.bam.CountReads.bedgraph  Genes.TSS1kb.DHS.ENCFF134DLD.bam.CountReads.bedgraph
> GeneList.bed                                           Genes.TSS1kb.DHS.ENCFF156LGK.bam.CountReads.bedgraph
> GeneList.TSS1kb.bed                                    Genes.TSS1kb.H3K27ac.ENCFF301TVL.bam.CountReads.bedgraph
> ```



### Step 3: computing the ABC score

If experimentally derived contact data is not available, one can run the ABC model using the powerlaw estimate only. In this case the ```--HiCdir``` argument should be excluded from ```predict.py``` and the ```--score_column powerlaw.Score``` argument should be included in ```predict.py```. In this case the ```ABC.Score``` column of the predictions file will be set to ```NaN```. The ```powerlaw.Score``` column of the output prediction files will be the relevant Score column to use.

```bash
srun --pty bash
python src/predict.py \
--enhancers K562_whole_cell/ABC_output/Neighborhoods/EnhancerList.txt \
--genes K562_whole_cell/ABC_output/Neighborhoods/GeneList.txt \
--score_column powerlaw.Score \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType K562 \
--outdir K562_whole_cell/ABC_output/Predictions/ \
--make_all_putative
```

> ```bash
> WARNING: Hi-C directory not provided. Model will only compute ABC score using powerlaw!
> reading genes
> reading enhancers
> Making predictions for chromosome: chr3
> Making putative predictions table...
> Completed chromosome: chr3. Elapsed time: 0.07643914222717285 
> 
> Making predictions for chromosome: chr5
> Making putative predictions table...
> Completed chromosome: chr5. Elapsed time: 0.06132197380065918 
> 
> Making predictions for chromosome: chr2
> Making putative predictions table...
> Completed chromosome: chr2. Elapsed time: 0.062340736389160156 
> 
> Making predictions for chromosome: chr4
> Making putative predictions table...
> Completed chromosome: chr4. Elapsed time: 0.046315908432006836 
> 
> Making predictions for chromosome: chr7
> Making putative predictions table...
> Completed chromosome: chr7. Elapsed time: 0.06442737579345703 
> 
> Making predictions for chromosome: chr6
> Making putative predictions table...
> Completed chromosome: chr6. Elapsed time: 0.14240288734436035 
> 
> Making predictions for chromosome: chr1
> Making putative predictions table...
> Completed chromosome: chr1. Elapsed time: 0.18433189392089844 
> 
> Writing output files...
> Done.
> ```

> ```bash
> $ l K562_whole_cell/ABC_output/Predictions/
> EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz	EnhancerPredictions.bedpe    EnhancerPredictions.txt  parameters.predict.txt
> EnhancerPredictionsAllPutative.txt.gz			EnhancerPredictionsFull.txt  GenePredictionStats.txt
> ```







We create a `sumstats` sub-directory in `ABC_output/Predictions` then compute summary statistics using `bedpe.sumstats.sh` script. We first need to import R.

```bash
module load system/R-3.6.2
../../../../../scripts/bedpe.sumstats.sh ../EnhancerPredictions.bedpe.gz "0.02-1" "500-500"
```



Does not work... permission denied. Let's do it locally.