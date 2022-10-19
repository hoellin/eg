# Intersect enhancers suspected to regulate the 13 initial genes (involved in hemochromatosis or iron metabolism regulation), with the NHGRI-EBI GWAS Catalog

## Data acquisition

### Download the latest GWAS catalog

The version information on the GWAS website were the following when I downloaded the catalog:

> As of 2021-10-22, the GWAS Catalog contains 5419 publications and 316782 associations. GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p13 and dbSNP Build 154.

```bash
cd /work2/project/regenet/results/gwas/gwas.catalog
```

```bash
wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O /work2/project/regenet/results/gwas/gwas.catalog/v1.0.2_r2021-10-02/gwas_catalog_v1.0.2-associations_r2021-10-02.tsv
```

### List of enhancers suspected to regulate genes involved in hemochromatosis

We obtained it while computing networks with the R markdown.

```bash
cd /work2/project/regenet/results/multi/abc.model/Nasser2021/
mv GWAS/ GWAS_10_initial_genes_fall_2021/
mkdir GWAS_13_initial_genes/ && cd GWAS_13_initial_genes/
```

```bash
cp /work2/project/regenet/workspace/thoellinger/shared/2022/networks_hemochromatosis/results/enhancers.list .
```

### List of enhancers suspected to regulate the 457 original+new genes

```bash
cd /work2/project/regenet/results/multi/abc.model/Nasser2021/GWAS_13_initial_genes/
cp /work2/project/regenet/workspace/thoellinger/shared/2022/networks_hemochromatosis/results/enhancers.new_genes.list .
```

## Data pre-processing

### Enhancer lists

We want to convert them in the bed format so we do the following:

```bash
conda activate base && module load bioinfo/bedtools-2.27.1
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"}; {split($1,l1,":"); split(l1[2],l2,"-"); print l1[1], l2[1], l2[2]}' enhancers.list |bedtools sort -faidx ../chr_sizes -i stdin > enhancers.sorted.bed
```

```bash
awk 'BEGIN{FS="\t"; OFS="\t"}; {split($1,l1,":"); split(l1[2],l2,"-"); print l1[1], l2[1], l2[2]}' enhancers.new_genes.list |bedtools sort -faidx ../chr_sizes -i stdin > enhancers.new_genes.sorted.bed
```

### SNP related to hemochromatosis

```bash
cd /work2/project/regenet/results/gwas/gwas.catalog/v1.0.2_r2021-10-02/
```

We want to extract from the GWAS catalog only the SNP related to hemochromatosis. Let us first investigate the list of all diseases / traits in the GWAS catalog.

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {traits[$8]++} END{for(u in traits){print u, traits[u]}}' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |wc -l
> ```
>
> 5845

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {traits[$8]++} END{for(u in traits){print u, traits[u]}}' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |sort -nrk2,2 |head -n 30
> ```
>
> ```bash
> Height	6144
> Schizophrenia	2161
> Hematocrit	1913
> Triglycerides	1487
> Neuroticism	1171
> Asthma	1168
> Weight	1037
> Plateletcrit	1014
> Hemoglobin	973
> Intelligence	581
> Chronotype	479
> Psoriasis	422
> Insomnia	422
> Hypertension	308
> Hypothyroidism	291
> Depression	290
> Morningness	274
> Tonsillectomy	241
> Eczema	214
> FEV1	206
> Neurociticism	181
> Longevity	169
> Adventurousness	167
> Glaucoma	162
> Snoring	158
> Cataracts	154
> Melanoma	149
> Endometriosis	148
> Worry	147
> Migraine	146
> ```

> ```bash
> awk -F "\t" '$8~/([Hh]emochromatosis)/' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |wc -l
> ```
>
> 2

```bash
grep emochromatos gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |awk 'BEGIN{FS="\t"; OFS="\t"} {traits[$8]++} END{for(u in traits){print u, traits[u]}}' |sort -t $'\t' -nrk2,2
```

> ```bash
> Iron status biomarkers (transferrin levels)	9
> Iron status biomarkers (transferrin saturation)	6
> Iron status biomarkers (iron levels)	5
> Iron status biomarkers (ferritin levels)	5
> Hereditary hemochromatosis-related traits (HFE mutation homozygotes)	2
> ```

Only 2 entries involve hemochromatosis directly, and a few others (25) involve studies about hemochromatosis but not SNP directly related to hemochromatosis. We check that no entries correspond to the other name of hemochromatosis, namely iron overload:

> ```bash
> grep 'iron.*overload' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |wc -l
> ```
>
> 0

We conclude that there are too few data to focus on those involving hemochromatosis. We will rather focus on those for which the trait is `Iron status biomarkers`

### SNP related to iron in general

> ```bash
> grep '[iI]ron.*biomarkers' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |wc -l
> ```
>
> 250

Well we use the following filter to be sure not to miss any trait involving iron, then we will filter out by hand those which do not actually correspond to iron (such as "env**iron**ment") or that do not interest us:

> ```bash
> awk -F "\t" '$8~/([iI]ron)/' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |wc -l
> ```
>
> 378

```bash
awk -F "\t" 'if($8~/([iI]ron)/ {traits[$8]++}; END{for(u in traits){print u, traits[u]}}' gwas_catalog_v1.0.2-associations_r2021-10-02.tsv |sort -t $'\t' -nrk2,2
```

> ```bash
> Iron status biomarkers (ferritin levels)	55
> Iron status biomarkers	54
> Iron status biomarkers (total iron binding capacity)	48
> Sensitivity to environmental stress and adversity	47
> Iron status biomarkers (transferrin saturation)	47
> Iron status biomarkers (iron levels)	37
> Liver iron content	23
> Economic and political preferences (environmentalism)	17
> BMI x  environmental factors (excluding physical activity) interaction	11
> Iron status biomarkers (transferrin levels)	9
> Forced expiratory volume in 1 second (occupational environmental exposures interaction)	7
> BMI x environmental factors (including physical activity) interaction	7
> Iron deficiency anemia	5
> Type 2 diabetes (dietary heme iron intake interaction)	4
> Forced expiratory volume in 1 second (environmental tobacco smoke interaction)	3
> Pancreas iron content	1
> Iron levels	1
> Iron deficiency	1
> Cleft plate (environmental tobacco smoke interaction)	1
> ```

So the only ones that could interest us are the following:

> ```bash
> Iron status biomarkers (ferritin levels)    55
> Iron status biomarkers    54
> Iron status biomarkers (total iron binding capacity)    48
> Iron status biomarkers (transferrin saturation)    47
> Iron status biomarkers (iron levels)    37
> Liver iron content    23
> Iron status biomarkers (transferrin levels)    9
> Iron deficiency anemia    5
> Type 2 diabetes (dietary heme iron intake interaction)    4
> Pancreas iron content    1
> Iron levels    1
> Iron deficiency    1
> Hereditary hemochromatosis-related traits (HFE mutation homozygotes) 2
> ```

We remove those involved in diabetes and pancreas, leading to the following list:

> ```bash
> Iron status biomarkers (ferritin levels)
> Iron status biomarkers
> Iron status biomarkers (total iron binding capacity)
> Iron status biomarkers (transferrin saturation)
> Iron status biomarkers (iron levels)
> Liver iron content
> Iron status biomarkers (transferrin levels)
> Iron deficiency anemia
> Iron levels
> Iron deficiency
> Hereditary hemochromatosis-related traits (HFE mutation homozygotes)
> ```

So now let's go back to our working directory and extract only the lines of interest:

```bash
cd /work2/project/regenet/results/multi/abc.model/Nasser2021/GWAS_13_initial_genes
gwas="/work2/project/regenet/results/gwas/gwas.catalog/v1.0.2_r2021-10-02/gwas_catalog_v1.0.2-associations_r2021-10-02.tsv"
```

```bash
awk -F "\t" '$8~/(Iron status biomarkers)|(Liver iron content)|(Iron deficiency anemia)|(Iron levels)|(Iron deficiency)|(Hereditary hemochromatosis-related traits (HFE mutation homozygotes))/' $gwas > associations_full.tsv
```

`associations_full.tsv` contains 280 rows. Remark: we see in the following section that only 269 of those lines are valid.

## Perform intersection (will enhancers regulating the 10 initial genes)

### TSV to bed

Fields $12 and $13 in `associations_full.tsv` are <chr id> and <chr pos> of the SNP. The 3 first columns of the bed file shall be the position of the SNP, and in column 4, as feature name, we concatenate all the columns of `associations_full.tsv`  with a separator like "::", so that the output is strictly in bed format.

```bash
awk 'BEGIN{FS="\t"; OFS="\t"} {if($13 && $14){ORS="\t"; print "chr"$12, $13, $13+1; ORS="\n"; gsub(/\t/,"::",$0); print}}' associations_full.tsv |bedtools sort -faidx ../chr_sizes -i stdin > associations_full.bed
```

### Liftover

The Genome Assembly used in our GWAS catalog is GRCh38.p13 whereas it's GRCh37/hg19 in the Nasser2021 predictions, so we need to start with a liftover.

In order to do this we need several modules, that should contain the keywords `kent` and `utils`. Note that the output of `module av` is actually directed to stderr so we use `2>&1` to redirect stderr to stdout, then `>/dev/null` to redirect stdout to `/dev/null`.

```bash
(module av) 2>&1 >/dev/null |grep kent
```

> ```bash
> bioinfo/kentUtils-302.1.0
> bioinfo/kentUtils-v370
> ```

So we can finally load those modules:

```bash
module load bioinfo/kentUtils-302.1.0 bioinfo/kentUtils-v370
```

(actually `bioinfo/kentUtils-v370` should suffice, cf mail from Sarah 07/12/2020, 14:34).

Now we can simply use `liftOver`:

> ```bash
> liftOver
> ```
>
> ```
> liftOver - Move annotations from one assembly to another
> usage:
>    liftOver oldFile map.chain newFile unMapped
> oldFile and newFile are in bed format by default, but can be in GFF and
> maybe eventually others with the appropriate flags below.
> The map.chain file has the old genome as the target and the new genome
> as the query.
> ```

`unlifted.bed` file is not an optional argument and will contain all genome positions that cannot be lifted in case it happens for some. Full documentation here: https://genome.sph.umich.edu/wiki/LiftOver

One can find the `chain` file that goes from hg38 to hg19 here: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/. Let's download it:

```bash
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz .
```

> Not this chain file was last modified at the date: 2013-12-31 17:24

So now we can perform the liftover. `enhancers.sorted.bed` is already mapped over the hg19 assembly, so all we need to do is to convert `associations_full.bed` to hg19.

```bash
liftOver -bedPlus=1 -tab associations_full.bed hg38ToHg19.over.chain associations_full.hg19_liftover.bed unlifted.bed
```

OK.

Note that using `meld` we found that all rows are different between `associations_full.bed` and `associations_full.hg19_liftover.bed`.

### Intersect

```bash
srun --mem=32G --pty bash
conda activate base && module load bioinfo/bedtools-2.27.1
```

```bash
bedtools intersect -wo -a enhancers.sorted.bed -b associations_full.hg19_liftover.bed > enhancers.overlap_GWAS.bedpe
```

> ```bash
> wc -l enhancers.overlap_GWAS.bedpe
> ```
>
> 4
>
> ```bash
> cat enhancers.overlap_GWAS.bedpe
> ```
>
> ```bash
> chr2	190405339	190407623	chr2	190407500	190407501	2021-08-13::34128465::Liu Y::2021-06-15::Elife::www.ncbi.nlm.nih.gov/pubmed/34128465::Genetic architecture of 11 organ traits derived from abdominal MRI using deep learning.::Liver iron content::32,858 European ancestry individuals::NA::2q32.2::2::189542774::SLC40A1::KDM3AP1 - SLC40A1::ENSG00000233996::ENSG00000138449::::54782::17816::rs115380467-?::rs115380467::0::115380467::regulatory_region_variant::1::NR::3E-12::11.522878745280337::(conditional)::0.181::NR mg/g increase::Affymetrix [9390170] (imputed)::N::liver iron measurement::http://www.ebi.ac.uk/efo/EFO_0010056::GCST90016674::Genome-wide genotyping array	1
> chr3	195799723	195801054	chr3	195800811	195800812	2021-03-22::33536631::Bell S::2021-02-03::Commun Biol::www.ncbi.nlm.nih.gov/pubmed/33536631::A genome-wide meta-analysis yields 46 new loci associating with biomarkers of iron homeostasis.::Iron status biomarkers (total iron binding capacity)::135,430 European ancestry individuals::NA::3q29::3::196073940::TFRC::TFRC::::::ENSG00000072274::::::rs3817672-C::rs3817672::0::3817672::missense_variant::0::0.44::6E-11::10.221848749616356::::0.031::[0.04-0.022] SD decrease::Affymetrix, Illumina [40000000] (imputed)::N::total iron binding capacity::http://www.ebi.ac.uk/efo/EFO_0006334::GCST011368::Genome-wide genotyping array	1
> chr3	195799723	195801054	chr3	195800811	195800812	2021-03-22::33536631::Bell S::2021-02-03::Commun Biol::www.ncbi.nlm.nih.gov/pubmed/33536631::A genome-wide meta-analysis yields 46 new loci associating with biomarkers of iron homeostasis.::Iron status biomarkers (transferrin saturation)::131,471 European ancestry individuals::NA::3q29::3::196073940::TFRC::TFRC::::::ENSG00000072274::::::rs3817672-C::rs3817672::0::3817672::missense_variant::0::0.44::2E-8::7.698970004336019::::0.026::[0.017-0.034] SD increase::Affymetrix, Illumina [40000000] (imputed)::N::transferrin saturation measurement::http://www.ebi.ac.uk/efo/EFO_0006333::GCST011366::Genome-wide genotyping array	1
> chr22	37505175	37505880	chr22	37505552	37505553	2015-06-14::25352340::Benyamin B::2014-10-29::Nat Commun::www.ncbi.nlm.nih.gov/pubmed/25352340::Novel loci affecting iron homeostasis and their effects in individuals at risk for hemochromatosis.::Iron status biomarkers (iron levels)::23,986 European ancestry individuals::Up to 24,986 European ancestry individuals::22q12.3::22::37109512::TMPRSS6::TMPRSS6::::::ENSG00000187045::::::rs228916-T::rs228916::0::228916::5_prime_UTR_variant::0::0.875::3E-8::7.522878745280337::::0.086::[0.055-0.117] unit decrease::Affymetrix, Illumina [~ 2100000] (imputed)::N::iron biomarker measurement, serum iron measurement::http://www.ebi.ac.uk/efo/EFO_0004461, http://www.ebi.ac.uk/efo/EFO_0006332::GCST002679::Genome-wide genotyping array	1
> ```

Well, there are only 4 overlap, involving only 3 unique enhancers. This is expected, as we are searching for overlap between a list of only 137 enhancer over the whole genome, and a list of only 269 SNP over the whole genome.

## Perform intersection (with enhancers regulating the 457 original + newly found genes)

```bash
srun --mem=32G --pty bash
conda activate base && module load bioinfo/bedtools-2.27.1
```

```bash
bedtools intersect -wo -a enhancers.new_genes.sorted.bed -b associations_full.hg19_liftover.bed > enhancers.new_genes.overlap_GWAS.bedpe
```

There is only one more overlap, namely:

> ```bash
> chr6	26103318	26105346	chr6	26104632	26104633	2021-08-13::34128465::Liu Y::2021-06-15::Elife::www.ncbi.nlm.nih.gov/pubmed/34128465::Genetic architecture of 11 organ traits derived from abdominal MRI using deep learning.::Liver iron content::32,858 European ancestry individuals::NA::6p22.2::6::26104404::HIST1H4C::H4C3 - H1-6::ENSG00000197061::ENSG00000187475::::67::3008::rs198851-?::rs198851::0::198851::TF_binding_site_variant::1::NR::3E-73::72.52287874528034::(conditional)::0.194::NR mg/g decrease::Affymetrix [9390170] (imputed)::N::liver iron measurement::http://www.ebi.ac.uk/efo/EFO_0010056::GCST90016674::Genome-wide genotyping array	1
> ```

## Are the genes regulated by the enhancers found to intersect a GWAS SNP, the same as those given in GWAS as associated to the SNP?

We did this analysis starting from 10 instead of 13 genes, but as it conduces to only 1 more overlap.