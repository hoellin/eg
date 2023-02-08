#!/usr/bin/env python

"""BENGI_to_bedpe.py: This script converts BENGI datasets (see "A curated benchmark of enhancer-gene interactions for evaluating enhancer-target gene prediction methods" 2020 Moore et al.)
to a more comprehensive bedpe-like format better-suited for basic statistical analysis."""

import pandas as pd

# =====================
# Set working directory
# =====================

# Uncomment "Inserm", "Personnal" or "Genotoul" -marked lines for pre-define environment,
# or use custom paths.

## Inserm
#root_dir = "/home/thoellinger/Documents/"

## Personal
#root_dir = "/home/hoellinger/Documents/INSERM/"

## MacOS
root_dir = "/Users/hoellinger/Documents/INSERM/"

## Genotoul
#TODO

name_gene_annotation = "gencode.v19.annotation.gtf"

work_dir = root_dir+"BENGI/Benchmark/"
path_to_TSSs_annotation = work_dir+"Annotations/GENCODEv19-TSSs.bed" # hg19v19
path_to_ccREs_annotation = work_dir+"Annotations/hg19-cCREs.bed" # hg19
path_to_genes_annotation = root_dir+"data/homo_sapiens/hg19.gencv19/"+name_gene_annotation
path_to_benchmarks = work_dir+"All-Pairs.Natural-Ratio/"

nb_benchmarks = 6
file_names = list()

# Names (without extension) of BENGI benchmarks to reprocess (will not be overwritten)
file_names.append("GM12878.CHiC-Benchmark.v3")
file_names.append("GM12878.CTCF-ChIAPET-Benchmark.v3")
file_names.append("GM12878.GEUVADIS-Benchmark.v3")
file_names.append("GM12878.GTEx-Benchmark.v3")
file_names.append("GM12878.HiC-Benchmark.v3")
file_names.append("GM12878.RNAPII-ChIAPET-Benchmark.v3")

# Short custom names for benchmarks, same order as above
names = ["CHiC", "CTCF", "GEUVADIS", "GTEx", "HiC", "RNAPII"]

# Should have nothing to change below this line
col_names = ["ccRE", "gene", "interaction", "CV"] # column names in benchmark dataframes
names_TSSs = ["chr", "start", "end", "transcript_id", "score", "strand", "gene_id"] # column names for TSSs annotation. "transcript", "score",  "strand" not used
names_ccREs = ["chr", "start", "end", "rDHS", "ccRE", "group"] # column names for ccREs. neither "rDHS" nor "group" are used.

files_list = list()
for k in range(nb_benchmarks):
	files_list.append(path_to_benchmarks+file_names[k]+".txt")

benchmarks = list()
for k in range(nb_benchmarks):
	benchmarks.append(pd.read_csv(files_list[k], sep='\t', header=None, names=col_names, dtype='str', engine='c'))
TSSs = pd.read_csv(path_to_TSSs_annotation, sep='\t', header=None, names=names_TSSs, dtype='str', engine='c')
ccREs = pd.read_csv(path_to_ccREs_annotation, sep='\t', header=None, names=names_ccREs, dtype='str', engine='c')
unprocessed_genes = pd.read_csv(path_to_genes_annotation, sep='\t', header=None, dtype='str', engine='c', skiprows=5)

unprocessed_genes = unprocessed_genes[unprocessed_genes.iloc[:,2]=="gene"]
unprocessed_genes.iloc[:,3] = (unprocessed_genes.iloc[:,3].astype(int) - 1).astype(str)
TSSs.iloc[:,[1, 2]] = (TSSs.iloc[:,[1, 2]].astype(int) - 1).astype(str)

list_of_gene_types = ["protein_coding", "pseudogene", "lincRNA", "antisense", "processed_transcript", "miRNA", "misc_RNA", "snRNA", "snoRNA", "polymorphic_pseudogene", "sense_intronic", "rRNA", "sense_overlapping", "IG_V_gene", "TR_V_gene", "IG_V_pseudogene", "TR_J_gene", "IG_C_gene", "IG_D_gene", "3prime_overlapping_ncrna", "TR_V_pseudogene", "IG_J_gene", "Mt_tRNA", "TR_C_gene", "IG_C_pseudogene", "TR_J_pseudogene", "TR_J_pseudogene", "TR_D_gene", "IG_J_pseudogene", "Mt_rRNA"]
filter_valid_genes = unprocessed_genes[8].str.split('; ').str[2].str.split(' \"').str[1].str.rstrip('\"').isin(list_of_gene_types)
unprocessed_genes = unprocessed_genes[filter_valid_genes]

col_genes = ["chr","start", "end", "gene_id", "strand"]
last_cols_unprocessed_genes = unprocessed_genes[8].str.split('; ').str[0].str.split(' \"')

genes = pd.concat([unprocessed_genes.iloc[:, [0,3,4]], last_cols_unprocessed_genes.str[1].str.rstrip('\"').to_frame(), unprocessed_genes.iloc[:,6]], axis=1)
genes.columns = col_genes

names_compacted_TSSs = [names_TSSs[6], "TSSs"]
TSSs_lists = TSSs.iloc[:,[1, 6]]
TSSs_lists = TSSs_lists.groupby(TSSs_lists.columns[1])[TSSs_lists.columns[0]].apply(lambda x: ','.join(x)).to_frame().reset_index()
TSSs_lists.columns = names_compacted_TSSs

names_bedpe=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'interaction', 'strand1', 'strand2', 'TSSs']

new_names = {col_genes[0]: names_bedpe[0], col_genes[1]: names_bedpe[1], col_genes[2]: names_bedpe[2], col_genes[3]: names_bedpe[6], col_genes[4]: names_bedpe[8]}
bedpe = [None]*nb_benchmarks    
for k in range(nb_benchmarks):
	bedpe[k] = genes.copy().merge(benchmarks[k].drop(col_names[3], axis=1), how='inner', left_on=col_genes[3], right_on=col_names[1]).merge(TSSs_lists, how='inner', left_on=col_genes[3], right_on=names_TSSs[-1]).merge(ccREs.iloc[:,[0,1,2,4]].rename(columns={names_ccREs[0]: names_bedpe[3], names_ccREs[1]: names_bedpe[4], names_ccREs[2]: names_bedpe[5]}), how='left', left_on=col_names[0], right_on=names_ccREs[4])
	bedpe[k].iloc[:,3] += ':'+bedpe[k].iloc[:,5].map(str)
	bedpe[k].drop([col_names[0], col_names[1]], axis=1, inplace=True)
	bedpe[k][col_names[2]] = bedpe[k][col_names[2]].str.rstrip(' ')
	bedpe[k].insert(5, names_bedpe[9], bedpe[k].iloc[:,4])
	bedpe[k].rename(columns=new_names, inplace=True)
	cols = bedpe[k].columns.tolist()
	new_cols = cols[:3] + cols[-3:] + [cols[3], cols[6]] + cols[4:6] + [cols[7]]
	bedpe[k] = bedpe[k][new_cols]

for k in range(nb_benchmarks):
	bedpe[k].to_csv(path_to_benchmarks+file_names[k]+".new.bedpe",sep='\t',header=False,index=False)
	with open(path_to_benchmarks+file_names[k]+'.header.new.bedpe.txt', 'w') as f:
		f.write("%s\n" % ', '.join(list(bedpe[k].columns)))

__author__ = "Tristan Hoellinger"
__credits__ = ["Sarah Djebali"]
__email__ = "tristan.hoellinger@inserm.fr"