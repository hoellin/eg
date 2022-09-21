# NOTEBOOK: Convert BENGI txt to custom bedpe containing enhancer and gene loci + TSS lists

*Jill E. Moore, Henry E. Pratt, Michael J. Purcaro et Zhiping Weng*

[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1924-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1924-8)

In this notebook we reprocess BENGI datasets into `.bedpe` containing comprehensive yet synthetic information for further statistical analysis with R.

## How to use this notebook?

First, make a copy of this notebook on your computer / cluster. The notebook can be found at [http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/](http://genoweb.toulouse.inra.fr/~thoellinger/notebooks/ipynb/).

To use Jupyter Notebook through ssh tunneling, one can refer to the following tutorial: [http://genoweb.toulouse.inra.fr/~thoellinger/notes/guidebooks/notebooks_ssh.html](http://genoweb.toulouse.inra.fr/~thoellinger/notes/guidebooks/notebooks_ssh.html)

Then, to use this notebook, one should only have to carefully modify the content of the "Set working directory" section, then to execute the notebook cell by cell, in the correct order. After execution of each cell, remember to check for errors before executing the next one.

## Introduction

In this notebook we reprocess BENGI datasets into `.bedpe` containing comprehensive yet synthetic information for further statistical analysis with R.


Note: Moore et al. used only 54,846 genes (or at least, provided TSS annotation for 54,846 genes), out of our 57,820 genes (coding + non-coding - of which 20,345 are protein-coding) in our whole GENCODE annotation. Yet we did not find any straightforward way to retrieve what filter they applied to evict 2,974 genes. Moreover, even restricting to these 54,846 genes, we find 175,554 TSS in our Gencode annotation (out of a total of 178,758 TSS), whereas the TSS list provided by Moore et al contains only 167,147 TSS. Again, we did not manage to find where this difference comes from.

In this notebook, we use the TSS list provided by Moore et al, not our own.

## Dependencies

[BENGI: Benchmark of Enhancer Gene Interactions](https://github.com/weng-lab/BENGI).

One need to download BENGI datasets from Weng lab, then to unzip all zipped files in `BENGI/Benchmark/Annotations/`.

```bash
gzip -cd GENCODEv19-TSSs.bed.gz > GENCODEv19-TSSs.bed
gzip -cd hg19-cCREs.bed.gz > hg19-cCREs.bed
```

## Data importation
### Set working directory


```python
# Inserm:
#root_dir = "/home/thoellinger/Documents/"
# Personal:
#root_dir = "/home/hoellinger/Documents/INSERM/"
#annotation_dir = 
# Genotoul:
annotation_dir = "/work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/"
work_dir = "/work2/project/regenet/results/multi/bengi/"

# When not on Genotoul:
#work_dir = root_dir+"BENGI/Benchmark/"
#annotation_dir = root_dir+"data/homo_sapiens/hg19.gencv19/"

path_to_TSSs_annotation = work_dir+"Annotations/GENCODEv19-TSSs.bed" # provided by Moore et al
path_to_ccREs_annotation = work_dir+"Annotations/hg19-cCREs.bed" # provided by Moore et al
path_to_genes_annotation = annotation_dir+"homo_sapiens.gtf" # full hg19 gencv19 annotation. Shall be processed.
path_to_benchmarks = work_dir+"All-Pairs.Natural-Ratio/"

nb_benchmarks = 6
file_names = list()
# benchmark files names without extensions
file_names.append("GM12878.CHiC-Benchmark.v3")
file_names.append("GM12878.CTCF-ChIAPET-Benchmark.v3")
file_names.append("GM12878.GEUVADIS-Benchmark.v3")
file_names.append("GM12878.GTEx-Benchmark.v3")
file_names.append("GM12878.HiC-Benchmark.v3")
file_names.append("GM12878.RNAPII-ChIAPET-Benchmark.v3")

# short custom names for benchmarks, same order as above
names = ["CHiC", "CTCF", "GEUVADIS", "GTEx", "HiC", "RNAPII"]

# Should have nothing to change below this line
# --------------------------------------------- 
col_names = ["ccRE", "gene", "interaction", "CV"] # column names in benchmarks
names_TSSs = ["chr", "start", "end", "transcript_id", "score", "strand", "gene_id"] # column names for TSSs annotation. "score" column is not used
names_ccREs = ["chr", "start", "end", "rDHS", "ccRE", "group"] # column names for ccREs. neither "rDHS" nor "group" are used.

files_list = list()
for k in range(nb_benchmarks):
    files_list.append(path_to_benchmarks+file_names[k]+".txt")
```

### Import standard libraries


```python
import pandas as pd
```

### Import files


```python
benchmarks = list()
for k in range(nb_benchmarks):
    benchmarks.append(pd.read_csv(files_list[k], sep='\t', header=None, names=col_names, dtype='str', engine='c'))
TSSs = pd.read_csv(path_to_TSSs_annotation, sep='\t', header=None, names=names_TSSs, dtype='str', engine='c')
ccREs = pd.read_csv(path_to_ccREs_annotation, sep='\t', header=None, names=names_ccREs, dtype='str', engine='c')
unprocessed_genes = pd.read_csv(path_to_genes_annotation, sep='\t', header=None, dtype='str', engine='c', skiprows=5)

unprocessed_genes = unprocessed_genes[unprocessed_genes.iloc[:,2]=="gene"]
unprocessed_genes.iloc[:,3] = (unprocessed_genes.iloc[:,3].astype(int) - 1).astype(str)
TSSs.iloc[:,[1, 2]] = (TSSs.iloc[:,[1, 2]].astype(int) - 1).astype(str)
```

## Data reprocessing
### `genes` annotation


```python
# The following filter keeps all gene types found in `hg19.gencv19/homo_sapiens.gtf`
# One can easily remove any unwanted gene type
list_of_gene_types = ["protein_coding", "pseudogene", "lincRNA", "antisense", "processed_transcript", "miRNA", "misc_RNA", "snRNA", "snoRNA", "polymorphic_pseudogene", "sense_intronic", "rRNA", "sense_overlapping", "IG_V_gene", "TR_V_gene", "IG_V_pseudogene", "TR_J_gene", "IG_C_gene", "IG_D_gene", "3prime_overlapping_ncrna", "TR_V_pseudogene", "IG_J_gene", "Mt_tRNA", "TR_C_gene", "IG_C_pseudogene", "TR_J_pseudogene", "TR_J_pseudogene", "TR_D_gene", "IG_J_pseudogene", "Mt_rRNA"]
filter_valid_genes = unprocessed_genes[8].str.split('; ').str[2].str.split(' \"').str[1].str.rstrip('\"').isin(list_of_gene_types)
unprocessed_genes = unprocessed_genes[filter_valid_genes]
```


```python
len(unprocessed_genes)
```


    57820


```python
col_genes = ["chr","start", "end", "gene_id", "strand"]
last_cols_unprocessed_genes = unprocessed_genes[8].str.split('; ').str[0].str.split(' \"')

genes = pd.concat([unprocessed_genes.iloc[:, [0,3,4]], last_cols_unprocessed_genes.str[1].str.rstrip('\"').to_frame(), unprocessed_genes.iloc[:,6]], axis=1)
genes.columns = col_genes
```


```python
genes.head(1)
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
      <th>chr</th>
      <th>start</th>
      <th>end</th>
      <th>gene_id</th>
      <th>strand</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>11868</td>
      <td>14412</td>
      <td>ENSG00000223972.4</td>
      <td>+</td>
    </tr>
  </tbody>
</table>
</div>


```python
len(genes)
```


    57820


```python
# make sure there are no scaffolds
#print(list(genes['chr'].drop_duplicates()))
# make sure gene ids are not degenerate
#print(len(genes['gene_id'].drop_duplicates()))
#genes.iloc[test] # find degenerate ones
```

### TSSs annotation -> 1 list / gene id


```python
names_compacted_TSSs = [names_TSSs[6], "TSSs"]
TSSs_lists = TSSs.iloc[:,[1, 6]]
TSSs_lists = TSSs_lists.groupby(TSSs_lists.columns[1])[TSSs_lists.columns[0]].apply(lambda x: ','.join(x)).to_frame().reset_index()
TSSs_lists.columns = names_compacted_TSSs
TSSs_lists.head(1)
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
      <th>gene_id</th>
      <th>TSSs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000003.10</td>
      <td>99891685,99891802,99894987</td>
    </tr>
  </tbody>
</table>
</div>


```python
len(TSSs_lists)
```


    54846

### Gather all useful information in 1 dataframe per file


```python
names_bedpe=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'interaction', 'strand1', 'strand2', 'TSSs']
# chrom1, start1 and end1 are relative to genes ; chrom2, start2, end2 to candidate regulatory elements

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
```


```python
bedpe[0]
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
      <th>chrom1</th>
      <th>start1</th>
      <th>end1</th>
      <th>chrom2</th>
      <th>start2</th>
      <th>end2</th>
      <th>name</th>
      <th>interaction</th>
      <th>strand1</th>
      <th>strand2</th>
      <th>TSSs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>536815</td>
      <td>659930</td>
      <td>chr1</td>
      <td>927848</td>
      <td>928157</td>
      <td>ENSG00000230021.3:EH37E0064164</td>
      <td>0</td>
      <td>-</td>
      <td>-</td>
      <td>655529,655573,655579,655579,659929</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>621058</td>
      <td>622053</td>
      <td>chr1</td>
      <td>927848</td>
      <td>928157</td>
      <td>ENSG00000185097.2:EH37E0064164</td>
      <td>0</td>
      <td>-</td>
      <td>-</td>
      <td>622052</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>657471</td>
      <td>660283</td>
      <td>chr1</td>
      <td>927848</td>
      <td>928157</td>
      <td>ENSG00000229376.3:EH37E0064164</td>
      <td>0</td>
      <td>+</td>
      <td>+</td>
      <td>657471</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>677192</td>
      <td>685396</td>
      <td>chr1</td>
      <td>927848</td>
      <td>928157</td>
      <td>ENSG00000235373.1:EH37E0064164</td>
      <td>0</td>
      <td>-</td>
      <td>-</td>
      <td>685395</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1</td>
      <td>677192</td>
      <td>685396</td>
      <td>chr1</td>
      <td>974131</td>
      <td>974485</td>
      <td>ENSG00000235373.1:EH37E0064194</td>
      <td>0</td>
      <td>-</td>
      <td>-</td>
      <td>685395</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>375723</th>
      <td>chrX</td>
      <td>154687177</td>
      <td>154687276</td>
      <td>chrX</td>
      <td>154459477</td>
      <td>154459854</td>
      <td>ENSG00000221603.1:EH37E1053818</td>
      <td>0</td>
      <td>+</td>
      <td>+</td>
      <td>154687177</td>
    </tr>
    <tr>
      <th>375724</th>
      <td>chrX</td>
      <td>154689079</td>
      <td>154689596</td>
      <td>chrX</td>
      <td>154459477</td>
      <td>154459854</td>
      <td>ENSG00000185978.4:EH37E1053818</td>
      <td>0</td>
      <td>-</td>
      <td>-</td>
      <td>154689595</td>
    </tr>
    <tr>
      <th>375725</th>
      <td>chrX</td>
      <td>154695630</td>
      <td>154841277</td>
      <td>chrX</td>
      <td>154459477</td>
      <td>154459854</td>
      <td>ENSG00000224533.3:EH37E1053818</td>
      <td>0</td>
      <td>+</td>
      <td>+</td>
      <td>154695630,154696200,154718983</td>
    </tr>
    <tr>
      <th>375726</th>
      <td>chrX</td>
      <td>154697946</td>
      <td>154716707</td>
      <td>chrX</td>
      <td>154459477</td>
      <td>154459854</td>
      <td>ENSG00000225393.1:EH37E1053818</td>
      <td>0</td>
      <td>+</td>
      <td>+</td>
      <td>154697946</td>
    </tr>
    <tr>
      <th>375727</th>
      <td>chrX</td>
      <td>154719775</td>
      <td>154899605</td>
      <td>chrX</td>
      <td>154459477</td>
      <td>154459854</td>
      <td>ENSG00000185973.6:EH37E1053818</td>
      <td>0</td>
      <td>-</td>
      <td>-</td>
      <td>154722149,154774937,154842537,154842555,154842...</td>
    </tr>
  </tbody>
</table>
<p>375728 rows × 11 columns</p>
</div>

### Export new dataframes to bedpe


```python
#for k in range(nb_benchmarks):
#    bedpe[k].to_csv(path_to_benchmarks+file_names[k]+".new.bedpe",sep='\t',header=False,index=False)
#    with open(path_to_benchmarks+file_names[k]+'.header.new.bedpe.txt', 'w') as f:
#        f.write("%s\n" % ', '.join(list(bedpe[k].columns)))
```

