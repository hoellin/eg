# Summary statistics on Nasser et al. 2021 ABC predictions in our 5 biosamples of interest

## Resources

### Data

All data we use are available in `/work2/project/regenet/results/multi/abc.model/Nasser2021/` (private).

### Code

We use Sarah Djebali's `bedpe.sumstats.sh` script. See example [here](../../../guidebooks/script_bedpe_summary/script_bedpe_summary/#using-bedpesumstatssh).

R packages `optparse` and `reshape2` are required. It suffices to load the following module:

```bash
module load system/R-4.0.4_gcc-9.3.0
```

### Miscellaneous

The script `bedpe.sumstats.sh` requires a lot of memory. By default, with `srun --pty bash`, we have too little memory for it not to be killed, so better run:

```bash
srun --mem=32G --pty bash
```

## Compute summary statistics

The idea is to do the following:

> ```bash
> srun --mem=32G --pty bash
> ```
>
> ```bash
> module load system/R-4.0.4_gcc-9.3.0
> ```
>
> ```bash
> /work2/project/regenet/workspace/thoellinger/scripts/bedpe.sumstats.sh ../Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.bedpe.gz $(awk 'BEGIN{FS="\t"; min=1; max=0;} {if($8>max){max=$8}; if($8<min){min=$8}} END{print min"-"max;}' ../Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.bedpe) "500-500"
> ```

for each file of interest. That would be quite long, so, we do the following once and for all:

```bash
srun --mem=32G --pty bash
module load system/R-4.0.4_gcc-9.3.0
```

```bash
declare -a files=("Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_intestine_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_liver_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.hepatocyte-ENCODE.all_putative_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.liver-ENCODE.all_putative_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.HepG2-Roadmap.all_putative_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.large_intestine_fetal-Roadmap.all_putative_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.small_intestine_fetal-Roadmap.all_putative_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.hepatocyte-ENCODE.all_putative_enhancers.merged_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.liver-ENCODE.all_putative_enhancers.merged_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.HepG2-Roadmap.all_putative_enhancers.merged_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.large_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers.sorted.bedpe"
"Nasser2021ABCPredictions.small_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers.sorted.bedpe")

cd /work2/project/regenet/results/multi/abc.model/Nasser2021/sumstats
for file in "${files[@]}"
do
	gzip -c "../$file" > "../$file.gz"
   	dir_name=$(echo "$file" | sed -e 's/Nasser2021ABCPredictions.//g' -e 's/.sorted.bedpe//g')
   	mkdir "../sumstats_results/$dir_name"
   	/work2/project/regenet/workspace/thoellinger/scripts/bedpe.sumstats.sh "../$file.gz" $(awk 'BEGIN{FS="\t"; min=1; max=0;} {if($8>max){max=$8}; if($8<min){min=$8}} END{print min"-"max;}' "../$file") "500-500"
   	mv ./* "../sumstats_results/$dir_name"
done
```

For testing purpose, whenever necessary, one can use the following:

> ```bash
> declare -a files=("Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.all_biosamples.all_putative_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.intestine.all_putative_enhancers.merged_intestine_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.intestine.all_putative_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver.all_putative_enhancers.merged_liver_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver.all_putative_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers.sorted.bedpe" "Nasser2021ABCPredictions.liver_and_intestine.all_putative_enhancers.sorted.bedpe")
> 
> 
> for file in "${files[@]}"
> do
> 	dir=$(echo "$file" | sed -e 's/Nasser2021ABCPredictions.//g' -e 's/.sorted.bedpe//g')
> 	echo $dir
> #   echo "${file//.sorted.bedpe/}"
> done
> ```



## Results

### Distances of Enhancer-TSS pairs

| Biosample Name \ Enhancers type |                      Original Enhancers                      | Merged enhancers (across considered biosamples)              | Merged enhancers (across all 131 biosamples)                 |
| :-----------------------------: | :----------------------------------------------------------: | ------------------------------------------------------------ | ------------------------------------------------------------ |
|       All 131 biosamples        | ![](sumstats_results/all_biosamples.all_putative_enhancers/Distance.png) | *Same as column 3*                                           | ![](sumstats_results/all_biosamples.all_putative_enhancers.merged_enhancers/Distance.png) |
|       Liver and Intestine       | ![](sumstats_results/liver_and_intestine.all_putative_enhancers/Distance.png) | ![](sumstats_results/liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers/Distance.png) | ![](sumstats_results/liver_and_intestine.all_putative_enhancers.merged_enhancers/Distance.png)) |
|              Liver              | ![](sumstats_results/liver.all_putative_enhancers/Distance.png) | ![](sumstats_results/liver.all_putative_enhancers.merged_liver_enhancers/Distance.png) | ![](sumstats_results/liver.all_putative_enhancers.merged_enhancers/Distance.png) |
|            Intestine            | ![](sumstats_results/intestine.all_putative_enhancers/Distance.png) | ![](sumstats_results/intestine.all_putative_enhancers.merged_intestine_enhancers/Distance.png) | ![](sumstats_results/intestine.all_putative_enhancers.merged_enhancers/Distance.png) |
|        hepatocyte-ENCODE        | ![](sumstats_results/hepatocyte-ENCODE.all_putative_enhancers/Distance.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/hepatocyte-ENCODE.all_putative_enhancers.merged_enhancers/Distance.png) |
|          HepG2-Roadmap          | ![](sumstats_results/HepG2-Roadmap.all_putative_enhancers/Distance.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/HepG2-Roadmap.all_putative_enhancers.merged_enhancers/Distance.png) |
|          liver-ENCODE           | ![](sumstats_results/liver-ENCODE.all_putative_enhancers/Distance.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/liver-ENCODE.all_putative_enhancers.merged_enhancers/Distance.png) |
|  large_intestine_fetal-Roadmap  | ![](sumstats_results/large_intestine_fetal-Roadmap.all_putative_enhancers/Distance.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/large_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers/Distance.png) |
|  small_intestine_fetal-Roadmap  | ![](sumstats_results/small_intestine_fetal-Roadmap.all_putative_enhancers/Distance.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/small_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers/Distance.png) |



### Distribution of nb of connections of enhancers (red) / of TSS (green)

Number of connections of element 1 (enhancers) / of element 2 (TSS) => TSS usually makes more connections to   enhancers, than enhancers make connections to TSS (the four screens distinguish between the ABC.score -quartiles).

| Biosample Name \ Enhancers type |                      Original Enhancers                      | Merged enhancers (across considered biosamples)              | Merged enhancers (across all 131 biosamples)                 |
| :-----------------------------: | :----------------------------------------------------------: | ------------------------------------------------------------ | ------------------------------------------------------------ |
|       All 131 biosamples        | ![](sumstats_results/all_biosamples.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | *Same as column 3*                                           | ![](sumstats_results/all_biosamples.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|       Liver and Intestine       | ![](sumstats_results/liver_and_intestine.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | ![](sumstats_results/liver_and_intestine.all_putative_enhancers.merged_liver_and_intestine_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | ![](sumstats_results/liver_and_intestine.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png)) |
|              Liver              | ![](sumstats_results/liver.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | ![](sumstats_results/liver.all_putative_enhancers.merged_liver_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | ![](sumstats_results/liver.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|            Intestine            | ![](sumstats_results/intestine.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | ![](sumstats_results/intestine.all_putative_enhancers.merged_intestine_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | ![](sumstats_results/intestine.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|        hepatocyte-ENCODE        | ![](sumstats_results/hepatocyte-ENCODE.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/hepatocyte-ENCODE.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|          HepG2-Roadmap          | ![](sumstats_results/HepG2-Roadmap.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/HepG2-Roadmap.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|          liver-ENCODE           | ![](sumstats_results/liver-ENCODE.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/liver-ENCODE.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|  large_intestine_fetal-Roadmap  | ![](sumstats_results/large_intestine_fetal-Roadmap.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/large_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |
|  small_intestine_fetal-Roadmap  | ![](sumstats_results/small_intestine_fetal-Roadmap.all_putative_enhancers/refelt.scorequantile.nbconn.nbtimes.png) | *No overlap between enhancers of a single biosample. Same as column 1.* | ![](sumstats_results/small_intestine_fetal-Roadmap.all_putative_enhancers.merged_enhancers/refelt.scorequantile.nbconn.nbtimes.png) |



### Summary statistics on enhancer lengths

|           Enhancers            | 1st Qu. | Median | Mean | 3rd Qu. | Max |
| :-----------------------------: | :-----: | ------ | ---- | :-----: | ------- |
|       All 131 biosamples        | 200 | 308 | 482 | 637 |6991|
|       All 131 biosamples (merged)       | 288 | 631 | 808 | 1078 |11616|
|              Liver + intestine              | 200 | 353 | 491 | 644 |4626|
|            Liver + intestine (merged)            | 200 | 423 | 575 | 756 |6298|
|        Liver + intestine (merged 131)        | 520 | 982 | 1231 | 1668 |11616|
|          Liver          | 200 | 296 | 455 | 585 |4626|
|          Liver (merged)          | 200 | 348 | 507 | 650 |4737|
| Liver (merged 131) | 528 | 1027 | 1292 | 1780 |11616|
| Intestine | 200 | 428 | 543 | 724 |3679|
| Intestine (merged) |   221   | 488    | 607  |   814   | 4461  |
| Intestine (merged 131) | 710 | 1209 | 1469 | 1978 |11073|


