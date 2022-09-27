# Differences between results `example_chr22` and `test_chr22`

We compared our results on the small test example to the expected results using `meld`, and found that the results are identical except for least significants digits in some cases.

## In `ABC Output/Neighborhoods`

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

Changes list:

- ``ABC Output/Neighborhoods`EnhancerList.bed`: 4 minor differences (one-digit change/line) among 3330 lines
- `EnhancerList.txt`: a lot of minor changes (lists reordered) among 3331 lines
- `Enhancers.DHS.wgEncodeUwDnaseK562AlnRep1.chr22.bam.CountReads.bedgraph`: 4 minor changes (one-digit change/line) among 3277 lines
- `Enhancers.DHS.wgEncodeUwDnaseK562AlnRep2.chr22.bam.CountReads.bedgraph`: 4 minor changes (one-digit change/line) among 3330 lines
- `Enhancers.H3K27ac.ENCFF384ZZM.chr22.bam.CountReads.bedgraph`: 4 minor changes (one-digit change/line) among 3330 lines
- `GeneList.TSS1kb.bed`: 2*2 minor changes (lines with identical loci swapped)

The remaining 8 files are identical to their other version. Note that almost all digit-changes were simply a value incremented or decremented by 1 (and the others were in/decremented by 2).

## In `ABC Output/Peaks`

> ```shell
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

Changes list:

- `params.txt` is different from its original version but the differences arincremented be expected path differences
- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_model.r`: A lot of difference, but all of them are due to the increased number of significative digits in our test version

This fact is very likely to explain the observed one-digit changes in `ABC Output/Neighborhoods`.

- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak`: a lot of changes (generally 1, 2  or 3 significative-digits changes ; sometimes more)
- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted`: idem
- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.candidateRegions.bed`: 4 minor differences (one-digit change/line) among 3330 lines

Indeed, this explains all one-digit change in `ABC Output/Neighborhoods`.

- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.wgEncodeUwDnaseK562AlnRep1.chr22.bam.Counts.bed`: files are identical
- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.xls`: a lot of changes, same as `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak`. Another (expected) change is path in the header, and one comment deleted.
- `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_summits.bed`: a lot of changes, same as `wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak`.

## In `ABC_output/Predictions`

> ```bash
> $ l example_chr22/ABC_output/Predictions/
> EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz  EnhancerPredictionsFull.txt
> EnhancerPredictionsAllPutative.txt                      EnhancerPredictions.txt
> EnhancerPredictionsAllPutative.txt.gz                   GenePredictionStats.txt
> EnhancerPredictions.bedpe                               parameters.predict.txt
> ```

```bash
cd example_chr22/ABC_output/Predictions/
meld <name> ../../../test_chr22/ABC_output/Predictions/<name>
```

* `EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz`: a lot of minor changes (last or last few significant digits)
* `EnhancerPredictionsAllPutative.txt.gz`: a lot of minor changes (last or last few significant digits)
* `EnhancerPredictionsAllPutative.txt`: the file only exists in `example_chr22`, not in our test
* etc