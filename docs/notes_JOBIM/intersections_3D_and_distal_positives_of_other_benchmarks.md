# Intersection between positive e-g pairs from 3D-based dataset and distal positive e-g pairs from other datasets



## Annotation des gènes

### GENCODE

#### Tous les gènes

Nombre de gènes uniques (**pas** le nombre de transcrits) présents dans l'annotation :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} $3=="gene"' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 57,820

Nombre total de transcrits de gènes dans l'annotation :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} $3=="transcript"' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 196,520

#### Gènes codants

Nombre de gènes codants uniques dans l'annotation :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="gene"){split($9,parts,"\""); if(parts[6]~/(^protein_coding$)/){print $0}}}' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 20,345

Nombre de transcrits de gènes codants dans l'annotation :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="transcript"){split($9,parts,"\""); if(parts[6]~/(^protein_coding$)/){print $0}}}' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 145,641

### Moore et al

> ```bash
> wc -l ../Annotations/GENCODEv19-TSSs.bed
> ```
>
> 167,147

167 147 transcrits distincts apparaissent dans la liste de transcrits fournie par Moore et al.

### Intersection avec la liste des gènes trouvés dans les BENGI

#### Intersection avec les gènes associés aux TSS fournis par Moore et al

Nombre de gènes uniques dans l'annotation des transcrits/TSS fournie par Moore et al :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {genes[$7]++} END{for(u in genes){print u, genes[u]}}' ../Annotations/GENCODEv19-TSSs.bed |wc -l
> ```
>
> 54,846

Ces gènes sont-ils dans notre annotation Gencode ?

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){if($3=="gene"){split($9,parts,"\""); genes[parts[2]]++}; next}; if(genes[$7]){found[$7]++}} END{for(u in found){print u, found[u]}}' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf ../Annotations/GENCODEv19-TSSs.bed |wc -l
> ```
>
> 54,846

Fortunately all genes of which the TSS are provided by Moore et al, are included in our gene annotation! So our annotation is indeed the correct one. The difference must be due to Moore et al having performed some filters.

A l'inverse, regardons à quoi ressemblent les 2974 gènes qui sont dans notre annotation mais pas dans celle obtenue via la liste de transcrits/TSS donnée par Moore et al :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$7]++; next}; if($3=="gene"){split($9,parts,"\""); if(!genes[parts[2]]){not_found[parts[2]]++}}} END{for(u in not_found){print u, not_found[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 2974

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$7]++; next}; if($3=="gene"){split($9,parts,"\""); if(!genes[parts[2]]){types[parts[6]]++}}} END{for(u in types){print u, types[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf
> ```
>
> antisense	2
> TR_J_pseudogene	4
> IG_C_pseudogene	9
> TR_V_pseudogene	27
> processed_transcript	23
> polymorphic_pseudogene	19
> lincRNA	1
> pseudogene	2699
> IG_J_pseudogene	3
> IG_V_pseudogene	187

Dans les 54 846 gènes de Moore et al, on retrouve des gènes de chacun de ces types, qui n'expliquent donc pas les 2974 manquants.

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$7]++; next}; if($3=="gene"){split($9,parts,"\""); if(genes[parts[2]]){types[parts[6]]++}}} END{for(u in types){print u, types[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf
> ```
>
> 3prime_overlapping_ncrna	21
> TR_D_gene	3
> misc_RNA	2034
> Mt_rRNA	2
> IG_J_gene	18
> TR_C_gene	5
> rRNA	527
> antisense	5274
> sense_overlapping	202
> Mt_tRNA	22
> TR_J_gene	74
> IG_V_gene	138
> polymorphic_pseudogene	26
> processed_transcript	492
> lincRNA	7113
> pseudogene	11232
> IG_D_gene	37
> miRNA	3055
> sense_intronic	742
> snoRNA	1457
> protein_coding	20345
> TR_V_gene	97
> IG_C_gene	14
> snRNA	1916



Remarquons qu'on obtient exactement la même chose en s'intéressant à la colonne "type des transcrits" pour ces gènes (ie en remplaçant par types[parts[6]]) (logique : pour les lignes "gènes", les colonnes "transcrits" contiennent la même info que les colonnes "gènes"). Par contre, puisqu'on va voir dans la suite que même en se restreignant, dans notre annotation Gencode, aux 54 846 gènes intervenant dans la liste des TSS de Moore et al, afin de reconstruire la liste des TSS, on trouve encore une liste de TSS plus longue que celle proposée par Moore et al. Partant de là, on va regarder quels sont les types des transcrits des TSS proposés par Moore et al ; et les transcrits évincés correspondent sont tous ceux d'un ou plusieurs types particuliers.

#### Intersection avec les gènes présents dans les BENGI



## Annotation des TSS

### A partir de l'annotation `fragencode`

#### Tous les gènes

Nombre total de TSS de gènes dans l'annotation :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="transcript"){$7=="+"? TSS[$4]++:TSS[$5]++}} END{for(u in TSS){print u, TSS[u]}}' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 178,758

#### Gènes codants

Nombre de TSS de gènes codants distincts :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if($3=="transcript"){split($9,parts,"\""); if(parts[6]~/(^protein_coding$)/){$7=="+"? TSS[$4]++:TSS[$5]++}}} END{for(u in TSS){print u, TSS[u]}}' /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 129,785

#### Gènes qui apparaissent dans la liste des TSS fournie par Moore et al

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$7]++; next}; split($9,parts,"\""); if(genes[parts[2]]){if($3=="transcript"){$7=="+"? TSS[$4]++:TSS[$5]++}}} END{for(u in TSS){print u, TSS[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 175,554

Remarquons que 178 758 - 175 554 = 3220 n'est pas très éloigné du nombre de gènes de notre annotation qui ne sont pas inclus dans celle de Moore et al - à savoir 2974. Ces gènes non inclus dans la liste de TSS de Moore et al sont donc des gènes qui ont peu de transcrits différents.

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){genes[$7]++; next}; split($9,parts,"\""); if(!genes[parts[2]]){if($3=="transcript"){$7=="+"? TSS[$4]++:TSS[$5]++}}} END{for(u in TSS){print u, TSS[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 3,220

### Fournie par Moore et al

Nombre de transcrits distincts :

> ```bash
> wc -l /work2/project/regenet/results/multi/bengi/Annotations/GENCODEv19-TSSs.bed
> ```
>
> 167,147

Nombre de TSS distincts :

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {$6=="+"?TSS[$2]++:TSS[$3]++} END{for(u in TSS){print u, TSS[u]}}' ../Annotations/GENCODEv19-TSSs.bed |wc -l
> ```
>
> 151,722

On a donc 175,554 - 151,722 = 23,383 TSS de moins. Reste à voir si les autres sont les mêmes, et qu'est-ce qui a pu motiver l'exclusion des TSS manquants.

### Conclusion

* 57 820 gènes (codants + non codants) dans notre annotation, et 54 846 (inclus dans nos 57 820) dans celle obtenue à partir de la liste des transcrits/TSS fournie par Moore et al.
* Il n'est pas évident de trouver quel filtre a été appliqué pour conduire à avoir 2 974 gènes de moins.
* Le plus étonnant est qu'en utilisant notre annotation `fragencode` pour construire la liste des TSS, mais en se restreignant aux 54 846 gènes inclus dans celle donnée par Moore et al., on trouve quand même davantage de TSSs : 175 554 au lieu de 151 722, soit 23 383 de plus

### Les TSS évincés correspondent-ils à tous les transcrits de certains types ?

#### Quelques chiffres

- Nombre de TSS dans notre annotation pour les 57 820 gènes : 178,758 (pour 196 520 transcrits)
- Nombre de TSS dans notre annotation pour les 54 846 / 57 820 gènes : 175 554 (pour 193 245 transcrits)

57 820 - 54 846 = 2 974 gènes de moins. 196 520 - 193 245 = 3 275 transcrits de moins. 178,758 - 175 554 = 3,204 TSS de moins.

- 193 245 - 167 147 = 26 098 transcrits de moins dans notre annotation réduite que dans celle de Moore et al.
- 175,554 - 151,722 = 23,383 TSS de moins entre notre annotation réduite et la liste des TSS de Moore et al.
- Il y a en fait exactement 26 098 transcrits, parmi les 193 245 transcrits de notre annotation réduite, qui ne sont pas dans la liste de Moore et al (voir ci-dessous). Donc la liste de transcrits de Moore et al. est incluse dans notre liste réduite de transcrits, pour les 54 846 gènes.

Intéressons nous donc au type des 26 098 transcrits évincés de l'annotation de Moore et al.



> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){transcripts[$4]++; genes[$7]++; next}; if($3=="transcript"){split($9,parts,"\""); if(genes[parts[2]]){print $0}}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 193,245
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){transcripts[$4]++; genes[$7]++; next}; if($3=="transcript"){split($9,parts,"\""); if(genes[parts[2]]){if(!transcripts[parts[4]]){print $0}}}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf |wc -l
> ```
>
> 26,098
>
> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){transcripts[$4]++; genes[$7]++; next}; if($3=="transcript"){split($9,parts,"\""); if(genes[parts[2]]){if(!transcripts[parts[4]]){types[parts[12]]++}}}} END{for(u in types){print u, types[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf
> ```
>
> unitary_pseudogene	76
> pseudogene	63
> polymorphic_pseudogene	39
> unprocessed_pseudogene	3
> retained_intron	25917

Ces types de transcrits apparaissent-ils aussi dans des transcrits qui ne sont pas évincés ? Regardons.

> ```bash
> awk 'BEGIN{FS="\t"; OFS="\t"} {if(NR==FNR){transcripts[$4]++; genes[$7]++; next}; if($3=="transcript"){split($9,parts,"\""); if(genes[parts[2]]){if(transcripts[parts[4]]){types[parts[12]]++}}}} END{for(u in types){print u, types[u]}}' ../Annotations/GENCODEv19-TSSs.bed /work2/project/fragencode/data/species.bk/homo_sapiens/hg19.gencv19/homo_sapiens.gtf
> ```
>
> antisense	9710
> TR_V_gene	97
> transcribed_unprocessed_pseudogene	860
> protein_coding	81814
> sense_intronic	802
> TR_J_gene	74
> non_stop_decay	58
> rRNA	531
> TR_C_gene	5
> transcribed_processed_pseudogene	442
> IG_D_gene	37
> sense_overlapping	330
> miRNA	3116
> TR_D_gene	3
> processed_transcript	28082
> 3prime_overlapping_ncrna	25
> Mt_tRNA	22
> misc_RNA	2050
> Mt_rRNA	2
> snRNA	1923
> IG_C_gene	18
> IG_J_gene	18
> snoRNA	1529
> IG_V_gene	144
> nonsense_mediated_decay	13052
> lincRNA	11780
> processed_pseudogene	10623



Parfait ! C'est le mieux qu'on pouvait espérer. On ne sait pas quels filtres Moore et al ont appliqué pour exclure 2 974 gènes, mais on sait que pour les gènes restants, ils ont exactement exclu tous les transcrits de types "retained_intron", "unitary_pseudogene", "pseudogene", "polymorphic_pseudogene", "unprocessed_pseudogene".

## Conversion des jeux de données BENGI au format `bedpe`

Il y a essentiellement trois stratégies possibles :

* soit en utilisant exclusivement l'annotation utilisée par Moore et al, c'est-à-dire leur liste de 151 722 TSS uniques (167 147 transcrits uniques) et notre annotation pour les 54 846 gènes associés 
* soit en utilisant nos 175 554 TSS (pour 193 245 transcrits) associés aux 54 846 gènes de Moore et al
* soit en utilisant notre annotation complète : 57 820 gènes, 178 758 TSS (pour 196 520 transcrits)

### Avec l'annotation de Moore et al

Déjà fait.

### Avec l'annotation GENCODE complète

A faire.



## Choix de la distance minimum pour considérer les élements comme distaux

On propose 3 stratégies :

* distance min. des paires HiC
* a-percentile des distances des paires HiC, avec a=1 par exemple
* a-percentile des distances des paires CHiC, avec a=5 ou 10 par exemple



