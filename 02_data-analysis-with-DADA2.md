Dada2 tutorial
================

  - [On débute](#on-débute)
  - [On inspecte les profils de qualité de
    lecture](#on-inspecte-les-profils-de-qualité-de-lecture)
  - [On filtre et coupe](#on-filtre-et-coupe)
  - [On apprend les taux d’erreurs](#on-apprend-les-taux-derreurs)
  - [Inférence d’échantillon](#inférence-déchantillon)
  - [On lie les paires de reads](#on-lie-les-paires-de-reads)
  - [On construit la table
    d’observation](#on-construit-la-table-dobservation)
  - [On retire les chimères](#on-retire-les-chimères)
  - [Track reads through the
    pipeline](#track-reads-through-the-pipeline)
  - [On assigne la taxonomie](#on-assigne-la-taxonomie)
  - [Evaluate accuracy](#evaluate-accuracy)
  - [Bonus:handoff to phyloseq](#bonushandoff-to-phyloseq)
      - [Import dans phyloseq](#import-dans-phyloseq)
      - [Vizualisation de
        l’alpha-diversité](#vizualisation-de-lalpha-diversité)
      - [Ordination](#ordination)
      - [Histogramme](#histogramme)

# On débute

On charge dada2

``` r
library("dada2")
```

    ## Loading required package: Rcpp

On créé un chemin vers le répertoir MiSeq-SOP contenant les dossiers
fastq du séquençage.

``` r
path <- "~/MiSeq_SOP" 
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

On lit les fichiers fastq, dans ceux-ci les séquences Forward et
Reverses sont toutes ensembles, leurs noms ont un format spécifique.
Donc on met les séquences Forward à la variable fnFs, et les Reverses à
la variable fnRs. Pour se faire, avec la fonction “sort()” on cherche
tous les Forwards, qui ont "\_R1\_001.fastq" à la fin de leurs noms, et
on les affecte à fnFs. De même pour les Reverse, qui ont
"\_R2\_001.fastq" à la fin de leurs noms, et sont affectés à fnRs. Puis
on utilise la fonction “strsplit” pour découper le noms des séquences et
ne garder que leur identifiant. Et la fonction “sapply” permet
d’appliquer cela à toute le fichier. Ces noms sont affectés à la
variable “sample.names”.

``` r
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

# On inspecte les profils de qualité de lecture

On visualise les profils de qualité de lecture des reads Forward sous
forme de graphique avec la fonction “plotQualityProfile()”. Le “1:4”
prend les quatre premiers reads.

``` r
plotQualityProfile(fnFs[1:4])
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

On a donc les scores de qualité en ordonnée, et en abscisse la position
sur la séquence. En illumina on utilise 250 paires de base, c’est
pourquoi ça va jusqu’à 250. La fréquence est représentée en gris par une
Heatmap, plus haute est la fréquence, plus foncé est le gris. La ligne
verte montre la moyenne, la orange la médiane, et la orange pointillée
les premier et troisième quartiles.

On observe donc un bon score de qualité (\>30),à peu près constant,
excepté pour la dernière dizaine de nucléotides, où le score diminue
quelque peu. Il faudra peut-être couper ses derniers nucléotides pour
minimiser les erreurs.

# On filtre et coupe

Tout d’abord, on créé des variables pour les fichiers que l’on va
filtrer. Une variable pour les Forwards et une pour les Reverses. La
fonction “file.path()” indique le chemin à suivre pour que les séquences
soient assignées.

``` r
filtFs <- file.path(path,"filtered", paste0(sample.names, "_F_filt.fasta.gz"))
filtRs <- file.path(path,"filtered", paste0(sample.names, "_R_filt.fasta.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Ensuite, on utilise la fonction “filterAndTrim()” pour couper et faire
le tri. Elle est appliquée aux variables fnFs, filtFs, fnRs, et filtRs.
les reads Forward sont coupés au 240ème nucléotide, et les Reverses au
160ème nucléotide, leurs scores de qualités diminuant après. On applique
les paramètres de filtration standards, avec “maxEE” le nombre d’erreur
attendues maximum à 2. Tous ces fichiers filtrés sont assignés à out. On
montre ensuite le début de out avec la fonction “head”.

``` r
out <- filterAndTrim(fnFs,filtFs, fnRs, filtRs, truncLen = c(240,160),maxN = 0, maxEE = c(2,2),truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

# On apprend les taux d’erreurs

On va calculer le modèle de taux d’erreurs pour notre set de données
avec la fonction “learnErrors()”. Pour se faire, un taux maximal est
choisi initialement, et en fonction de l’échantillon il est modifié
jusqu’à atteindre un taux consistent. errF est la variable pour les
Forward.

``` r
errF <- learnErrors(filtFs, multithread = TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

Et errR, la variable pour les Reverses.

``` r
errR <- learnErrors(filtRs, multithread = TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

On vérifie visuellement :

``` r
plotErrors(errF, nominalQ = TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Ici, la probabilité de mutation d’une base à une autre est observée en
fonction du score de qualité. La ligne noir représente le taux d’erreurs
calculé précédement, la ligne rouge le taux d’erreurs initial d’après le
score de qualité, et les points sont les erreurs observées.

On observe que les taux d’erreurs calculés (lignes noires) correspondent
bien aux erreurs observées (points). De plus le taux d’erreur diminue
quand le score de qualité augmente. Donc on peut continuer.

# Inférence d’échantillon

Avec la fonction “dada()”, les séquences abondantes sont gardées pour
les reads, et les séquences en-dessous d’un seuil d’abondance, déterminé
par le taux d’erreur, ne sont pas gardées. Les séquences Forward gardées
sont affectées à la variable dadaFs.

``` r
dadaFs <- dada(filtFs, err=errF, multithread= TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

Les séquences Reverse gardées sont affectées à la variable dadaRs.

``` r
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

On appelle le premier objet de dadaFs. C’est le “dada-class”, qui montre
le résultat du débruitage effectué précédement. Ici, le débruitage a
donné 128 variants de séquence depuis 1979 séquences uniques dans le
premier échantillon.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# On lie les paires de reads

Avec la fonction “mergePairs()”, on aligne les reads Forward et Reverse,
qui ont été débruité, pour faire des contigs. Ces contigs sont assignés
à la variable mergers, dont on montre le début avec la fonction
“head()”.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
```

    ## 6551 paired-reads (in 106 unique pairings) successfully merged out of 6907 (in 199 pairings) input.

    ## 5025 paired-reads (in 100 unique pairings) successfully merged out of 5188 (in 156 pairings) input.

    ## 4973 paired-reads (in 80 unique pairings) successfully merged out of 5268 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2756 (in 109 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3622 paired-reads (in 53 unique pairings) successfully merged out of 4103 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6515 (in 198 pairings) input.

    ## 3961 paired-reads (in 90 unique pairings) successfully merged out of 4384 (in 188 pairings) input.

    ## 14231 paired-reads (in 143 unique pairings) successfully merged out of 15358 (in 351 pairings) input.

    ## 10526 paired-reads (in 120 unique pairings) successfully merged out of 11166 (in 279 pairings) input.

    ## 11156 paired-reads (in 137 unique pairings) successfully merged out of 11799 (in 298 pairings) input.

    ## 4329 paired-reads (in 84 unique pairings) successfully merged out of 4788 (in 180 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7193 (in 187 pairings) input.

    ## 4430 paired-reads (in 67 unique pairings) successfully merged out of 4605 (in 127 pairings) input.

    ## 4574 paired-reads (in 100 unique pairings) successfully merged out of 4736 (in 172 pairings) input.

    ## 6094 paired-reads (in 109 unique pairings) successfully merged out of 6314 (in 172 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

# On construit la table d’observation

On utilise la fonction “makeSequenceTable()” pour construire un tableau
de séquence variant d’amplicon, à partir des contigs formé plutôt. Ce
qui est plus précis qu’un tableau d’OTU. Le tableau est dans la variable
seqtab. Et on utilise “dim()” pour avoir la dimension du tableau.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

On observe la distribution de la longueur des séquences avec la fonction
“nchar()”.

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

``` r
#inspect distribution os sequences lengths
```

La matrice contient 293 ASVs, et la longueur de séquences fusionnées est
entre 252 et 253, ce qui était attendu.

# On retire les chimères

Malgré le fait qu’on ai filtré les erreurs, des chimères peuvent rester.
En utilisant la fonction “remove BimeraDenovo()” avec la méthode
“consensus”, les chimères sont détectées quand il y a des bimeras.
C’est à dire quand en assemblant les séquences les plus abondantes, on
peut reconstuire à l’identique la séquence chimérique.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread= TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

61 bimeras ont été identifiées, il reste 20 232 séquences.

On calcul la fréquence de séquence chimérique.

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.964263

Celle-ci est de O.964, ce qui me semble élevé.

# Track reads through the pipeline

On regarde le nombre de reads après chaque étape.

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6996      6978   6551    6539
    ## F3D1    5869     5299      5227      5239   5025    5014
    ## F3D141  5958     5463      5339      5351   4973    4850
    ## F3D142  3183     2914      2799      2833   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4146      4224   3622    3483

Il n’y a pas eu trop de perte.C’est bon signe.

# On assigne la taxonomie

On charge les données silva.

Puis on utilise la fonction “assignTaxonomy()” pour assigner la
taxonomie du fichier téléchargé à la table de séquences sans chimère.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread = TRUE)
```

Et avec la fonction “addSpecies()” on assigne aussi les espèces aux ASVs
qui correspondent exactement aux séquences de référence.

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

On regarde le résultat de l’assignation.

``` r
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

Les Bacteroides sont dans le taxon le plus abondant. Beaucoup d’espèces
n’ont pas été assignées car il est difficile avec du 16s de ne pas avoir
des assignations ambigues, de plus le microbiote des intestin de souris
n’est pas assez bien référencé dans les bases de données.

# Evaluate accuracy

On va utiliser la communauté fictive pour vérifier les variants obtenus
avec DADA2. D’abord on enlève de la communauté fictive les séquences qui
ne sont pas présentes dans les données brutes avec la fonction “sort()”.
Puis on concatène avec la fonction “cat()”.

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

Maintenant, on compare les deux sets de séquences.

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

DADA2 a donc identifié 20 séquences qui correspondent à celles du génome
de référence, et le taux d’erreur pour cet échantillon est 0%. Bien sûr
celà n’est pas représentatif d’un cas réel, le taux d’erreur sera plus
élevé et toutes les séquences ne corresponderont pas au génome
référence.

# Bonus:handoff to phyloseq

On installe tous les packages nécéssaires.

``` r
BiocManager::install("phyloseq")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phyloseq'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.32.0'

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## [1] '2.56.0'

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.2'

## Import dans phyloseq

Après avoir chargé les packages nécéssaires, on applique le thème “bw”.

``` r
theme_set(theme_bw())
```

Puis on construit un échantillon avec la fonction “data.frame()”.
Normalement on utilise les données d’un fichier.

``` r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

Avec la fonction “phyloseq()” on construit l’objet “ps” qui va contenir
la table d’OTU, les données des échantillons, et le tableau taxonomique.
Puis on enlève les échantillons “Mocks” pour simuler un filtrage avec la
fonction “prune\_samples()”.

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) 
```

On créé des diminutifs aux noms des séquences en liant les noms de
taxons aux séquences avec la fonction “merge\_phyloseq()”, puis en
raccouricissant ces noms.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 232 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 232 reference sequences ]

## Vizualisation de l’alpha-diversité

On trace les graphiques d’alpha-diversité, avec des indices de Shannon
et Simpson, en utilisant la fonction “plot\_richness()”.

``` r
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Les deux indices ne montrent pas de différence notable entre les
échantillons à 0 jour et à 150 jours.

## Ordination

On transforme les données à des proportions adaptées aux distances de
Bray-Curtis avec la fonction “transform\_sample\_counts()”. Ensuite on
fait une ordination NMDS avec la fonction “ordinate()”.

``` r
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.08574537 
    ## Run 1 stress 0.374426 
    ## Run 2 stress 0.08574537 
    ## ... New best solution
    ## ... Procrustes: rmse 2.253804e-06  max resid 5.71696e-06 
    ## ... Similar to previous best
    ## Run 3 stress 0.08574537 
    ## ... New best solution
    ## ... Procrustes: rmse 1.19247e-06  max resid 2.593456e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.1233096 
    ## Run 5 stress 0.08574537 
    ## ... Procrustes: rmse 3.14537e-06  max resid 1.031066e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04283499  max resid 0.1433426 
    ## Run 7 stress 0.08574537 
    ## Run 8 stress 0.08574537 
    ## Run 9 stress 0.08942874 
    ## Run 10 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 6.048282e-06  max resid 1.836699e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.3554612 
    ## Run 12 stress 0.1216669 
    ## Run 13 stress 0.08574537 
    ## Run 14 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 3.15448e-06  max resid 8.600243e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.144585 
    ## Run 16 stress 0.09421601 
    ## Run 17 stress 0.08942867 
    ## Run 18 stress 0.08942865 
    ## Run 19 stress 0.126472 
    ## Run 20 stress 0.09421601 
    ## *** Solution reached

Puis on trace le graphique avec la fonction “plot\_ordination()”.

``` r
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

On observe un nette séparation entre les échantillons tôts et tards.

## Histogramme

On trace l’histogramme des familles en fonction de l’abondance et des
jours avec la fonction “plot\_bar()”.

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

![](02_data-analysis-with-DADA2_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

On ne peut rien observer dans la distribution taxonomique des 20
premières séquences qui pourrait expliquer le résultat de l’ordination.
