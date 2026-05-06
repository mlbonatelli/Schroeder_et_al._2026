Paper_comparison
================
2026-04

# Loading libraries. If needed, use BiocManager to install them

    ## Package: phyloseq
    ## Version: 1.48.0
    ## 
    ## -- File: C:/Users/maria/AppData/Local/R/win-library/4.4/phyloseq/Meta/package.rds 
    ## -- Fields read: Package, Version

    ## Package: ggplot2
    ## Version: 4.0.0
    ## 
    ## -- File: C:/Users/maria/AppData/Local/R/win-library/4.4/ggplot2/Meta/package.rds 
    ## -- Fields read: Package, Version

    ## Package: vegan
    ## Version: 2.7-1
    ## 
    ## -- File: C:/Users/maria/AppData/Local/R/win-library/4.4/vegan/Meta/package.rds 
    ## -- Fields read: Package, Version

    ## Package: readr
    ## Version: 2.2.0
    ## 
    ## -- File: C:/Users/maria/AppData/Local/R/win-library/4.4/readr/Meta/package.rds 
    ## -- Fields read: Package, Version

    ## Package: tidyverse
    ## Version: 2.0.0
    ## 
    ## -- File: C:/Users/maria/AppData/Local/R/win-library/4.4/tidyverse/Meta/package.rds 
    ## -- Fields read: Package, Version

    ## Package: RColorBrewer
    ## Version: 1.1-3
    ## 
    ## -- File: C:/Users/maria/AppData/Local/R/win-library/4.4/RColorBrewer/Meta/package.rds 
    ## -- Fields read: Package, Version

# Importing files, renaming ASVs and creating phyloseq object - AUER

``` r
seqtab.nochim_auer <- readr::read_csv("Auer_seqtab.nochim_16S.csv")
```

    ## New names:
    ## Rows: 16 Columns: 4546
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (4545):
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTT...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
seqtab.nochim_auer <- as.data.frame(seqtab.nochim_auer)
rownames(seqtab.nochim_auer) <- seqtab.nochim_auer[, 1]
seqtab.nochim_auer <- seqtab.nochim_auer[, -1]
colnames(seqtab.nochim_auer) = paste0("ASV", 1:ncol(seqtab.nochim_auer)) 

taxa_auer <- read.csv("Auer_taxa_16S.csv", row.names = 1)
row.names(taxa_auer) = colnames(seqtab.nochim_auer)
taxa_auer <- as.matrix(taxa_auer)

metadata_auer <- read.csv("Metadata_AUER.csv", stringsAsFactors = TRUE, row.names = 1)

ps_auer <- phyloseq(
  otu_table(seqtab.nochim_auer, taxa_are_rows=FALSE),
  sample_data(metadata_auer),
  tax_table(taxa_auer))

ps_auer
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 4545 taxa and 16 samples ]
    ## sample_data() Sample Data:       [ 16 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 4545 taxa by 6 taxonomic ranks ]

# Importing files, renaming ASVs and creating phyloseq object - DENG_2017

``` r
seqtab.nochim_DENG7 <- readr::read_csv("Deng_seqtab.nochim_16S.csv")
```

    ## New names:
    ## Rows: 12 Columns: 409
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (408):
    ## TCGAGAATAGTCTACAATGGACGGAAGTCTGATAGTGCGACGCCGCGTGAACGAAGAATCCCTTC...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
seqtab.nochim_DENG7 <- as.data.frame(seqtab.nochim_DENG7)
rownames(seqtab.nochim_DENG7) <- seqtab.nochim_DENG7[, 1]
seqtab.nochim_DENG7 <- seqtab.nochim_DENG7[, -1]
colnames(seqtab.nochim_DENG7) = paste0("ASV", 1:ncol(seqtab.nochim_DENG7)) 

taxa_DENG7 <- read.csv("Deng_taxa_16S.csv", row.names = 1)
row.names(taxa_DENG7) = colnames(seqtab.nochim_DENG7)
taxa_DENG7 <- as.matrix(taxa_DENG7)

metadata_DENG7 <- read.csv("Metadata_DENG_2017.csv", stringsAsFactors = TRUE, row.names = 1)

ps_DENG7 <- phyloseq(
  otu_table(seqtab.nochim_DENG7, taxa_are_rows=FALSE),
  sample_data(metadata_DENG7),
  tax_table(taxa_DENG7))

ps_DENG7
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 408 taxa and 12 samples ]
    ## sample_data() Sample Data:       [ 12 samples by 6 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 408 taxa by 6 taxonomic ranks ]

# Importing files, renaming ASVs and creating phyloseq object - DENG_2018

``` r
seqtab.nochim_DENG8 <- readr::read_csv("Deng_2_seqtab.nochim_16S.csv")
```

    ## New names:
    ## Rows: 10 Columns: 304
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (303):
    ## TGAGGAATATTGGACAATGGCCGGAAGGCTGATCCAGCCATGCCGCGTGCGGGAGGACGGCCCTA...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
seqtab.nochim_DENG8 <- as.data.frame(seqtab.nochim_DENG8)
rownames(seqtab.nochim_DENG8) <- seqtab.nochim_DENG8[, 1]
seqtab.nochim_DENG8 <- seqtab.nochim_DENG8[, -1]
colnames(seqtab.nochim_DENG8) = paste0("ASV", 1:ncol(seqtab.nochim_DENG8)) 

taxa_DENG8 <- read.csv("Deng_2_taxa_16S.csv", row.names = 1)
row.names(taxa_DENG8) = colnames(seqtab.nochim_DENG8)
taxa_DENG8 <- as.matrix(taxa_DENG8)

metadata_DENG8 <- read.csv("Metadata_DENG_2018.csv", stringsAsFactors = TRUE, row.names = 1)

ps_DENG8 <- phyloseq(
  otu_table(seqtab.nochim_DENG8, taxa_are_rows=FALSE),
  sample_data(metadata_DENG8),
  tax_table(taxa_DENG8))

ps_DENG8
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 303 taxa and 10 samples ]
    ## sample_data() Sample Data:       [ 10 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 303 taxa by 6 taxonomic ranks ]

# Importing files, renaming ASVs and creating phyloseq object - FENG

``` r
seqtab.nochim_FENG <- readr::read_csv("Feng_seqtab.nochim_16S.csv")
```

    ## New names:
    ## Rows: 18 Columns: 954
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (953):
    ## TGGGGAATATTGCGCAATGGGGGGAACCCTGACGCAGCGACGCCGCGTGAAGGAAGAAGGCCTTC...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
seqtab.nochim_FENG <- as.data.frame(seqtab.nochim_FENG)
rownames(seqtab.nochim_FENG) <- seqtab.nochim_FENG[, 1]
seqtab.nochim_FENG <- seqtab.nochim_FENG[, -1]
colnames(seqtab.nochim_FENG) = paste0("ASV", 1:ncol(seqtab.nochim_FENG)) 

taxa_FENG <- read.csv("Feng_taxa_16S.csv", row.names = 1)
row.names(taxa_FENG) = colnames(seqtab.nochim_FENG)
taxa_FENG <- as.matrix(taxa_FENG)

metadata_FENG <- read.csv("Metadata_FENG.csv", stringsAsFactors = TRUE, row.names = 1)

ps_FENG <- phyloseq(
  otu_table(seqtab.nochim_FENG, taxa_are_rows=FALSE),
  sample_data(metadata_FENG),
  tax_table(taxa_FENG))

ps_FENG
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 953 taxa and 18 samples ]
    ## sample_data() Sample Data:       [ 18 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 953 taxa by 6 taxonomic ranks ]

# Importing files, renaming ASVs and creating phyloseq object - Current work

``` r
seqtab.nochim <- readr::read_csv("seqtab.nochim.csv")
```

    ## New names:
    ## Rows: 95 Columns: 1587
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (1586):
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTT...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
seqtab.nochim <- as.data.frame(seqtab.nochim)
rownames(seqtab.nochim) <- seqtab.nochim[, 1]
seqtab.nochim <- seqtab.nochim[, -1]
colnames(seqtab.nochim) = paste0("ASV", 1:ncol(seqtab.nochim)) 

taxa <- read.csv("taxa.csv", row.names = 1)
row.names(taxa) = colnames(seqtab.nochim)
taxa <- as.matrix(taxa)

metadata <- read.csv("metadata.csv", stringsAsFactors = TRUE, row.names = 1)

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  sample_data(metadata),
  tax_table(taxa))

# excluding foam samples that were not used in this work
ps_Without_FOAM <- subset_samples(ps, Phase != "R1_Foam")
ps_Without_FOAM <- subset_samples(ps_Without_FOAM, Phase != "R2_Foam")
ps_Without_FOAM <- subset_samples(ps_Without_FOAM, Phase != "R3_Foam")
ps_Without_FOAM
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1586 taxa and 92 samples ]
    ## sample_data() Sample Data:       [ 92 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1586 taxa by 6 taxonomic ranks ]

# Joining phyloseq objects

``` r
ps_all <- merge_phyloseq(ps_Without_FOAM, ps_auer, ps_DENG7, ps_DENG8, ps_FENG)
ps_all
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 4545 taxa and 148 samples ]
    ## sample_data() Sample Data:       [ 148 samples by 24 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 4545 taxa by 6 taxonomic ranks ]

# Rarefying

``` r
set.seed(1)

counts <- sample_sums(ps_all)
sort(counts)
```

    ## SRR6151783 SRR6151787 SRR6151791         B6 SRR6151782 SRR6151793 SRR6151789 
    ##        391       2867       8941       9421      11738      12193      13390 
    ## SRR5028698 ERR3445453 SRR6151785 SRR5028318 SRR5028641 SRR6151780 SRR5028701 
    ##      15660      16915      18557      18729      19987      20978      24219 
    ## SRR5028640 ERR3445450 ERR3445456 SRR5028307 SRR5028694 ERR3445459 ERR3445454 
    ##      24856      25202      25785      25969      26097      26147      26666 
    ## ERR3445444 SRR5028696 ERR3445462 ERR3445463 SRR6122451 ERR3445465 SRR5028645 
    ##      27102      27505      28236      28250      28286      29637      29839 
    ## ERR3445447 SRR6122513 SRR5028692 SRR5028658 ERR3445451 ERR3445457 ERR3445445 
    ##      31024      31971      33480      34013      35190      35325      35544 
    ## SRR6122514 SRR6122511 ERR3445469 SRR6151781 SRR5028661 ERR3445448 ERR3445468 
    ##      38037      38252      38857      40260      40338      40390      40728 
    ##        B66 SRR6122509 SRR6151792        B82        B17 SRR6151779 ERR3445460 
    ##      41839      41965      43589      43613      44403      48218      48355 
    ## SRR6151794 SRR6151788        B41 ERR3445466 SRR6122510         B2        B55 
    ##      49120      49892      50224      50722      52148      53501      53925 
    ##        B33 SRR6151786        B84 SRR6151790        B74        B27        B49 
    ##      54140      54141      55187      56048      56729      57028      57145 
    ## SRR6122507        B81 SRR6151784        B73         B3        B51         B4 
    ##      57655      58060      59722      61032      61150      61272      61952 
    ##        B32        B57        B76        B25         B7 SRR6122506        B23 
    ##      63046      63408      65244      65629      65640      65761      65791 
    ##        B89        B52        B60        B39        B91        B34        B86 
    ##      66136      66816      67461      68218      68727      68746      69202 
    ##        B31        B80        B22         B5        B37        B19        B24 
    ##      69711      69749      70397      70742      71139      71155      71412 
    ##        B67        B90        B13        B88        B83        B64        B30 
    ##      72099      72403      72833      73423      73682      73728      74551 
    ##        B46        B70        B12        B68        B56        B18        B72 
    ##      74734      74994      75072      75101      75160      75237      75603 
    ##        B92        B58 SRR6122508        B28        B26        B79        B21 
    ##      75836      75840      76102      76181      76194      76856      77026 
    ##        B54        B78        B43        B87         B8         B1        B16 
    ##      77377      78910      79545      80355      80545      80610      80610 
    ##         B9        B11        B15        B40 SRR6122515        B29        B48 
    ##      81230      81925      83416      84621      85526      86536      86658 
    ##        B75        B44        B62        B71        B59        B63        B50 
    ##      87107      88618      89237      89894      90450      91672      92145 
    ##        B10        B47        B35        B14        B65        B61        B53 
    ##      92546      93026      95549      97532     116686     191399     213724 
    ##        B77        B20        B69        B36        B45        B42        B85 
    ##     232114     241982     243451     262624     268189     393335     418839 
    ##        B93 
    ##     432096

``` r
ps_all_rar<- rarefy_even_depth(ps_all, sample.size = 11738,
                  rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 4 samples removedbecause they contained fewer reads than `sample.size`.

    ## Up to first five removed samples are:

    ## B6SRR6151783SRR6151787SRR6151791

    ## ...

    ## 756OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps_all_rar
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3789 taxa and 144 samples ]
    ## sample_data() Sample Data:       [ 144 samples by 24 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3789 taxa by 6 taxonomic ranks ]

# Permanova

``` r
dist <- phyloseq::distance(ps_all_rar, method="bray")
head(dist)
```

    ## [1] 0.7098313 0.7278071 0.7285739 0.7084682 0.7085534 0.6758392

``` r
metadata_dist <- data.frame(sample_data(ps_all_rar))
head(metadata_dist)
```

    ##     geo_loc_name_country_continent     Paper Sample_Name geo_loc_name_country
    ## B1                          Europe This work        <NA>              Germany
    ## B10                         Europe This work        <NA>              Germany
    ## B11                         Europe This work        <NA>              Germany
    ## B12                         Europe This work        <NA>              Germany
    ## B13                         Europe This work        <NA>              Germany
    ## B14                         Europe This work        <NA>              Germany
    ##     Sample.Name isolation_source Reactor_tax Reactor_div Reactor_nmd
    ## B1         <NA>             <NA>          R1          R1          R1
    ## B10        <NA>             <NA>          R1          R1          R1
    ## B11        <NA>             <NA>          R1          R1          R1
    ## B12        <NA>             <NA>          R1          R1          R1
    ## B13        <NA>             <NA>          R1          R1          R1
    ## B14        <NA>             <NA>          R1          R1          R1
    ##     Sampling_day Sampling_point Days_operation OLR HRT Sampling_interval Code
    ## B1      13.11.20             T1              2   5   8                 4   B1
    ## B10     14.01.21            T32             64   5   8                12  B10
    ## B11     26.01.21            T38             76   5   8                12  B11
    ## B12     07.02.21            T44             88   5   8                12  B12
    ## B13     19.02.21            T50            100   5   8                12  B13
    ## B14     03.03.21            T56            112   5   8                12  B14
    ##           Phase Library_Name Host Source Organism rel_to_oxygen Library.Name
    ## B1  Adapt_phase         <NA> <NA>   <NA>     <NA>          <NA>         <NA>
    ## B10      Cond_1         <NA> <NA>   <NA>     <NA>          <NA>         <NA>
    ## B11      Cond_1         <NA> <NA>   <NA>     <NA>          <NA>         <NA>
    ## B12      Cond_1         <NA> <NA>   <NA>     <NA>          <NA>         <NA>
    ## B13      Cond_1         <NA> <NA>   <NA>     <NA>          <NA>         <NA>
    ## B14      Cond_1         <NA> <NA>   <NA>     <NA>          <NA>         <NA>
    ##     geographic_location_.country_and.or_sea.
    ## B1                                      <NA>
    ## B10                                     <NA>
    ## B11                                     <NA>
    ## B12                                     <NA>
    ## B13                                     <NA>
    ## B14                                     <NA>

``` r
test.adonis_studies <- adonis2(dist ~ Paper, data = metadata_dist)
test.adonis_studies
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dist ~ Paper, data = metadata_dist)
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## Model      4    4.541 0.13943 5.6303  0.001 ***
    ## Residual 139   28.029 0.86057                  
    ## Total    143   32.571 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Taxonomy plot - Family

``` r
custom_col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7")

y1p <- tax_glom(ps_all_rar, taxrank = 'Family') # agglomerate taxa
y3p <- transform_sample_counts(y1p, function(x) x/sum(x)) #get abundance in %
y4p <- psmelt(y3p) # create dataframe from phyloseq object
y4p$Family <- as.character(y4p$Family) #convert to character
y4p$Family[y4p$Abundance < 0.05] <- " Family < 5% abund." #rename genera with < 1% abundance

#set color palette to accommodate the number of genera
colourCount <- length(unique(y4p$Family))
getPalette = colorRampPalette(custom_col42, interpolate="spline")

# all samples
#png("figures_article/Tax_Family_days.png", width = 10, height = 6, units = "in", res=300)
p_p <- ggplot(data=y4p, aes(x=Sample, y=Abundance, fill=Family))
p_p + geom_bar(aes(), stat="identity", position="stack") + 
  facet_wrap(~ Paper , scales = "free", labeller = labeller(.multi_line = FALSE)) +
  scale_fill_manual(values=getPalette(colourCount)) + 
  xlab("Samples") + ylab("Relative Abundance") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="bottom")
```

![](Paper_comparison_files/figure-gfm/tax%20plot%20family-1.png)<!-- -->

``` r
#dev.off()
```
