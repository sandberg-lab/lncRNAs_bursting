``` r
options(stringsAsFactors = FALSE)
`%notin%` <- Negate(`%in%`)

# Basic folders
  home.dir <- '/home/perj/'

# Project path
  prj_name <- 'prj_lncRNAs'
  prj_path = file.path(home.dir, 'projects', prj_name)

# Data
  rpkm_path = file.path(prj_path, 'data/ss2_n751_RPKM.rds')
  counts_path = file.path(prj_path, 'data/ss2_n751_readCounts.rds')
  cast_path <- file.path(prj_path, 'data/ss2_n751_castCounts.rds')
  c57_path <- file.path(prj_path, 'data/ss2_n751_c57Counts.rds')
  meta_path = file.path(prj_path, 'data/ss2_n751_meta.rds')
  biomart_path = file.path(prj_path, 'data/ss2_n533_biomart.rds')
    
# Resources
  imprint_path = file.path(home.dir, 'resources/GeneImprint/mouse.csv')
  escape_path = file.path(home.dir, 'resources/GeneEscape/mouse/escape_NatGen2016.txt')

# Lib
  library(tidyr)
  library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
  library(ggplot2)
  library(scales) # scale on X/Y axis
  library(ggpubr) # ggarange
  library(VennDiagram)
```

    ## Loading required package: grid

    ## Loading required package: futile.logger

    ## 
    ## Attaching package: 'VennDiagram'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     rotate

``` r
  library(zoo)
```

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

``` r
  source(file.path(prj_path, '/R_libs/Rlib_basic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_allelic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_allelicImbalance.R'))
  
# ggplot theme
  theme_man <- 
    theme_classic() +
    theme(
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size = 6),
      axis.title = element_text(colour='black', size=9),
      legend.title=element_blank(),
      legend.key=element_blank(),
      axis.line.x = element_line(colour = "black", size=0.25),
      axis.line.y = element_line(colour = "black", size=0.25),
      axis.ticks.x = element_line(colour = "black", size = 0.25),
      axis.ticks.y = element_line(colour = "black", size = 0.25),
      strip.background=element_blank(),
      strip.text=element_text(size=8))
  
  # Nice colors
  col_sel <- c("#39B54A", "#27AAE1", "#21409A", "#EF4136",  "violet", "grey")
    col.coding <- col_sel[1]
    col.nc <- col_sel[3]
```

``` r
  counts <- readRDS(counts_path)
  rpkm <- readRDS(rpkm_path)
  cast <- as.matrix(readRDS(cast_path))
  c57 <- as.matrix(readRDS(c57_path))

  biomart <- readRDS(biomart_path)
  imprint <- read.csv(imprint_path)
  meta <- readRDS(meta_path)
```

``` r
# Filter genes with few allelic counts
  nCounts = 3  # same as QC
  nCells = 20   # more strict then QC (2 cells)

  genesPass <- filter_genes_4allelicCounts(cast, c57, nCounts, nCells)
    cast_filt <- cast[genesPass, ]
    c57_filt <- c57[genesPass, ]
    biomart_filt <- biomart %>% filter(GeneStableID %in% genesPass)

# Autosomal genes only
  genesPass <- splitGenes_autosomeXimprint(rownames(cast_filt), biomart_filt, imprint)$autosomes
    cast_filt <- cast_filt[genesPass, ]
    c57_filt <- c57_filt[genesPass, ]
    biomart_filt <- biomart_filt %>% filter(GeneStableID %in% genesPass)

# Exlcude some GeneTypes
  biomart_filt <- 
    biomart_filt %>%
    filter(
      GeneType %notin% c('miRNA', 'snoRNA', 'TR_J_gene', 'rRNA', 'polymorphic_pseudogene', 'TR_C_gene', 'TEC', 'processed_pseudogene'),
      GeneStableID %in% rownames(cast_filt)) %>%
    select(GeneStableID, GeneName, GeneType, Strand, Chromosome, TSS, coding, lncRNA)
  
  rpkm_filt <- rpkm[biomart_filt$GeneStableID, ]
  cast_filt <- cast[biomart_filt$GeneStableID, ]
  c57_filt <- c57[biomart_filt$GeneStableID, ]
  
### Add allelicImbalance to biomart
  biomart_filt <-
    biomart_filt %>%
    mutate(
      allelicImbalance = dscore_gene(cast_filt, c57_filt, biomart_filt$GeneStableID),
      meanExpression = apply(rpkm[rownames(cast_filt), ], 1, mean)
  )
```

``` r
### Binomial test allelic imbalance

# Level of significance
  pBinomSig = 0.01
  
# Binomial test
  binomStats <- binomTest_AllelicImbalance(cast_filt, c57_filt)

  biomart_filt <- 
    biomart_filt %>%
    mutate(
      binomTest = binomStats[['pVal']],
      binomTestBH = binomStats[['pVal_adj']],
    ) %>%
    mutate(
      allelicImbalanceSig = ifelse(binomTestBH < pBinomSig, TRUE, FALSE)
    )

# Sig vs noSig by GeneType 
  biomart$allelicImbalanceSigBiotype <- FALSE
  
  lncRNA_pass <-
    biomart_filt %>%
    filter(!coding) %>%
    filter(binomTestBH < pBinomSig)
  coding_pass <-
    biomart_filt %>%
    filter(coding) %>%
    filter(binomTestBH < pBinomSig)
  
  biomart_filt <-
    biomart_filt %>%
    mutate(
      allelicImbalanceBiotype = ifelse(GeneStableID %in% coding_pass$GeneStableID, 'sig_coding', 'nonSig')
    )
  biomart_filt[match(lncRNA_pass$GeneStableID, biomart_filt$GeneStableID), ]$allelicImbalanceBiotype <- 'sig_lncRNA'
  
  # allelic imbalance vs p-values
  ggImbalance2pvalue <- 
    ggplot(biomart_filt, aes(x=allelicImbalance, y=-log10(binomTestBH), color=allelicImbalanceBiotype, size=allelicImbalanceBiotype, alpha=allelicImbalanceBiotype)) +
    geom_point(shape=16) + 
    scale_color_manual(values=c('grey', col.coding, col.nc)) +
    scale_size_manual(values=c(0.1, 0.15, 0.3)) +
    scale_alpha_manual(values=c(0.7, 0.7, 0.95)) +
    ylab('-log10(binomTest (BH adjusted)') +
    geom_hline(yintercept=-log10(pBinomSig), linetype="dashed", color = 'red', size=0.25) +
    geom_vline(xintercept=0, linetype="dashed", color = 'grey', size=0.25) +
    theme_man +
    theme(
      legend.position = 'none')
  
  # Density allelic imbalance
  ggImbalanceDensity <- 
    ggplot(biomart_filt, aes(x=allelicImbalance, color=coding)) +
    geom_density() +
    scale_color_manual(values=c(col.nc, col.coding)) +
    geom_vline(xintercept=0, linetype="dashed", color = 'grey', size=0.5) +
    xlim(-0.5, 0.5) +
    theme_man +
    theme(
      legend.position = 'none')
  
  
### Allelic imbalance vs mean expression  

  df.lncRNA <-
    biomart_filt %>%
    filter(lncRNA) %>%
    arrange(meanExpression) %>%
    mutate(
      allelicImbalance_abs = abs(.$allelicImbalance)
    ) %>%
    mutate(
      rollingMedian = zoo::rollmedian(allelicImbalance_abs, k=25, na.pad=TRUE)
    )
  
  df.coding <-
    biomart_filt %>%
    filter(coding) %>%
    arrange(meanExpression) %>%
    mutate(
      allelicImbalance_abs = abs(.$allelicImbalance)
    ) %>%
    mutate(
      rollingMedian = zoo::rollmedian(allelicImbalance_abs, k=25, na.pad=TRUE)
    )

  ggRolling <- 
    ggplot(rbind(df.lncRNA, df.coding), aes(x=meanExpression, y=rollingMedian, color=GeneType, size=GeneType)) + 
    geom_point(alpha=0.25) +
    scale_size_manual(values = c(0.1, 0.1)) +
    scale_color_manual(values=c(col.nc, col.coding)) +
    scale_x_log10(limits=c(0.01, 10)) + 
    ylim(0, 0.2) +
    annotation_logticks(sides = "b", size=0.25) +
    geom_smooth(data = df.lncRNA, aes(x=meanExpression, y=rollingMedian), size=0.25, colour=col.nc, fill = col.nc) +
    geom_smooth(data = df.coding, aes(x=meanExpression, y=rollingMedian), size=0.25, colour=col.coding, fill = col.coding) +
    ylab('| Allelic Imbalance |') +
    xlab('Mean expression (RPKM)') +
    ggtitle('Sliding window, k=25 genes') +
    theme_man +
    theme(
      legend.position = 'none')

  
  ggarrange(ggImbalance2pvalue, ggImbalanceDensity, ggRolling, ncol = 3)
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

    ## Warning: Removed 51 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 1360 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1411 rows containing missing values (geom_point).

![](alllelicImbalance_files/figure-markdown_github/Allelic%20imbalance%20for%20lncRNAs%20and%20mRNAs-1.png)

``` r
# Input parameters
  dist = 500000
  nRand = 1000
  minGenesInProx = 3

# Biomart input data (coding and lncRNAs)
  biomart_sel <- biomart_filt %>% select(GeneStableID, GeneName, GeneType, Chromosome, TSS, Strand, allelicImbalance, lncRNA, coding)

  biomart_lncRNAs <- biomart_sel %>% filter(!coding)
  biomart_coding <-  biomart_sel %>% filter(coding)

# Genes in proximity (of lncRNAs)
  lncRNAstats.l <- lapply(biomart_lncRNAs$GeneStableID, get_genesInProx_wStats, biomart_sel, dist)
  names(lncRNAstats.l) <- biomart_lncRNAs$GeneStableID
  
  # Exclude lncRNAs with no genes in prox
  genes.tmp <- vector()
  for(i in 1:length(lncRNAstats.l)){
    genes.tmp[i] <- nrow(lncRNAstats.l[[i]])
  }
  lncRNAstats.l <- lncRNAstats.l[which(genes.tmp >1)]
  

# Random sampling (moving lncRNA genes) to coding gene loci
  genesInProx.coding.l <- parallel::mclapply(biomart_coding$GeneStableID, get_genesInProx, biomart_sel, dist, mc.cores = 20)
  names(genesInProx.coding.l) <- biomart_coding$GeneStableID

    # Exclude genes with no/few genes in prox
    nGenes <- vector()
    
    for(i in 1:length(genesInProx.coding.l)){
      nGenes[i] <- nrow(genesInProx.coding.l[[i]])
    }
    
    genesInProx.coding.l <- genesInProx.coding.l[which(nGenes > minGenesInProx)]
    sampled_Genes <- sample(names(genesInProx.coding.l), nRand)
      
# Run test
  testRandom <- parallel::mclapply(biomart_lncRNAs$GeneStableID, randomSampling_allelicImbalance, biomart_sel, dist, sampled_Genes, mc.cores = 20)
  names(testRandom) <- biomart_lncRNAs$GeneStableID
  testRandom.df <- do.call(rbind.data.frame, testRandom)
  
# Remove NAs
  testRandom.df <- testRandom.df[rowSums(is.na(testRandom.df)) == 0, ]

# Calculate p-values and convert to df
  for(i in names(lncRNAstats.l)){
    tmp <- lncRNAstats.l[[i]]
    pTest = vector()
    
    for(j in 1:nrow(tmp)){
      pTest[j] <- length(which(testRandom[[i]]$allelicImbalanceTest >= tmp[j, ]$allelicImbalanceTest)) / nrow(testRandom[[i]])
    }
    lncRNAstats.l[[i]]$pTest <- pTest
  }
  
  lncRNAstats.df <- do.call(rbind.data.frame, lncRNAstats.l)
  
  lncRNAstats.df <-
    lncRNAstats.df %>%
    arrange(pTest) %>%
    mutate(
      pTest.sig = ifelse(pTest < 0.05, TRUE, FALSE),
      rank = seq(nrow(lncRNAstats.df))
    )

# Dump stats
  saveRDS(lncRNAstats.df, file=file.path(prj_path, 'lncRNAstats_allelicImbalance.rds'))
```

``` r
# Candidates
  selected_candidates <- lncRNAstats.df %>% filter(GeneName %in% c('Txnrd1', 'Gsta4', 'Tmc7', 'Fam78b'))

# Filter significant genes and extract cancidates
  tmp.df <-
    lncRNAstats.df %>%
    filter(pTest.sig) %>%
    arrange(allelicDirection*allelicImbalanceTest) %>%
    mutate(nGene = seq(nrow(.)))
  
  tmp.df.sel <- tmp.df %>% filter(GeneName %in% selected_candidates$GeneName)
  
# Plot
  ggRank_v1 <- 
    ggplot(tmp.df, aes(x=allelicDirection*allelicImbalanceTest, y=nGene)) +
    geom_point(shape=16, color='red', size=0.5, alpha=0.95) +
    xlim(-1, 1) +
    geom_text(data=tmp.df.sel, aes(x=allelicDirection*allelicImbalanceTest, y=nGene, label=GeneName), color='black', size = 4, vjust=-0.5) + 
    geom_point(data=tmp.df.sel, aes(x=allelicDirection*allelicImbalanceTest, y=nGene), color='black', shape=1, size=0.75, alpha=0.95) + 
    geom_vline(xintercept=0, linetype="dashed", color = 'black', size=0.25) +
    ylab('Gene-pairs, sorted by Ranking') +
    xlab('Gene-pairs, sorted by allelic imbalance') +
    theme_man

  ggRank_v2 <- 
    ggplot(lncRNAstats.df, aes(x=-log10(pTest), y=abs(distToTSS), color=pTest.sig, size=pTest.sig)) +
    geom_point(shape=16, alpha=0.95) +
    scale_y_continuous(breaks=c(0, dist/2, dist)) +
    scale_color_manual(values = c('grey', 'red')) +
    scale_size_manual(values = c(0.15, 0.5)) +
    geom_text(data=tmp.df.sel, aes(x=-log10(pTest), y=abs(distToTSS), label=GeneName), color='black', size = 3, vjust=-0) + 
    geom_point(data=tmp.df.sel, aes(x=-log10(pTest), y=abs(distToTSS)), color='black', shape=1, alpha=0.95) + 
    geom_vline(xintercept=-log10(0.05), linetype="dashed", color = 'red', size=0.25) +
    ylab('| Distance between TSSs (mRNA-lncRNA) |') +
    xlab('-log10(p-value)') +
    theme_man +
    theme(
      legend.position = 'none')

  
  ggarrange(ggRank_v1, ggRank_v2, ncol = 2)
```

![](alllelicImbalance_files/figure-markdown_github/Plot%20significant%20candidates%20and%20highlight%20validated%20candidates-1.png)

``` r
# Distribution of gene density to TSS
  lncRNAstats.df$absDistToTSS <- abs(lncRNAstats.df$distToTSS)
  testRandom.df$absDistToTSS <- abs(testRandom.df$distToTSS)
  
  tmp <- seq(0, dist, dist/20)
  gg1 <- 
    ggplot(lncRNAstats.df %>% filter(pTest.sig), aes(x=absDistToTSS)) + 
    geom_histogram(color='red', fill=NA, size=0.25, breaks=tmp) + 
    scale_x_continuous(breaks=c(0, dist/2, dist)) +
    xlab('Distance to TSS') +
    ylab('Occurence') +
    ggtitle(paste0('nSig lncRNA-gene interactions : ', nrow(lncRNAstats.df %>% filter(pTest.sig)))) +
    theme_man +
    theme(
      legend.position = 'none',
      plot.title = element_text(size = 8))
  
  gg2 <- 
    ggplot(lncRNAstats.df, aes(x=absDistToTSS)) + 
    geom_histogram(color=col.nc, fill=NA, size=0.25, breaks=tmp) + 
    scale_x_continuous(breaks=c(0, dist/2, dist)) +
    xlab('Distance to TSS') +
    ylab('Occurence') +
    ggtitle(paste0('lncRNA-gene interactions; ', nrow(lncRNAstats.df))) +
    theme_man +
    theme(
      legend.position = 'none',
      plot.title = element_text(size = 8))
  
  gg3 <- 
    ggplot(testRandom.df, aes(x=absDistToTSS)) + 
    geom_histogram(color='grey', fill=NA, size=0.25, breaks=tmp) + 
    scale_x_continuous(breaks=c(0, dist/2, dist)) +
    xlab('Distance to TSS') +
    ylab('Occurence') +
    ggtitle(paste0('Randomly generated lncRNA-mRNA gene pairs: ', nrow(testRandom.df))) +
    theme_man +
    theme(
      legend.position = 'none',
      plot.title = element_text(size = 8))
  
  
    ggarrange(gg1, gg2, gg3, ncol = 3, nrow = 1)
```

![](alllelicImbalance_files/figure-markdown_github/Plot%20distribution%20of%20lncRNA-mRNA%20gene%20pairs-1.png)

``` r
# Fnc to make nice plots
  make_nicePlot <- function(geneSel='Txnrd1', biomart, dist){
    # Selected gene
    biomartSel <- biomart %>% filter(GeneStableID %in% geneSel | GeneName %in% geneSel)
  
    # Genes in proximity
    df <-
      biomart %>%
      filter(
        Chromosome %in% biomartSel$Chromosome
      ) %>%
      mutate(
        distToTSS = TSS - biomartSel$TSS,
        distToTSS_abs =  abs(TSS - biomartSel$TSS)
      ) %>%
      filter(distToTSS_abs < dist)
  
    gg <- 
      ggplot(df, aes(x=distToTSS, y=allelicImbalance, label=GeneName)) +
      geom_point(alpha=0.95, shape=16, size=0.5) +
      geom_hline(yintercept=0, linetype="dashed", color = 'black', size=0.25) +
      geom_text(aes(x=distToTSS, y=allelicImbalance, label=GeneName), color='black', size = 3, vjust=-0.5) + 
      geom_hline(yintercept=c(-0.25, 0.25), linetype="dashed", color = 'grey40', size=0.1) +
      ylim(-0.5, 0.5) +
      xlim(-dist, dist) +
      xlab('Distance to lncRNA-TSS') +
      ylab('Allelic imbalance') +
      ggtitle(geneSel) +
      theme_man +
      theme(
        plot.title = element_text(size = 8)
            )
    return(gg)
  }

  gg1 <- make_nicePlot(gene='1700028I16Rik', biomart_filt, dist=500000)
  gg2 <- make_nicePlot(gene='B230311B06Rik', biomart_filt, dist=500000)
  gg3 <- make_nicePlot(gene='C920006O11Rik', biomart_filt, dist=500000)
  gg4 <- make_nicePlot(gene='Gm16701', biomart_filt, dist=500000)

  ggarrange(gg1, gg2, gg3, gg4, nrow=1, ncol = 4)
```

![](alllelicImbalance_files/figure-markdown_github/Plot%20allelic%20imbalance%20of%20validated%20candidates-1.png)
