``` r
# Load data, packages and set output
  options(stringsAsFactors = FALSE)
  `%notin%` <- Negate(`%in%`)

# Basic folders
  home.dir <- '/home/perj/'

# Project path
  prj_name <- 'prj_lncRNAs'
  prj_path = file.path(home.dir, 'projects', prj_name)

# Data
  rpkm_path = file.path(prj_path, 'data/ActD_fibroblasts_rpkm.rds')
  counts_path = file.path(prj_path, 'data/ActD_fibroblasts_counts.rds')
  meta_path = file.path(prj_path, 'data/ActD_fibroblasts_meta.rds')
  biomart_path = file.path(prj_path, '/data/GRCm38.p6.csv')
  
  halfLives_ES_path <- file.path(prj_path, 'data/halfLife_mES_Herzog2017.csv')

# Time points  
  timePoints = c(0, 1, 2, 4, 7, 10)
  delayUptate = 5/60
  timePoints = timePoints - delayUptate
  timePoints[which(timePoints<0)] <- 0

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
  library(scales)
  library(ggpubr)
  library(hexbin)
  library(ggpointdensity)
  library(viridis)
```

    ## Loading required package: viridisLite

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:viridisLite':
    ## 
    ##     viridis.map

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
  source(file.path(prj_path, '/R_libs/Rlib_basic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_RNAdecay.R'))

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
    col.nc  <- col_sel[3]
    col.coding <- col_sel[1]
```

``` r
  rpkm <- readRDS(rpkm_path)
  counts <- readRDS(counts_path)
  meta <- readRDS(meta_path)
  biomart <- read.csv(biomart_path, header=TRUE)

  dataES <- read.csv(halfLives_ES_path, header=TRUE, sep=';')
```

``` r
# Thresholds
  minCounts = 8
  minSamples = 4
  
  minCountsSample = 50000
  
# Apply filters
  
  # meta
  meta.qc <- 
    meta %>%
    mutate(
      readCounts = colSums(counts[ ,.$XC])
    ) %>%
    filter(readCounts > minCountsSample)

  # QC genes (for expression in t0)
  counts_t0 <-  counts[ ,(meta.qc %>% filter(treatment %in% 't0'))$XC]
  genes_pass <- filter_genes(counts_t0, minCounts, minSamples)
  
# QC biomart
  biomart.qc <-
    biomart %>%
    filter(GeneStableID %in% genes_pass) %>%
    distinct(GeneStableID, .keep_all = TRUE) %>%
    filter(GeneType %notin% c('TR_C_gene', 'TEC', 'Mt_tRNA', 'Mt_rRNA', 'misc_RNA')) %>%
    select(GeneStableID, GeneName, GeneType, Chromosome)
  
  rpkm.qc <- rpkm[biomart.qc$GeneStableID, meta.qc$XC]
  counts.qc <- counts[biomart.qc$GeneStableID, meta.qc$XC]
```

``` r
# Genes to use for normalization:
  # 50 read counts in all four t0 samples
  # 10 read counts in all other samples
  # Genes with 1h < half-life < 8h in mES data set (Herzog2017)
  # protein coding only

# Thresholds
  minSamplesForNorm_t0 = 4  
  
  minCountsForNorm_t0 = 50
  minCountsForNorm = 10
  
  halfLife_limits <- c(1, 8)

# Identify genes to use for normalization
  # Limit to genes with detection in all 4 samples with at least 10 (minCountsForNorm_t0) counts
    genes_pass_1 <- 
      counts.qc[ ,(meta.qc %>% filter(treatment %in% 't0'))$XC] %>%
      filter_genes(minCountsForNorm_t0, minSamplesForNorm_t0)
  
  # Limit to genes with detection in all other samples (20) with at least 10 (minCountsForNorm) counts
    meta.tmp <- meta.qc %>% filter(treatment %notin% 't0')
  
    genes_pass_2 <- 
      counts.qc[ ,meta.tmp$XC] %>%
      filter_genes(minCountsForNorm, nrow(meta.tmp))
  
  # Use genes passing both thresholds
  genes_pass <- intersect(genes_pass_1, genes_pass_2)

# Identify genes with half-life estimates in ES cell data set
  genesToUseForNorm.df <- 
    dataES %>%
    filter(GeneName %in% (biomart.qc %>% filter(GeneType %in% 'protein_coding' & GeneStableID %in% genes_pass))$GeneName) %>%
    filter(HalfLife_h > halfLife_limits[1] & HalfLife_h < halfLife_limits[2]) %>%
    select(GeneName, HalfLife_h) %>%
    mutate(
      k = getDecayRate(.$HalfLife_h)
    )
    

# Calculate expected expression levels (given decay / half-lifes)
  genesToUseForNorm.df <- 
    genesToUseForNorm.df %>%
    mutate(
      t1 = apply(genesToUseForNorm.df %>% select(k), 1, getExpectedExpression, t=timePoints[2], expression.t0=1),
      t2 = apply(genesToUseForNorm.df %>% select(k), 1, getExpectedExpression, t=timePoints[3], expression.t0=1),
      t3 = apply(genesToUseForNorm.df %>% select(k), 1, getExpectedExpression, t=timePoints[4], expression.t0=1),
      t4 = apply(genesToUseForNorm.df %>% select(k), 1, getExpectedExpression, t=timePoints[5], expression.t0=1),
      t5 = apply(genesToUseForNorm.df %>% select(k), 1, getExpectedExpression, t=timePoints[6], expression.t0=1)
    ) %>%
    mutate(
      t0 = 1,
      GeneStableID = biomart.qc[match(.$GeneName, biomart.qc$GeneName), ]$GeneStableID
    ) %>%
    select(GeneStableID, GeneName, HalfLife_h, k, t0, t1, t2, t3, t4, t5)
```

``` r
# Normalize the expression to t0 for each replicate
  normToT0 <- function(sampleID, meta, expr){
    meta.tmp <- meta %>% filter(ID %in% sampleID)
    meta.t0 <- meta %>% filter(ID %in% sampleID) %>% filter(treatment %in% 't0')
    expr_norm <- expr[ ,meta.tmp$XC] / expr[ ,meta.t0$XC]
    return(expr_norm)
  }

  rpkm.t0 <-
    cbind(
      normToT0('7_1', meta.qc, rpkm.qc),
      normToT0('7_2', meta.qc, rpkm.qc),
      normToT0('8_1', meta.qc, rpkm.qc),
      normToT0('8_2', meta.qc, rpkm.qc)
      )

### Normalize the expression over the individual time points
  samples <- c('7_1', '7_2', '8_1', '8_2')
  
  rpkmNorm.l = list()
  rpkmNorm.l[[1]] <- normData_bySubSetGene(rpkm.t0, meta.qc, biomart.qc, samples[1], genesToUseForNorm.df)
  rpkmNorm.l[[2]] <- normData_bySubSetGene(rpkm.t0, meta.qc, biomart.qc, samples[2], genesToUseForNorm.df)
  rpkmNorm.l[[3]] <- normData_bySubSetGene(rpkm.t0, meta.qc, biomart.qc, samples[3], genesToUseForNorm.df)
  rpkmNorm.l[[4]] <- normData_bySubSetGene(rpkm.t0, meta.qc, biomart.qc, samples[4], genesToUseForNorm.df)
  names(rpkmNorm.l) <- c('7.1', '7.2', '8.1', '8.2')
  rpkmNorm <- do.call(cbind, rpkmNorm.l)
  
# Observed half-life (decay) parameters
  # Note: Genes with half-life < 2h are not used to normalize the data for later times points (7 and 10h)

  library(aomisc)
```

    ## Loading required package: drc

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 
    ## 'drc' has been loaded.

    ## Please cite R and 'drc' if used for a publication,

    ## for references type 'citation()' and 'citation('drc')'.

    ## 
    ## Attaching package: 'drc'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     gaussian, getInitial

    ## Loading required package: plyr

    ## ------------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## ------------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     mutate

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## Loading required package: car

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

    ## Loading required package: multcompView

``` r
  decayParameters.df <- 
    do.call(
      rbind,
      parallel::mclapply(rownames(rpkmNorm.l[[1]]), getDecayParameters, rpkmNorm.l, meta.qc, mc.cores = 20)) %>%
    mutate(
      GeneName = biomart.qc[match(.$GeneStableID, biomart.qc$GeneStableID), ]$GeneName
    )
  detach("package:aomisc", unload=TRUE)
  detach("package:drc", unload=TRUE)
  detach("package:MASS", unload=TRUE)
```

    ## Warning: 'MASS' namespace cannot be unloaded:
    ##   namespace 'MASS' is imported by 'TH.data' so cannot be unloaded

``` r
# Filter genes (fibs)
  fib.filt <-
    decayParameters.df %>%
    filter(HalfLife_h > 0 & HalfLife_h < 10 & k < 2) %>%
    mutate(
      cellLine = 'Fibs',
      GeneName = biomart.qc[match(.$GeneStableID, biomart.qc$GeneStableID), ]$GeneName
      )

# Filter genes (mESC, Herzog et al 2017)
  ES.filt <- 
    dataES %>%
    distinct(GeneName, .keep_all = TRUE) %>%
    filter(HalfLife_h < 10 & HalfLife_h > 0 & k < 2) %>%
    mutate(
      cellLine = 'mES'
    )
  
# Compare estimate from fibroblast and mES (for genes passing QC)
  # Filtered genes present for fibs & mESC
  genes_pass <- intersect(ES.filt$GeneName, fib.filt$GeneName)
  
  df <- 
    rbind(
    fib.filt %>% filter(GeneName %in% genes_pass) %>% select(GeneName, HalfLife_h, cellLine),
    ES.filt %>% filter(GeneName %in% genes_pass) %>% select(GeneName, HalfLife_h, cellLine)
    ) %>%
    spread(cellLine, HalfLife_h)

  tmp.cor <- round(cor(df$Fibs, df$mES), 3)

# Plot Extended Data Fig 4d
  ggplot(df, aes(x=Fibs, y=mES)) +
    geom_pointdensity(shape=16, size=0.5, alpha=0.95) +
    scale_color_viridis() +
    geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.25) +
    ggtitle(tmp.cor) +
    theme_man +
    theme(legend.position = c(0.8, 0.95))
```

![](RNAdecay_ActD_fibroblast_files/figure-markdown_github/Filter%20and%20compare%20to%20previous%20estimates%20from%20ES%20cells%20(Herzog%20et%20al%202017)-1.png)

``` r
# Add biotype info to estimates
  fib.filt.sel <- 
    fib.filt %>%
    mutate(
      biotype = biomart.qc[match(.$GeneStableID, biomart.qc$GeneStableID), ]$GeneType
    ) %>%
    filter(.$biotype %in% c('protein_coding', 'lncRNA'))

# Half.lives
  p = round(
    t.test(
      (fib.filt.sel %>% filter(biotype %in% 'lncRNA'))$HalfLife_h, 
      (fib.filt.sel %>% filter(biotype %in% 'protein_coding'))$HalfLife_h
    )$p.value, 3)
  
  ggHalf <- 
    ggplot(fib.filt.sel, aes(y=HalfLife_h, x=biotype, color=biotype)) +
    geom_violin(size=0.25) +
    scale_color_manual(values=c(col.nc, col.coding)) +
    scale_y_continuous(breaks = c(0, 5, 10, 15), limits=c(0, 15)) +
    geom_boxplot(lwd=0.25, outlier.colour = NA, width=0.35, fatten=1) + 
    ggtitle(paste0('p: ', p)) +
    theme_man +
    theme(legend.title = element_text(size = 6), 
          legend.text  = element_text(size = 6),
          legend.position = c(0.8, 0.95))
  
# Decay rate
  p = round(
    t.test(
      (fib.filt.sel %>% filter(biotype %in% 'lncRNA'))$k, 
      (fib.filt.sel %>% filter(biotype %in% 'protein_coding'))$k
    )$p.value, 4)
  
  ggK <- 
    ggplot(fib.filt.sel, aes(y=k, x=biotype, color=biotype)) +
    geom_violin(size=0.25) +
    scale_color_manual(values=c(col.nc, col.coding)) +
    scale_y_log10() +
    annotation_logticks(sides='l', size=0.25) +
    geom_boxplot(lwd=0.25, outlier.colour = NA, width=0.35, fatten=1) + 
    ylab('Decay rate') +
    ggtitle(paste0('p: ', p)) +
    theme_man +
    theme(legend.title = element_text(size = 6), 
          legend.text  = element_text(size = 6),
          legend.position = c(0.8, 0.95))

    ggarrange(ggHalf, ggK, ncol = 2, nrow=1)
```

![](RNAdecay_ActD_fibroblast_files/figure-markdown_github/Compare%20mRNAs%20to%20lncRNAs-1.png)

``` r
controlPlots <- lapply(c('Fos', 'Myc', 'Actb', 'Gapdh'), plotModelToRealData, rpkmNorm, meta.qc, biomart.qc, fib.filt)

ggarrange(controlPlots[[1]], controlPlots[[2]], controlPlots[[3]], controlPlots[[4]], ncol = 4, nrow=1)
```

![](RNAdecay_ActD_fibroblast_files/figure-markdown_github/Conttrol%20plots,%20fit%20to%20observed%20expression-1.png)
