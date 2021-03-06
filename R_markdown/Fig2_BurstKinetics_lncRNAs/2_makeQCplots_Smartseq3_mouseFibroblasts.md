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
  rpkm_path = file.path(prj_path, 'data/ss3_n682_fibs_RPKM.rds')
  counts_path = file.path(prj_path, 'data/ss3_n682_fibs_readCounts.rds')
  umi_path = file.path(prj_path, 'data/ss3_n682_fibs_umiCounts.rds')
  cast_path <- file.path(prj_path, 'data/ss3_n682_fibs_castCounts.rds')
  c57_path <- file.path(prj_path, 'data/ss3_n682_fibs_c57Counts.rds')
  meta_path = file.path(prj_path, 'data/ss3_n682_fibs_meta.rds')
  biomart_path = file.path(prj_path, '/data_dump/ss3_fibs_biomart_wGeneAnnotations.rds')
    
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
  library(scales)
  library(ggpubr)
  library(hexbin)
  
  source(file.path(prj_path, '/R_libs/Rlib_basic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_genetypes.R'))
  source(file.path(prj_path, '/R_libs/Rlib_allelic.R'))

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
```

``` r
# Load data
  rpkm <- readRDS(rpkm_path)
  counts <- readRDS(counts_path)
  umi <- readRDS(umi_path)
  cast <- readRDS(cast_path)
  c57 <- readRDS(c57_path)
  meta <- readRDS(meta_path)
  biomart <- readRDS(biomart_path)
  
  imprint <- read.csv(imprint_path)
  escape <- read.csv(escape_path, header=TRUE)
  
  
# QC genes (note that data input is already QC)
  min.counts = 1
  min.cells = 1
  
  counts_qc <- counts[filter_genes(counts, min.counts, min.cells), ]
  
  biomart_qc <- 
    biomart %>%
    filter(GeneStableID %in% rownames(counts_qc)) %>%
    distinct(GeneName, .keep_all = TRUE)
```

``` r
# QC threshold
  minReadsTot = 5e5
  minUMIsTot = 1e5

# Add stats to meta
  meta <- 
    meta %>%
    mutate(
      nUMIs = colSums(umi[ , .$Sample]),
      nGenesUMI = colSums(umi[ , .$Sample] >= 1),
      nReadCounts = colSums(counts_qc[ ,.$Sample])
      )
  
  gg_S3A <- 
    ggplot(meta, aes(x=nReadCounts, y=nUMIs)) +
    geom_hex(bins=75) +
    annotation_logticks(sides='bl', size = 0.25, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    scale_x_log10(limits = c(100, 1e7)) + 
    scale_y_log10(limits = c(100, 1e6)) +
    scale_fill_continuous(type = "viridis") +
    geom_hline(yintercept = minUMIsTot, color="red", linetype="dashed", size=0.25) +
    geom_vline(xintercept = minReadsTot, color="red", linetype="dashed", size=0.25) +
    theme_man + 
    theme(legend.position = c(0.2, 0.3))

  ggarrange(gg_S3A)
```

![](2_makeQCplots_Smartseq3_mouseFibroblasts_files/figure-markdown_github/QC%20read%20counts%20to%20UMIs,%20Fig%20S3A-1.png)

``` r
#  QC thresholds
  n.counts = 1
  n.cells = 1
  auto_pass = c(-0.1, 0.1)
  X_pass = c(-0.4, 0.4)
  
# Filter genes
  passed_genes <- filter_genes_4allelicCounts(cast, c57, n.counts, n.cells)
    cast_qc <- cast[passed_genes, ]
    c57_qc <- c57[passed_genes, ]

  biomart_qc_allelic <-
    biomart_qc %>%
    filter(GeneStableID %in% passed_genes) %>%
    distinct(GeneStableID, .keep_all = TRUE)

# Retrieve autosomal/X/imprinted genes
  tmp <- splitGenes_autosomeXimprint(rownames(cast_qc), biomart_qc_allelic, imprint)

# Exclude escapee genes on X
  biomart_X <- 
    biomart_qc_allelic %>%
    filter(GeneStableID %in% tmp$X) %>%
    filter(GeneName %notin% escape$GeneName)

# Autosomes (imprinted genes already excluded)
  biomart_auto <- 
    biomart_qc_allelic %>%
    filter(GeneStableID %in% tmp$autosomes)
  
# Distribution of allelic read counts (dScore)
  dScore.df <- 
    data.frame(
      auto = allelic_balance_perCell(cast_qc, c57_qc, biomart_auto$GeneStableID),
      X = allelic_balance_perCell(cast_qc, c57_qc, biomart_X$GeneStableID),
      cell.random = sample(ncol(cast_qc))
    )
  
# Add to meta
  meta <- 
    meta %>%
    mutate(
      dScore.X = dScore.df[.$Sample, ]$X,
      dScore.auto = dScore.df[.$Sample, ]$auto
      ) %>%
    mutate(
      dScore.X.pass = ifelse(.$dScore.X < X_pass[1] | .$dScore.X > X_pass[2], TRUE, FALSE),
      dScore.auto.pass = ifelse(.$dScore.auto > auto_pass[1] & .$dScore.auto < auto_pass[2], TRUE, FALSE)
    )
      
# Plot
gg_S3B <-   
  ggplot(dScore.df, aes(x=cell.random, y=auto), main='test') + 
  geom_point(alpha=0.75, shape=16, colour = "black", size=2) + 
  ylim(-0.5,0.5) +
  theme_classic() +
  xlab('Sorted by cell') +
  ylab('Allelic distribution (autosomes)') +
  geom_hline(yintercept=c(-0.1, 0.1), linetype="dashed", color = "red", size=0.25) +
  geom_hline(yintercept=0, linetype="dashed", color = "lightgrey", size=0.25) +
  theme(legend.position = 'none') +
  theme_man

gg_S3C <- 
  ggplot(dScore.df, aes(x=cell.random, y=X)) + 
  geom_point(alpha=0.75, shape=16, colour = "black", size=2) + 
  ylim(-0.5, 0.5) +
  theme_classic() +
  xlab('Sorted by cell') +
  ylab('Allelic distribution (X)') +
  geom_hline(yintercept=c(-0.4, 0.4), linetype="dashed", color = "red", size=0.25) +
  geom_hline(yintercept=0, linetype="dashed", color = "lightgrey", size=0.25) +
  theme(legend.position = 'none') +  
  theme_man

ggarrange(gg_S3B, gg_S3C)
```

![](2_makeQCplots_Smartseq3_mouseFibroblasts_files/figure-markdown_github/QC%20allelic%20Distribution%20of%20SS3%20read%20counts,%20Extended%20Data%20Fig.%20S3B,%20S3C-1.png)
