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
  rpkm_path = file.path(prj_path, 'data/ss2_n533_RPKM.rds')
  counts_path = file.path(prj_path, 'data/ss2_n533_readCounts.rds')
  cast_path <- file.path(prj_path, 'data/ss2_n533_castCounts.rds')
  c57_path <- file.path(prj_path, 'data/ss2_n533_c57Counts.rds')
  meta_path = file.path(prj_path, 'data/ss2_n533_meta.rds')
  biomart_path = file.path(prj_path, '/data_dump/ss2_fibs_biomart_wGeneAnnotations.rds')
    
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
  library(zoo)
```

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

``` r
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
  cast <- readRDS(cast_path)
  c57 <- readRDS(c57_path)
  meta <- readRDS(meta_path)
  biomart <- readRDS(biomart_path)
  
  imprint <- read.csv(imprint_path)
  escape <- read.csv(escape_path, header=TRUE)
  
  
# QC genes (note that data input is already QC)
  min.counts = 5
  min.cells = 2
  
  counts_qc <- counts[filter_genes(counts, min.counts, min.cells), ]
  
  biomart_qc <- 
    biomart %>%
    filter(GeneStableID %in% rownames(counts_qc)) %>%
    distinct(GeneName, .keep_all = TRUE)
```

``` r
# QC threshold
minReads = 5e5

# Extended Data Fig. 1a
  ggplot(meta, aes(y=Reads, x='')) +
  geom_boxplot(outlier.shape = NA, lwd=0.5, fatten=1) +
  geom_jitter(position = position_jitter(0.25), shape=20, size=0.5, alpha=0.25) +
  scale_color_manual(aes(colour = "lightgrey")) +
  scale_y_log10(limits=c(1e5, 1e8), label=scientific_10) +
  annotation_logticks(sides='l', size=0.25) +
  ylab('Read counts mapped to exons (log10)') +
  geom_hline(yintercept=minReads, linetype="dashed", color = "red", size=0.5) +
  theme_man +
  theme(
    legend.position = "none",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
```

![](2_makeQCplots_SS2_mouseFibroblasts_files/figure-markdown_github/QC%20read%20counts%20Fig%20S1A-1.png)

``` r
# Filter genes (for allelic reads)
  n.counts = 3
  n.cells = 2
  
  auto_pass = c(-0.1, 0.1) 
  X_pass = c(-0.4, 0.4)

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
ggAuto <-   
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

ggX <- 
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

ggarrange(ggAuto, ggX)
```

![](2_makeQCplots_SS2_mouseFibroblasts_files/figure-markdown_github/Allelic%20Distribution%20of%20read%20counts,%20Extended%20Data%20Fig.%201BC-1.png)

``` r
# Threshold for independent gene promoter
  dist_independent = 4000

# Rolling median
  k=51 

  df.tmp <-
    biomart_qc %>%
    mutate(
      mean = apply(rpkm[biomart_qc$GeneStableID, ], 1, mean)
      ) %>%
    arrange(minDist2TSS) %>%
    mutate(
      roll.med = rollmedian(mean, k=k, na.pad=TRUE)
      )

  ggplot(df.tmp, aes(x=minDist2TSS, y=roll.med)) +
    geom_point(shape=16, size=0.01,  colour='black') +
    geom_smooth(size=0.25) +
    xlim(0, 25000) +
    ylab('Median expression (sliding window, RPKM)') +
    geom_vline(xintercept=dist_independent, linetype="dashed", size=0.25, color='red') +
    theme_man
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 10195 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 10195 rows containing missing values (geom_point).

![](2_makeQCplots_SS2_mouseFibroblasts_files/figure-markdown_github/Expression%20(RPKM)%20vs%20Distance%20between%20two%20TSSs,%20Extended%20Data%20Fig.%201d-1.png)

``` r
# Fnc to count genes per cell
  genesPerCell <- function(biomart_sel, exp, meta, threshold) {
    nGenes <- colSums(exp[(biomart_sel$GeneStableID), ] >= threshold)
    return(nGenes)
  }

# Genes per cell
  countThr = 3  
  
  meta <- 
    meta %>%
    mutate(
      coding = genesPerCell(biomart_qc %>% filter(GeneStableID %in% rownames(counts_qc) & coding), counts, meta=., countThr),
      lncRNA = genesPerCell(biomart_qc %>% filter(GeneStableID %in% rownames(counts_qc) & lncRNA), counts, meta=., countThr),
      intergenic = genesPerCell(biomart_qc %>% filter(GeneStableID %in% rownames(counts_qc) & lncRNA & intergenic), counts, meta=., countThr),
      antisense = genesPerCell(biomart_qc %>% filter(GeneStableID %in% rownames(counts_qc) & lncRNA) %>% filter(convergent | divergent), counts, meta=., countThr)
    )

# Plot
gg_coding <-
  ggplot(meta, aes(y=coding, x='')) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1 ), lwd=0.5, fatten=1) + 
  geom_jitter(width = 0.25, alpha=0.1, size=0.1, color='black') +
  scale_color_manual(aes(colour = "lightgrey")) +
  scale_y_continuous(limits=c(0, 12000), breaks=seq(0,12000,2000)) +
  ylab('coding per cell (>= 3 read counts)') +
  theme_man +
  theme(
    axis.ticks.x=element_blank(),
    legend.title = element_blank(),
    axis.title.x=element_blank(),
    legend.position = 'none')


gg_noncoding <-
  meta %>%
  select(lncRNA, intergenic, antisense) %>%
  gather('biotype', 'genes') %>%
  
  ggplot(aes(y=genes, x=factor(biotype, level = c('lncRNA', 'antisense', 'intergenic')), color=biotype)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1), lwd=0.5, fatten=1) + 
    geom_jitter(width = 0.25, alpha=0.1, size=0.1, aes(colour = biotype)) +
    scale_color_manual(values=col_sel[1:3]) +
    scale_y_continuous(limits=c(0, 1050), breaks=seq(0,1050,250)) +
    ylab('LncRNA per cell (>= 3 read counts)') +
    theme_man +
    theme(
      axis.ticks.x=element_blank(),
      legend.title = element_blank(),
      axis.text.x=element_text(size=6),
      axis.title.x=element_blank(),
      legend.position = 'none') 


ggarrange(gg_coding, gg_noncoding, widths=c(1/3, 2/3))
```

![](2_makeQCplots_SS2_mouseFibroblasts_files/figure-markdown_github/Genes%20per%20cell,%20coding%20and%20lncRNAs-1.png)
