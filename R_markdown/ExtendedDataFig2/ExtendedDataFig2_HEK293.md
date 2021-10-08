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
  rpkm_path = file.path(prj_path, 'data/ss3_n112_HEK293_rpkm.rds')
  counts_path = file.path(prj_path, 'data/ss3_n112_HEK293_counts.rds')
  meta_path = file.path(prj_path, 'data/ss3_n112_HEK293_meta.rds')
  biomart_path = file.path(prj_path, '/data/GRCh38.p12.csv')
  gtf_path = file.path(home.dir, 'resources/Biomart/Homo_sapiens.GRCh38.91.chr.biotypes.txt')
  
# Resources
  imprint_path = file.path(home.dir, 'resources/GeneImprint/human.csv')
  
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
  library(GenomicRanges)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

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
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomeInfoDb

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
  source(file.path(prj_path, '/R_libs/Rlib_genetypes.R'))

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
    col_coding = col_sel[1]
    col_nc = col_sel[3]
```

``` r
# Thresholds (Note: input data here is already post QC)
  minCells = 3
  minCounts = 1
  minGenes = 7500
  minReads = 5e5

# Meta
  meta <-
    meta %>%
    mutate(
      Reads = colSums (counts[ ,.$Sample]),
      nGenes = colSums(counts[ ,.$Sample] >= minCounts)
    ) %>%
    mutate(
      Read.pass = ifelse(.$Reads >= minReads, TRUE, FALSE),
      nGenes.pass = ifelse(.$nGenes >= minGenes, TRUE, FALSE)
    )
  
  meta.qc <- meta %>% filter(Reads.pass & nGenes.pass)  
  
# Filter genes
  genes_passed <- filter_genes(counts[ ,meta.qc$Sample], minCounts, minCells)
  
  counts.qc <- counts[genes_passed, meta.qc$Sample]
  rpkm.qc <- rpkm[genes_passed, ]
  
# Use annotations from gtf file
  biomart.qc <- 
    biomart %>%
    filter(GeneStableID %in% genes_passed) %>%
    distinct(GeneStableID, .keep_all = TRUE) %>%
    mutate(
      TranscriptType_gtf = gtf_geneTypes[match(.$GeneStableID, gtf_geneTypes$V1), ]$V2  
    ) %>%
    mutate(
      coding = ifelse(.$TranscriptType_gtf == 'protein_coding', TRUE, FALSE),
      lncRNA = ifelse(
        .$TranscriptType_gtf == 'lincRNA' | 
          .$TranscriptType_gtf == 'antisense_RNA' |
          .$TranscriptType_gtf == 'bidirectional_promoter_lncRNA',
          TRUE, FALSE)
    )
```

``` r
# Define transcript types
  dist_divergent = 500
  dist_convergent = 2000
  dist_independent = 4000
  dist_intergenic = 4000

# Identify distance to gene in proximity (limited to QC expressed genes)
  genesInProx.l <- parallel::mclapply(biomart.qc$GeneStableID, find_minDist2TSS, counts, biomart.qc, mc.cores = 20)
  genesInProx.df <- as.data.frame(do.call(rbind, genesInProx.l))

  biomart.qc <- cbind(biomart.qc, genesInProx.df)

# Annotate transcript type
  
  # Independent
  biomart.qc <- 
    biomart.qc %>%
    mutate(
      independent = ifelse(minDist2TSS > dist_independent, TRUE, FALSE)
    )
  
  # Divergent
  biomart.qc <- 
    biomart.qc %>%
    mutate(
      divergent =
        as.vector(
          do.call(rbind,
                  parallel::mclapply(biomart.qc$GeneStableID, find_Divergent, biomart.qc, dist_divergent, mc.cores = 20)))
      )
  
  # Convergent  
  biomart.qc <- 
    biomart.qc %>%
    mutate(
      convergent = 
        as.vector(
          do.call(rbind,
                  parallel::mclapply(biomart.qc$GeneStableID, find_Convergent, biomart.qc, dist_convergent, mc.cores = 20)))
    )
  
# Intergenic genes
  GR <-
    GRanges(
    seqnames = biomart.qc$Chromosome,
    ranges = IRanges(start = biomart.qc$GeneStart,
                     end = biomart.qc$GeneEnd),
    strand = biomart.qc$GeneStrand
  )
  names(GR) <- biomart.qc$GeneName

  ol1 <- findOverlaps(GR, select='first', maxgap=dist_intergenic)
  ol2 <- findOverlaps(GR, select='last', maxgap=dist_intergenic)

  tmp <- data.frame(GeneName = biomart.qc$GeneName, ol1, ol2)
  biomart.qc$GeneOverlapEnsemble <- biomart.qc[tmp$ol2, ]$GeneStableID
  biomart.qc$GeneOverlapRefSeq <- biomart.qc[tmp$ol2, ]$GeneName

  biomart.qc$intergenic <-
    ifelse(
      biomart.qc$GeneStableID == biomart.qc$GeneOverlapEnsemble &
        biomart.qc$minDist2TSS > dist_intergenic, TRUE, FALSE)
```

``` r
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

  minReads = 5e5

  ggReads <- 
    ggplot(meta, aes(y=Reads, x='')) +
    geom_boxplot(outlier.shape = NA, lwd=0.5, fatten=1) +
    geom_jitter(position = position_jitter(0.25), shape=20, size=0.5, alpha=0.25) +
    scale_color_manual(aes(colour = "lightgrey")) +
    scale_y_log10(limits=c(1e5, 1e7), label=scientific_10) +
    annotation_logticks(sides='l', size=0.25) +
    ylab('Total reads (log10)') +
    geom_hline(yintercept=minReads, linetype="dashed", color = "red", size=0.5) +
    theme_man +
      theme(
      legend.position = "none",
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
```

``` r
# Filter for non-imprinted and autosomal genes
  tmp <- splitGenes_autosomeXimprint(rownames(counts.qc), biomart.qc, imprint)
    counts.qc <- counts.qc[tmp$autosomes, ]
    rpkm.qc <- rpkm.qc[tmp$autosomes, ]
    biomart.qc <- biomart.qc %>% filter(GeneStableID %in% tmp$autosomes)

# Select subset of genes
# coding / noncoding independent genes (promoters separated by at least 4kb)
  genes.sel <- list()
  genes.sel[[1]] <- (biomart.qc %>% filter(coding & independent))$GeneStableID
  genes.sel[[2]] <- (biomart.qc %>% filter(lncRNA & independent))$GeneStableID

  names(genes.sel) <- c('protein_coding', 'lncRNA')

# Expression limits for making fit (CV2 vs mean)
  mean.thr <- c(0.01, 100) 
  names(mean.thr) <- c('min', 'max')

# Gene stats
  # Coding
  geneStats.coding <-
    basic_summary_stats(rpkm, genes.sel[[1]]) %>%
    mutate(biotype=names(genes.sel)[1]) %>%
    arrange(mean)

# Noncoding
  geneStats.nc <-
    basic_summary_stats(rpkm, genes.sel[[2]]) %>%
    mutate(biotype=names(genes.sel)[2]) %>%
    arrange(mean)
  
  df.tmp <- rbind(geneStats.coding, geneStats.nc)

# Number of events
  n.coding <- paste0('coding (n=', length(genes.sel[['protein_coding']]), ')') 
  n.lncRNA <- paste0('lncRNA (n=', length(genes.sel[['lncRNA']]), ')') 

# Apply upper and lower expression limits for making rolling mean
  geneStats.coding.lmt <-
    geneStats.coding %>%
    filter(mean <= mean.thr['max'] & mean >= mean.thr['min']) %>%
    arrange(mean) %>%
    mutate(
      roll = rollmean(.$cv2, k=15, na.pad=TRUE)
    )
  
  geneStats.nc.lmt <-
    geneStats.nc %>%
    filter(mean <= mean.thr['max'] & mean >= mean.thr['min']) %>%
    arrange(mean) %>%
    mutate(
      roll = rollmean(.$cv2, k=15, na.pad=TRUE)
    )

# Plot CV2 to the mean with loess fit to the rolling mean
  gg_CV2mean <-
    ggplot(df.tmp, aes(x=mean, y=cv2, size = biotype, color = biotype, alpha = biotype)) +
    geom_point() +
    scale_x_log10(labels=comma) +
    scale_y_log10(labels=comma) +
    annotation_logticks(sides = "bl", size=0.25) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.lncRNA, n.coding)) +
    scale_alpha_manual(values=c(0.75, .35), guide = 'none') +
    scale_size_manual(values=c(0.15, .0005), guide = 'none') +
    xlab("Mean expression (RPKM)") + 
    ylab("CV2") +
    geom_vline(xintercept = mean.thr, linetype="dashed", color = "red", size=0.25) +
    geom_smooth(method='loess', level=0, data = geneStats.coding.lmt, aes(x=mean, y=cv2), inherit.aes = F, color='darkgreen', size=0.75) +
    geom_smooth(method='loess', level=0, data = geneStats.nc.lmt, aes(x=mean, y=cv2), inherit.aes = F, color=col_nc, size=0.75) +
    theme_man +
    theme(legend.position = c(0.8, 0.95),
          legend.title = element_blank(),
          legend.text=element_text(size=6))


# Density and boxplot and for mean expression levels
  median.coding <-
    df.tmp %>%
    filter(biotype %in% 'protein_coding') %>% 
    summarise(median = median(mean))

  median.lncRNA <-
    df.tmp %>%
    filter(biotype %in% 'lncRNA') %>% 
    summarise(median = median(mean))

  p_mean <- wilcox.test(geneStats.coding$mean, geneStats.nc$mean)$p.value

  # Density mean
  limX = c(0.0005, 2000)  # to make nice plot
  gg_mean_density <-
    ggplot(df.tmp, aes(x=mean, color = biotype)) +
    stat_density(geom="line") +
    scale_x_log10(labels=comma, limits=limX) +
    annotation_logticks(sides = "b", size=0.25) +
    scale_color_manual(values=c(col_nc, col_coding)) +
    xlab('Mean expression (RPKM)') +
    geom_vline(xintercept=median.coding[1,1], linetype="dashed", size=0.25, color=col_coding) +
    geom_vline(xintercept=median.lncRNA[1,1], linetype="dashed", size=0.25, color=col_nc) +
    theme_man +
    theme(
      legend.position = 'none')

  # Boxplot mean
  gg_mean_box <-
    ggplot(df.tmp, aes(y=mean, x=biotype, color = biotype)) +
    geom_boxplot(outlier.size = 0.2, lwd=0.25) + 
    scale_y_log10(limits=limX) +
    annotation_logticks(sides = "l", size=0.25) +
    scale_color_manual(values=c(col_nc, col_coding)) +
    ylab('RPKM') +
    xlab('') +
    expand_limits(x = -.75) +
    theme_man +
    theme(
      axis.ticks.x=element_blank(),
      legend.title = element_blank(),
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Density and boxplot and for CV2
  # Set limits to make nice plot
  tmp.lmt <- c(0.01, 1000)

  p_cv2 <- wilcox.test(geneStats.coding$cv2, geneStats.nc$cv2)$p.value
  
  gg_cv2 <- 
    ggplot(df.tmp, aes(y=cv2, x=biotype, color = biotype)) +
    geom_violin(lwd=0.25) +
    geom_boxplot(lwd=0.25, outlier.colour = NA, width=0.15, fatten=1) + 
    scale_y_log10(labels=comma, limits = tmp.lmt) +
    annotation_logticks(sides = "l", size=0.25) +
    scale_color_manual(values=c(col_nc, col_coding)) +
    ylab('CV2') +
    xlab("") +
    expand_limits(x = -0.01) +
    theme_man +
    theme(
      axis.ticks.x=element_blank(),
      legend.title = element_blank(),
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

ggarrange(ggReads, gg_cv2, gg_mean_box, gg_mean_density, gg_CV2mean, ncol=5, widths = c(0.5/4, 0.5/4, 1/4, 1/4, 1/4), heights = c(1, 0.5, 1, 1))
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](ExtendedDataFig2_HEK293_files/figure-markdown_github/CV2%20vs%20mean%20expression%20(RPKM)-1.png)

``` r
print(p_mean)
```

    ## [1] 3.991817e-275

``` r
print(p_cv2)
```

    ## [1] 0
