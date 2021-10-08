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
  meta_path = file.path(prj_path, 'data_dump/ss2_fibs_n533_meta.rds')
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
  source(file.path(prj_path, '/R_libs/Rlib_matchMean_rankVar.R'))

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
# Load data
  rpkm <- readRDS(rpkm_path)
  counts <- readRDS(counts_path)
  meta <- readRDS(meta_path)
  biomart <- readRDS(biomart_path)

  imprint <- read.csv(imprint_path)
  escape <- read.csv(escape_path, header=TRUE)
```

``` r
# QC genes (note that data input is already QC)
  min.counts = 3
  min.cells = 5
  
  passed_genes <- filter_genes(counts, min.counts, min.cells)
    counts_qc <- counts[passed_genes, ]
    rpkm_qc <- rpkm[passed_genes, ]
    biomart_qc <- biomart %>% filter(GeneStableID %in% passed_genes)

# Filter for non-imprinted and autosomal genes
  tmp <- splitGenes_autosomeXimprint(rownames(counts_qc), biomart, imprint)
    counts_qc <- counts_qc[tmp$autosomes, ]
    rpkm_qc <- rpkm_qc[tmp$autosomes, ]
    biomart_qc <- biomart_qc %>% filter(GeneStableID %in% tmp$autosomes)

# Select subset of genes
# coding / noncoding independent genes (promoters separated by at least 4kb)
  genes.sel <- list()
  genes.sel[[1]] <- (biomart_qc %>% filter(coding & independent))$GeneStableID
  genes.sel[[2]] <- (biomart_qc %>% filter(lncRNA & independent))$GeneStableID

  names(genes.sel) <- c('protein_coding', 'lncRNA')
```

``` r
# Expression limits for making fit (CV2 vs mean)
  mean.thr <- c(0.001, 1) 
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
  gg_main1E <-
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
  limX = c(0.001, 1)  # to make nice plot
  gg_main1B_density <-
    ggplot(df.tmp, aes(x=mean, color = biotype)) +
    stat_density(geom="line") +
    scale_x_log10(limits=limX) +
    annotation_logticks(sides = "b", size=0.25) +
    scale_color_manual(values=c(col_nc, col_coding)) +
    xlab('Mean expression (RPKM)') +
    geom_vline(xintercept=median.coding[1,1], linetype="dashed", size=0.25, color=col_coding) +
    geom_vline(xintercept=median.lncRNA[1,1], linetype="dashed", size=0.25, color=col_nc) +
    theme_man +
    theme(
      legend.position = 'none')

  # Boxplot mean
  gg_main1B_box <-
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
  
  gg_main1D <- 
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

ggarrange(gg_main1E, gg_main1B_density, gg_main1B_box, gg_main1D, ncol=4, widths = c(1/3, 1/3, 0.5/3, 0.5/3), heights = c(1, 0.5, 1, 1))
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning: Removed 4037 rows containing non-finite values (stat_density).

    ## Warning: Removed 4037 rows containing non-finite values (stat_boxplot).

![](3_cv2_mRNA.vs.lncRNA_SS2_mouseFibroblasts_files/figure-markdown_github/CV2%20vs%20mean%20expression%20(RPKM)-1.png)

``` r
print(p_mean)
```

    ## [1] 2.718409e-282

``` r
print(p_cv2)
```

    ## [1] 3.80927e-293

``` r
# Set Basic parameters for ranking
  # Expression limits (lncRNAs) for ranking test (Low expressed genes have few other genes to be matched against)
  mean.perm.thr <- c(0.001, 100) 
  names(mean.perm.thr) <- c('min', 'max')
  
  # Number of genes to Rank against
  nRank = 100

# Filter for expression limits
  geneStats.nc.lmt <-
    geneStats.nc %>%
    filter(
      mean <= mean.perm.thr['max'] &
      mean >= mean.perm.thr['min']
      )
  
  geneStats.coding.lmt <-
    geneStats.coding %>%
    filter(
      mean <= mean.perm.thr['max'] &
      mean >= mean.perm.thr['min']
      )
  
# Rank cv2 of lncRNAs to 100 mRNAs of similar mean expression levels
  # Note: Using lncRNA (within min/max limits) to match match with all coding --> indirect limits also for coding
  
  # Gene stats for coding and noncoding into one df
  geneStats_nc.coding <- rbind(geneStats.nc, geneStats.coding)
  
  # Find genes (mRNAs) of most similar mean expression
  matched_lncRNA2coding =
    lapply(
      geneStats.nc.lmt$GeneStableID,
      match_mean,
      geneStats=geneStats_nc.coding,
      subset_genes=geneStats.coding$GeneStableID,
      n_matched=nRank)
  
  names(matched_lncRNA2coding) <- geneStats.nc.lmt$GeneStableID
  
  # Rank cv2 of lncRNAs to matched genes
  df.rank.lncRNA <- 
    do.call(
      rbind.data.frame,
      lapply(
        names(matched_lncRNA2coding),
        rank_var,
        geneStats=geneStats_nc.coding,
        exprMatchedGenes=matched_lncRNA2coding))

# Rank CV2 of randomly chosen mRNAs (background test)

  # Match all mRNAs to other mRNAs of similar mean expression
  matched_coding2coding <-
    lapply(
      geneStats.coding.lmt$GeneStableID,
      match_mean, 
      geneStats = geneStats.coding,
      subset_genes = geneStats.coding$GeneStableID,
      n_matched = nRank
    )
  names(matched_coding2coding) <- geneStats.coding.lmt$GeneStableID

  # Randomly sample mRNAs (equal number of lncRNAs) and rank to other mRNAs
  nSampling = 100 # Number of sampling (of coding genes)

  perm_coding <- 
    parallel::mclapply(
      seq(nSampling),
      rank_var_perm,
      geneStats=geneStats.coding,
      exprMatchedGenes = matched_coding2coding,
      nGenes_sample=nrow(geneStats.nc.lmt),
      mc.cores = 20)
  
  perm_coding.df <- do.call(rbind.data.frame, perm_coding)

  gg_main1G <- 
    ggplot(perm_coding.df, aes(x=rank, group=i)) +
    geom_density(size=0.05, colour=col_coding, alpha=0.25) +
    xlim(c(0,1)) +
    xlab('Rank: cv2(lncRNA) > cv2(coding)') +
    geom_density(data=df.rank.lncRNA, aes(x=rank), colour=col_nc, size=0.25, inherit.aes = FALSE) +
    geom_vline(xintercept=median(df.rank.lncRNA$rank), linetype="dashed", color = col_nc, size=0.25) +
    geom_vline(xintercept=median(perm_coding.df$rank), linetype="dashed", color = col_coding, size=0.25) +
    theme_man +
    theme(legend.position = 'none')

  plot(gg_main1G)
```

![](3_cv2_mRNA.vs.lncRNA_SS2_mouseFibroblasts_files/figure-markdown_github/match%20each%20lncRNA%20to%20100%20mRNAs%20of%20most%20similar%20mean%20expression%20and%20rank%20CV2-1.png)

``` r
# Basic input parameters
  # Permutations
  nPerm = 10000

  # Gene state for coding and noncoding genes (calculated above)
  geneStats_nc.coding <- rbind(geneStats.nc, geneStats.coding)

# Match mean expression of lncRNA to mRNAs
  matched_lncRNA2coding.l <-
    lapply(
      geneStats.nc.lmt$GeneStableID,
      match_mean,
      geneStats_nc.coding,
      geneStats.coding$GeneStableID,
      10)
  
  names(matched_lncRNA2coding.l) <- geneStats.nc.lmt$GeneStableID
  
# df (long) for downstream analysis
  matched_lncRNA2coding.df <-
    as.data.frame(do.call(cbind, matched_lncRNA2coding.l)) %>%
    gather('geneToMatch', 'geneMatched', 1:length(matched_lncRNA2coding.l))

# Make permutation test
  # w/wo replacement 
  # Note wo replacement: many permutations are needed to get successful matching for all lncRNA genes
  
  permTest <- 
    parallel::mclapply(
    seq(nPerm),
    sample_matchedMean,
    matched.df = matched_lncRNA2coding.df,
    expressionStats = geneStats.coding,
    wReplacement = TRUE,
    mc.cores=20)
  
  permTest.df <- do.call(rbind, permTest)
  
# Plot CV2 of expression matched mRNAs (contrasted lncRNAs)
  p_frequency <- length(which(permTest.df$cv2 > median(geneStats.nc.lmt$cv2))) / nPerm
    
  gg_main1F <- 
    ggplot(permTest.df, aes(cv2)) +
    geom_histogram(fill=col_coding, colour='lightgrey', size=0.1) +
    xlab('CV2 (median)') +
    ylab('Freq') +
    geom_vline(xintercept=median(geneStats.nc.lmt$cv2), linetype="dashed", color = col_nc, size=0.25) +
    theme_man
  
  gg_control <- 
    ggplot(permTest.df, aes(meanExp_mean)) +
    geom_histogram(fill=col_coding, colour='lightgrey', size=0.1) +
    scale_x_log10() +
    xlab('RPKM (mean)') +
    ylab('Freq') +
    geom_vline(xintercept=mean(geneStats.nc.lmt$mean), linetype="dashed", color = col_nc, size=0.25) +
    theme_man

  
  
  ggarrange(gg_main1F, gg_control, ncol = 2, nrow = 1)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](3_cv2_mRNA.vs.lncRNA_SS2_mouseFibroblasts_files/figure-markdown_github/permutation%20test%20of%20CV2%20(lncRNA%20vs%20mRNAs)%20for%20expression%20matched%20mRNAs-1.png)

``` r
###########
# Some control plots for the matched expression

        
# Some representative plots of randomly expression matched mRNAs (vs lncRNA)
# Comparing the distribution of mean expression of lncRNA and expression matched mRNAs

  fncToPlot <- function(matched.df, geneStats, geneStats_lncRNAs){
    matched_genes_sampled <- 
      (matched.df %>%
      group_by(geneToMatch) %>%
      sample_n(1))$geneMatched
    
    stats_sampled <- geneStats[match(matched_genes_sampled, geneStats$GeneStableID), ]
    
    tmp <- rbind(stats_sampled, geneStats_lncRNAs)
    
    gg <- 
      ggplot(tmp, aes(mean, color=biotype)) +
      geom_density() +
      scale_x_log10() + 
      annotation_logticks(sides = 'b') +
      theme_man
    
    return(gg)

  }
  
  gg1 <- fncToPlot(matched_lncRNA2coding.df, geneStats.coding, geneStats.nc.lmt)
  gg2 <- fncToPlot(matched_lncRNA2coding.df, geneStats.coding, geneStats.nc.lmt)
  gg3 <- fncToPlot(matched_lncRNA2coding.df, geneStats.coding, geneStats.nc.lmt)
  gg4 <- fncToPlot(matched_lncRNA2coding.df, geneStats.coding, geneStats.nc.lmt)
  gg5 <- fncToPlot(matched_lncRNA2coding.df, geneStats.coding, geneStats.nc.lmt)
  gg6 <- fncToPlot(matched_lncRNA2coding.df, geneStats.coding, geneStats.nc.lmt)

  ggarrange(gg1, gg2, gg3, gg4, gg5, gg6, ncol = 3, nrow = 2)
```

![](3_cv2_mRNA.vs.lncRNA_SS2_mouseFibroblasts_files/figure-markdown_github/Some%20more%20control%20plots%20for%20the%20expression%20matched%20mRNAs-1.png)
