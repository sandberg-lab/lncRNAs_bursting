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
  
  umiCast_path = file.path(prj_path, 'data/ss3_n682_fibs_umiCountsCast.csv')
  umiC57_path = file.path(prj_path, 'data/ss3_n682_fibs_umiCountsC57.csv')
  
  castBurst_path <- file.path(prj_path, 'data/ss3_n682_fibs_bursting_cast.qc.rds')
  c57Burst_path <- file.path(prj_path, 'data/ss3_n682_fibs_bursting_c57.qc.rds')
  
  meta_path = file.path(prj_path, 'data/ss3_n682_fibs_meta.rds')
  biomart_path = file.path(prj_path, 'data_dump/ss3_fibs_biomart_wGeneAnnotations.rds')
  
  halfLives_path <- file.path(prj_path, 'data_dump/primaryFibs_decayRates.rds')
  
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
  library(grid)
  
  source(file.path(prj_path, '/R_libs/Rlib_basic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_bursting.R'))

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
    col_coding <- col_sel[1]
    col_nc <- col_sel[3]
```

``` r
# Load data
  umi <- readRDS(umi_path)
  meta <- readRDS(meta_path)
  biomart <- readRDS(biomart_path)
  
  umiCast <- read.csv(umiCast_path, header = TRUE, row.names = 1) 
  umiC57 <- read.csv(umiC57_path, header = TRUE, row.names = 1)
  
  burst.cast <- readRDS(castBurst_path) # postQC
  burst.c57 <- readRDS(c57Burst_path)   # postQC

  halfLives <- readRDS(halfLives_path)
  
  imprint <- read.csv(imprint_path)
  escape <- read.csv(escape_path, header=TRUE)
```

``` r
# Note: Already filtered for QC. No further filter should be needed

# QC parameters (genes)
  minUMI = 1
  minUMIcells = 5

  passed_genes <- filter_genes(umi, minUMI, minUMIcells)
  
    umi <- umi[passed_genes, ]

    burst.cast <- burst.cast %>% filter(GeneStableID %in% passed_genes) %>% filter(biotype %in% c('lncRNA', 'protein_coding'))
    burst.c57 <- burst.c57 %>% filter(GeneStableID %in% passed_genes) %>% filter(biotype %in% c('lncRNA', 'protein_coding'))

    biomart <- biomart %>% filter(GeneStableID %in% passed_genes) %>% filter(GeneType %in% c('lncRNA', 'protein_coding'))
      
# Non-imprinted autosomal genes only
  tmp <- splitGenes_autosomeXimprint(burst.cast$GeneStableID, biomart, imprint)
    burst.cast <- burst.cast %>% filter(GeneStableID %in% tmp$autosomes)
    
  tmp <- splitGenes_autosomeXimprint(burst.c57$GeneStableID, biomart, imprint)
    burst.c57 <- burst.c57 %>% filter(GeneStableID %in% tmp$autosomes)
    
    umiCast <- umiCast[burst.cast$GeneStableID, ]
    umiC57 <- umiC57[burst.c57$GeneStableID, ]
    
    biomart <- biomart %>% filter(GeneStableID %in% rownames(umiCast) | GeneStableID %in% rownames(umiC57))
```

``` r
# QC parameters
  maxHalfLife_h = 10
  maxBurstDuration_h = 72

# Cast
  burst.cast <- 
    burst.cast %>%
    mutate(
      k_h = halfLives[match(.$GeneStableID, halfLives$GeneStableID), ]$k,
      HalfLife_h = halfLives[match(.$GeneStableID, halfLives$GeneStableID), ]$HalfLife_h,
      bf.ci.width = abs(.$bfhigh / .$bflow),
      bs.ci.width = abs(.$bshigh / .$bslow)
      ) %>%
    mutate(
      timeToBurst = 1 / (.$kon * .$k_h)
    ) %>%
    mutate(
      timeToBurst.pass = ifelse(timeToBurst < maxBurstDuration_h & !is.na(timeToBurst), TRUE, FALSE),
      HalfLife_h.pass = ifelse(HalfLife_h < maxHalfLife_h & !is.na(HalfLife_h), TRUE, FALSE)
    )

# C57
  burst.c57 <- 
    burst.c57 %>%
    mutate(
      k_h = halfLives[match(.$GeneStableID, halfLives$GeneStableID), ]$k,
      HalfLife_h = halfLives[match(.$GeneStableID, halfLives$GeneStableID), ]$HalfLife_h,
      bf.ci.width = abs(.$bfhigh / .$bflow),
      bs.ci.width = abs(.$bshigh / .$bslow)
      ) %>%
    mutate(
      timeToBurst = 1 / (.$kon * .$k_h)
    ) %>%
    mutate(
      timeToBurst.pass = ifelse(timeToBurst < maxBurstDuration_h & !is.na(timeToBurst), TRUE, FALSE),
      HalfLife_h.pass = ifelse(HalfLife_h < maxHalfLife_h & !is.na(HalfLife_h), TRUE, FALSE)
    )
```

``` r
# main Fig 2A-C

# QC parameters
  thr.size = c(0.2, 50)
  thr.kon = c(0.01, 30)
  thr.mean = c(0.01, 100)


# Make df with genes passing QC for both Cast and C57
  cast.tmp <-
    burst.cast %>%
    select(GeneStableID, kon, size, mean, biotype) %>%
    filter(
      GeneStableID %in% burst.c57$GeneStableID & 
      biotype %in% c('protein_coding', 'lncRNA')
    )
  
  c57.tmp <-
    burst.c57 %>%
    select(GeneStableID, kon, size, mean, biotype) %>%
    filter(
      GeneStableID %in% burst.cast$GeneStableID &
      biotype %in% c('protein_coding', 'lncRNA')
    )

  colnames(cast.tmp) <- c('GeneStableID', 'kon_cast', 'size_cast', 'mean_cast', 'biotype')
  colnames(c57.tmp) <- c('GeneStableID_c57', 'kon_c57', 'size_c57', 'mean_c57', 'biotype_c57')

  df.2alleles.filt <- cbind(cast.tmp, c57.tmp)

# Correlations Cast:C57

  # Coding
  cor_coding_freq = 
    round(cor(
      (df.2alleles.filt %>% filter(biotype %in% 'protein_coding'))$kon_cast,
      (df.2alleles.filt %>% filter(biotype %in% 'protein_coding'))$kon_c57,
      method=c('spearman')), digits = 2)
  cor_coding_size = 
    round(cor(
      (df.2alleles.filt %>% filter(biotype %in% 'protein_coding'))$size_cast,
      (df.2alleles.filt %>% filter(biotype %in% 'protein_coding'))$size_c57,
      method=c('spearman')), digits = 2)
  cor_coding_mean = 
    round(cor(
      (df.2alleles.filt %>% filter(biotype %in% 'protein_coding'))$mean_cast,
      (df.2alleles.filt %>% filter(biotype %in% 'protein_coding'))$mean_c57,
      method=c('spearman')), digits = 2)
  
  # lncRNA
  cor_lncRNA_freq = 
    round(cor(
      (df.2alleles.filt %>% filter(biotype %in% 'lncRNA'))$kon_cast,
      (df.2alleles.filt %>% filter(biotype %in% 'lncRNA'))$kon_c57,
      method=c('spearman')), digits = 2)
  cor_lncRNA_size = 
    round(cor(
      (df.2alleles.filt %>% filter(biotype %in% 'lncRNA'))$size_cast,
      (df.2alleles.filt %>% filter(biotype %in% 'lncRNA'))$size_c57,
      method=c('spearman')), digits = 2)
  cor_lncRNA_mean = 
    round(cor(
      (df.2alleles.filt %>% filter(biotype %in% 'lncRNA'))$mean_cast,
      (df.2alleles.filt %>% filter(biotype %in% 'lncRNA'))$mean_c57,
      method=c('spearman')), digits = 2)
  
# Legends (correlation in plot)
  # Frequency
  txt_freq_nc <- grobTree(textGrob(
    paste0("rho (lncRNA) = ", cor_lncRNA_freq),
    x=0.1,  y=0.95, hjust=0,
    gp=gpar(col=col_nc, fontsize=8)))
  txt_freq_coding <- grobTree(textGrob(
    paste0("rho (coding) = ", cor_coding_freq),
    x=0.1,  y=0.85, hjust=0,
    gp=gpar(col=col_coding, fontsize=8)))
  
  # Size
  txt_size_nc <- grobTree(textGrob(
    paste0("rho (lncRNA) = ", cor_lncRNA_size),
    x=0.1,  y=0.95, hjust=0,
    gp=gpar(col=col_nc, fontsize=6)))
  txt_size_coding <- grobTree(textGrob(
    paste0("rho (coding) = ", cor_coding_size),
    x=0.1,  y=0.85, hjust=0,
    gp=gpar(col=col_coding, fontsize=8)))

  # Mean
  txt_mean_nc <- grobTree(textGrob(
    paste0("rho (lncRNA) = ", cor_lncRNA_mean),
    x=0.1,  y=0.95, hjust=0,
    gp=gpar(col=col_nc, fontsize=8)))
  txt_mean_coding <- grobTree(textGrob(
    paste0("rho (coding) = ", cor_coding_mean),
    x=0.1,  y=0.85, hjust=0,
    gp=gpar(col=col_coding, fontsize=8)))
  
  # Number of genes
  n.coding <- table(df.2alleles.filt$biotype)[2]
  n.lncRNA <- table(df.2alleles.filt$biotype)[1]


  gg_main2A <- 
    ggplot(df.2alleles.filt, aes(x=kon_cast, y=kon_c57, color=biotype, size=biotype, alpha=biotype)) +
    geom_point(shape=16) +
    scale_x_log10(labels=comma, limit=thr.kon) +
    scale_y_log10(labels=comma, limit=thr.kon) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.lncRNA, n.coding)) +
    scale_size_manual(values = c(0.4, 0.25)) +
    scale_alpha_manual(values = c(0.95, 0.5)) +
    annotation_logticks(sides = "bl", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    geom_abline(intercept = 0, slope = 1, linetype="dashed", color = 'red', size=0.5) +
    xlab('Burst frequency (CAST)') +
    ylab('Burst frequency (C57)') +
    annotation_custom(txt_freq_nc) +
    annotation_custom(txt_freq_coding) +
    theme_man +
    theme(legend.position = 'none')
  
  gg_main2B <- 
    ggplot(df.2alleles.filt, aes(x=size_cast, y=size_c57, color=biotype, size=biotype, alpha=biotype)) +
    geom_point(shape=16) +
    scale_x_log10(labels=comma, limit=thr.size) +
    scale_y_log10(labels=comma, limit=thr.size) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.lncRNA, n.coding)) +
    scale_size_manual(values = c(0.4, 0.25)) +
    scale_alpha_manual(values = c(0.95, 0.5)) +
    annotation_logticks(sides = "bl", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    geom_abline(intercept = 0, slope = 1, linetype="dashed", color = 'red', size=0.5) +
    xlab('Burst size (CAST)') +
    ylab('Burst size (C57)') +
    annotation_custom(txt_size_nc) +
    annotation_custom(txt_size_coding) +
    theme_man +
    theme(legend.position = 'none')
  
  gg_main2C <- 
    ggplot(df.2alleles.filt, aes(x=mean_cast, y=mean_c57, color=biotype, size=biotype, alpha=biotype)) +
    geom_point(shape=16) +
    scale_x_log10(labels=comma, limit=thr.mean) +
    scale_y_log10(labels=comma, limit=thr.mean) +
    scale_size_manual(values = c(0.4, 0.25)) +
    scale_alpha_manual(values = c(0.95, 0.5)) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.lncRNA, n.coding)) +
    annotation_logticks(sides = "bl", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    geom_abline(intercept = 0, slope = 1, linetype="dashed", color = 'red', size=0.5) +
    xlab('Mean (CAST)') +
    ylab('Mean (C57)') +
    annotation_custom(txt_mean_nc) +
    annotation_custom(txt_mean_coding) +
    theme_man +
    theme(legend.position = c(0.8, 0.2),
          legend.key.size = unit(0.1, "cm"),
          legend.key.width = unit(0.1, "cm"))
  

ggarrange(gg_main2A, gg_main2B, gg_main2C, ncol = 3, nrow = 1)
```

![](3_ss3_n682_burstKinetics_lncRNAs2mRNAs_files/figure-markdown_github/Burst%20kinetics%20correlations%20CAST%20vs%20C57%20(lncRNAs%20vs%20mRNAs)-1.png)

``` r
# Generate plots
  c57_plots <- 
    burst.c57 %>%
    filter(independent) %>%
    filter(biotype %in% c('protein_coding', 'lncRNA')) %>%
    plot_coding2noncoding_density(biopmart, 1000, title = 'C57')

  cast_plots <- 
    burst.cast %>%
    filter(independent) %>%
    filter(biotype %in% c('protein_coding', 'lncRNA')) %>%
    plot_coding2noncoding_density(biopmart, 1000, title='CAST')


# Figs main 2d-f (C57), Extended Data Fig 4a-c (CAST) 
ggarrange(
  c57_plots[[1]], c57_plots[[2]], c57_plots[[3]],
  cast_plots[[1]], cast_plots[[2]], cast_plots[[3]],
  ncol = 3, nrow = 2
  )
```

![](3_ss3_n682_burstKinetics_lncRNAs2mRNAs_files/figure-markdown_github/burst%20size%20and%20frequencies%20for%20lncRNA%20and%20mRNAs-1.png)

``` r
# Main Fig 2g and Extended Data Fig. 4f

# Plot time to next burst; coding vs non-coding
  binSize_burstDuration = 3

  ggC57 <- 
    burst.c57 %>% filter(biotype %in% c('lncRNA', 'protein_coding')) %>%
    plot_burstDuration(binSize=binSize_burstDuration, xLim = 72, label='CAST')

  ggCast <- 
    burst.cast %>% filter(biotype %in% c('lncRNA', 'protein_coding')) %>%
    plot_burstDuration(binSize=binSize_burstDuration, xLim = 72, label='CAST')

  
ggarrange(ggC57, ggCast)
```

![](3_ss3_n682_burstKinetics_lncRNAs2mRNAs_files/figure-markdown_github/burst%20duration-1.png)

``` r
# Gene Stats
  # CAST
  genes_sel <- burst.cast %>% filter(independent) %>% filter(coding | lncRNA)
  tmp <- parallel::mclapply(genes_sel$GeneStableID, get_geneStats, umiCast, mc.cores = 15)
  geneStats_umiCast <-
    as.data.frame(do.call(rbind, tmp)) %>%
    mutate(
      GeneStableID = genes_sel$GeneStableID,
      biotype = genes_sel$biotype
      )

  # C57
  genes_sel <- burst.c57 %>% filter(independent) %>% filter(coding | lncRNA)
  tmp <- parallel::mclapply(genes_sel$GeneStableID, get_geneStats, umiC57, mc.cores = 15)
  geneStats_umiC57 <-
    as.data.frame(do.call(rbind, tmp)) %>%
    mutate(
      GeneStableID = genes_sel$GeneStableID,
      biotype = genes_sel$biotype
      )
    
# Rank cv2 of genes
  nTop = 50
  
  # CAST
  lncRNA_ranks_Cast <- 
    data.frame(
      GeneStableID = (burst.cast %>% filter(independent & lncRNA))$GeneStableID
    ) %>%
    mutate(
      rankingScore = sapply(.$GeneStableID, get_cv2rank, (burst.cast %>% filter(coding))$GeneStableID, geneStats_umiCast)
    ) %>%
    arrange(rankingScore) %>%
    mutate(
      Rank = seq(nrow(.))
    ) %>%
    mutate(
      highRank = ifelse(.$Rank <= nTop, TRUE, FALSE)
    )
  
  # C57
  lncRNA_ranks_C57 <- 
    data.frame(
      GeneStableID = (burst.c57 %>% filter(independent & lncRNA))$GeneStableID
    ) %>%
    mutate(
      rankingScore = sapply(.$GeneStableID, get_cv2rank,  (burst.c57 %>% filter(coding))$GeneStableID, geneStats_umiC57)
    ) %>%
    arrange(rankingScore) %>%
    mutate(
      Rank = seq(nrow(.))
    ) %>%
    mutate(
      highRank = ifelse(.$Rank <= nTop, TRUE, FALSE)
    )
```

``` r
# Plot top ranked genes
  gg_C57 <- plot_topRanks(geneStats_umiC57, lncRNA_ranks_C57, label='C57')
  gg_CAST <- plot_topRanks(geneStats_umiCast, lncRNA_ranks_Cast, label='CAST')

  ggarrange(gg_C57, gg_CAST)
```

![](3_ss3_n682_burstKinetics_lncRNAs2mRNAs_files/figure-markdown_github/Plot%20CV2%20to%20mean%20expression%20and%20mark%20lncRNAs%20with%20highly%20ranked%20CV2-1.png)

``` r
# Bursting stats (medians / means) for lncRNA
  # CAST
    # lncRNA medians / mean
    topGenes <- lncRNA_ranks_Cast %>% filter(highRank)
    cast_lncRNAs_median <- get_medians(topGenes$GeneStableID, burst.cast)
    
    # Match expression
    tmp.l <- 
      parallel::mclapply(
        topGenes$GeneStableID,
        match_mean,
        geneStats = geneStats_umiCast,
        subset_genes = (geneStats_umiCast %>% filter(biotype %in% 'protein_coding'))$GeneStableID,
        n_matched = 10,
        mc.cores = 20)
    
    names(tmp.l) <- topGenes$GeneStableID
   
    matched_cast_long <- 
      as.data.frame(do.call(cbind, tmp.l)) %>%
      gather('geneToMatch', 'geneMatched', 1:length(tmp.l))

  # C57
    # lncRNA medians / mean
    topGenes <- lncRNA_ranks_C57 %>% filter(highRank)
    c57_lncRNAs_median <- get_medians(topGenes$GeneStableID, burst.c57)
    
    # Match expression
    tmp.l <- 
      parallel::mclapply(
        topGenes$GeneStableID,
        match_mean,
        geneStats = geneStats_umiC57,
        subset_genes = (geneStats_umiC57 %>% filter(biotype %in% 'protein_coding'))$GeneStableID,
        n_matched = 10,
        mc.cores = 20)
    
    names(tmp.l) <- topGenes$GeneStableID
   
    matched_c57_long <- 
      as.data.frame(do.call(cbind, tmp.l)) %>%
      gather('geneToMatch', 'geneMatched', 1:length(tmp.l))

  
# Number of permutations
  nPerm = 10000

  matched_cast.l <-
    parallel::mclapply(
      seq(nPerm),
      sample_matchedMean_bursting,
      matched_Genes = matched_cast_long,
      burst_stats = burst.cast %>% select(GeneStableID, kon, size, mean),
      mc.cores = 20
      )
  matched_cast.df <- do.call(rbind, matched_cast.l)    
  
  matched_c57.l <-
    parallel::mclapply(
      seq(nPerm),
      sample_matchedMean_bursting,
      matched_Genes = matched_c57_long,
      burst_stats = burst.c57 %>% select(GeneStableID, kon, size, mean),
      mc.cores = 20
      )
  matched_c57.df <- do.call(rbind, matched_c57.l)

# Make plots (main Figs 2H-I and Extended Data Fig 4I,J)
  plot_C57 <- makePlots(matched_c57.df, c57_lncRNAs_median, nPerm=nPerm, 'C57')
  plot_Cast <- makePlots(matched_cast.df, cast_lncRNAs_median, nPerm=nPerm, 'CAST')

    
  ggarrange(
    plot_C57[[1]], plot_C57[[2]], plot_C57[[3]],
    plot_Cast[[1]], plot_Cast[[2]], plot_Cast[[3]],
    ncol = 3, nrow = 2)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](3_ss3_n682_burstKinetics_lncRNAs2mRNAs_files/figure-markdown_github/Permutation%20test%20for%20burst%20kinetics%20comparing%20lncRNA%20bursting%20to%20sxpression%20matched%20mRNAs-1.png)
