### Functions used for bursting

### Gene stats (mean, cv2, var) for allele-specific UMIs
get_geneStats <- function(gene, expr) {
  tmp <- expr[gene, ]
  tmp <- tmp[!is.na(tmp)]
  tmp[which(tmp==0)] <- 0
  
  gene.mean = mean(tmp)
  gene.var = var(tmp)
  gene.cv2 = gene.var / gene.mean^2
  geneStats <- c(gene.mean, gene.var, gene.cv2)
  names(geneStats) = c('mean', 'var', 'cv2')
  return(geneStats)  
}  


### Rank CV2 (i.e. of lncRNAs)
get_cv2rank <- function(geneToRank, genesToRankAgainst, geneStats, nGenes_useToRank=100, minMean=0.1){
  
  stats_geneToRank <- geneStats %>% filter(GeneStableID %in% geneToRank)
  
  if(stats_geneToRank$mean >= minMean){
    stats_genesToRankAgainst <-
      geneStats %>%
      filter(GeneStableID %in% genesToRankAgainst) %>%
      mutate(
        diff = .$mean - stats_geneToRank$mean
      )
    
    tmp.df <- 
      rbind(
        (stats_genesToRankAgainst %>% 
           filter(diff >= 0) %>%
           arrange(diff)
        )[1:(nGenes_useToRank/2), ],
        (stats_genesToRankAgainst %>% 
           filter(diff < 0) %>%
           arrange(abs(diff))
        )[1:(nGenes_useToRank/2), ]
      )
    return(length(which(tmp.df$cv2 > stats_geneToRank$cv2)))
  }
  
  if(stats_geneToRank$mean < minMean){
    return(NA)
  }
}

### Medians for bursting parameters
get_medians <- function(subsetGenes, burstingParameters){
  df <-
    data.frame(
      kon = median((burstingParameters %>% filter(GeneStableID %in% subsetGenes))$kon),
      size = median((burstingParameters %>% filter(GeneStableID %in% subsetGenes))$size),
      meanExp = median((burstingParameters %>% filter(GeneStableID %in% subsetGenes))$mean)
    )
  return(df)
}


### Find genes with similar mean expression
match_mean <- function(gene, geneStats, subset_genes, n_matched) {
  # gene: selected gene (to be matched)
  # geneStats: df with (mean and GeneStableID)
  # subset_genes: subset of genes to be used for matching
  # n_matched: number of matched genes to output
  
  gene_sel <- geneStats %>% filter(GeneStableID %in% gene)
  
  # Select genes and exclude gene to be matched
  geneStats_sel <-
    geneStats %>%
    filter(
      GeneStableID %in% subset_genes &
        GeneStableID %notin% gene) %>%
    mutate(diff = .$mean - gene_sel$mean) %>%
    arrange(abs(diff))
  
  # Dump matcehd genes
  return(geneStats_sel$GeneStableID[1:n_matched])
}


### Sample expression matched genes
sample_matchedMean_bursting <- function(nPerm, matched_Genes, burst_stats){
  # nPerm: Number of permutations
  # matched_Genes: matched genes (lncRNA), output from match_mean in long df format
  # burst_stats: bursting parameters with (mean, kon, size, GeneStableID)
  `%notin%` <- Negate(`%in%`)
  
  # Draw one gene for each öncRNA
  matched_genes_sampled <- 
    (matched_Genes %>%
       group_by(geneToMatch) %>%
       sample_n(1))$geneMatched
  
  # Get bursting parameters for sampled genes and claculate the medians
  burst_sampled <- burst_stats[match(matched_genes_sampled, burst_stats$GeneStableID), ]
  
  df.out <- 
    data.frame(
      kon = median(burst_sampled$kon),
      size = median(burst_sampled$size),
      meanExp_median = median(burst_sampled$mean)
    )
  # Dump
  return(df.out)
}



###########################
### To make nice plots
###########################

# Plot burst size, burst freq and mean expression (for allele distributed UMIs)
plot_coding2noncoding_density <- function(df, biomart, UMImax, title) {
  df <- df %>% filter(mean < UMImax)
  
  # Coding vs lncRNA
  df.tmp1 <- df %>% filter(biotype %in% 'protein_coding')
  df.tmp2 <- df %>% filter(biotype %in% 'lncRNA')
  
  # n.biotype
  n.coding <- nrow(df.tmp1)
  n.lncRNA <- nrow(df.tmp2)
  
  # p-values
  pSize = signif(wilcox.test(df.tmp1$size, df.tmp2$size)$p.value, 3)
  pFreq = signif(wilcox.test(df.tmp1$kon, df.tmp2$kon)$p.value, 3)
  pMean = signif(wilcox.test(df.tmp1$mean, df.tmp2$mean)$p.value, 3)
  
  # Legends for plot
  pSize <- 
    grobTree(textGrob(
      paste0('p: ', pSize),
      x=0.5,  y=0.95, hjust=0,
      gp=gpar(col='black', fontsize=8)))
  
  pFreq <- 
    grobTree(textGrob(
      paste0('p: ', pFreq),
      x=0.5,  y=0.95, hjust=0,
      gp=gpar(col='black', fontsize=8)))
  
  pMean <- 
    grobTree(textGrob(
      paste0('p: ', pMean),
      x=0.5,  y=0.95, hjust=0,
      gp=gpar(col='black', fontsize=8)))
  
  nCoding <- 
    grobTree(textGrob(
      paste0('n: ', n.coding),
      x=0.1,  y=0.95, hjust=0,
      gp=gpar(col=col_coding, fontsize=8)))
  
  nlncRNA <- 
    grobTree(textGrob(
      paste0('n: ', n.lncRNA),
      x=0.1,  y=0.85, hjust=0,
      gp=gpar(col=col_nc, fontsize=8)))
  
  # Stats
  bs_coding <- df.tmp1 %>% summarise(median(size))
  bs_lncRNA <- df.tmp2 %>% summarise(median(size))
  bf_coding <- df.tmp1 %>% summarise(median(kon))
  bf_lncRNA <- df.tmp2 %>% summarise(median(kon))
  median_coding <- df.tmp1 %>% summarise(median(mean))
  median_lncRNA <- df.tmp2 %>% summarise(median(mean))
  
  fold_mean = signif(median_coding / median_lncRNA, 3)
  fold_size = signif(bs_coding / bs_lncRNA, 3)
  fold_freq =  signif(bf_coding / bf_lncRNA, 3)
  
  # Plots
  
  dump.plots = list()
  
  dump.plots[[1]] <-
    ggplot(df, aes(x=kon, color=biotype)) +
    geom_density() +
    scale_x_log10(labels=comma, limits=c(0.01, 50)) +
    annotation_logticks(sides = "b", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.coding, n.lncRNA)) +
    geom_vline(xintercept=bf_coding[1,1], linetype="dashed", size=0.25, color=col_coding) +
    geom_vline(xintercept=bf_lncRNA[1,1], linetype="dashed", size=0.25, color=col_nc) +
    xlab("Burst frequency") + 
    annotation_custom(pFreq) +
    ggtitle(paste0(title, ' ', 'fold: ', fold_freq)) +
    theme_man +
    theme(legend.position = 'none')
  
  
  dump.plots[[2]] <-
    ggplot(df, aes(x=size, color=biotype)) +
    geom_density() +
    scale_x_log10(labels=comma, limits=c(0.1, 100)) +
    annotation_logticks(sides = "b", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.coding, n.lncRNA)) +
    geom_vline(xintercept=bs_coding[1,1], linetype="dashed", size=0.25, color=col_coding) +
    geom_vline(xintercept=bs_lncRNA[1,1], linetype="dashed", size=0.25, color=col_nc) +
    xlab("Burst size") + 
    annotation_custom(pSize) +
    annotation_custom(nCoding) +
    annotation_custom(nlncRNA) +
    ggtitle(paste0(title, ' ', 'fold: ', fold_size)) +
    theme_man +
    theme(legend.position = 'none')
  
  dump.plots[[3]] <-
    ggplot(df, aes(x=mean, color=biotype)) +
    geom_density() +
    scale_x_log10(labels=comma, limits=c(0.01, 200)) +
    annotation_logticks(sides = "b", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    scale_color_manual(values=c(col_nc, col_coding), labels = c(n.coding, n.lncRNA)) +
    geom_vline(xintercept=median_coding[1,1], linetype="dashed", size=0.25, color=col_coding) +
    geom_vline(xintercept=median_lncRNA[1,1], linetype="dashed", size=0.25, color=col_nc) +
    xlab("Mean expression") + 
    annotation_custom(pMean) +
    ggtitle(paste0(title, ' ', 'fold: ', fold_mean)) +
    theme_man +
    theme(legend.position = 'none')
  
  names(dump.plots) <- c('size', 'freq', 'mean')
  return(dump.plots)
}  


### Plot duration between two bursts on the same allele
plot_burstDuration <-
  function(
    burstInference,
    binSize = binSize_burstDuration,
    xLim = 72,
    label
  ){
    
    # Filter genes to plot
    burstInference <- 
      burstInference %>%
      filter(
        !is.na(timeToBurst) &
          timeToBurst.pass, HalfLife_h.pass
      )
    
    # Set labels 
    n.coding <- paste0('coding (n=', nrow(burstInference %>% filter(biotype %in% 'protein_coding')), ')') 
    n.lncRNA <- paste0('lncRNA (n=', nrow(burstInference %>% filter(biotype %in% 'lncRNA')), ')') 
    
    # Median
    medCoding <- median((burstInference %>% filter(biotype %in% 'protein_coding'))$timeToBurst)
    medNC <- median((burstInference %>% filter(biotype %in% 'lncRNA'))$timeToBurst)
    
    tmp <- seq(0, xLim, binSize)
    gg <- 
      ggplot(burstInference, aes(x=timeToBurst, color=biotype, fill=biotype)) +
      geom_histogram(aes(y=..density..), alpha = 0.25, breaks=tmp, position = 'identity', size=0.25) +
      scale_color_manual(values=c(col_nc, col_coding)) +
      scale_fill_manual(values=c(col_nc, col_coding), labels = c(n.lncRNA, n.coding)) +
      scale_x_continuous(breaks=seq(0, xLim, 12)) +
      xlab('Time to next burst (h)') +
      ggtitle(label) +
      geom_vline(xintercept = c(medCoding, medNC), linetype="dashed", color = c(col_coding, col_nc), size=0.25) +
      geom_vline(xintercept = 24, linetype="dashed", color = 'grey50', size=0.25) +
      theme_man +
      theme(
        #legend.position = c(0.5, 0.8),
        legend.position = 'top',
        legend.title = element_text(size=6),
        legend.text=element_text(size=6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        plot.title = element_text(size = 6)
      )
    return(gg)
  }  



### Plot CV2 to mean expression and highlight highly ranked genes (lncRNAs, ranked by CV2)
plot_topRanks <- function(geneStats, cvRankings, minMean=0.1, label = NA){
  
  topGenes <- geneStats %>% filter(GeneStableID %in% (cvRankings %>% filter(highRank))$GeneStableID)
  
  col_coding <- "#39B54A"
  col_nc <- "#21409A"
  
  gg <- 
    ggplot(geneStats, aes(x=mean, y=cv2, color=biotype, size=biotype, alpha=biotype)) +
    geom_point(shape=16) +
    annotation_logticks(sides = 'bl', size=0.25, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    scale_color_manual(values=c(col_nc, col_coding)) +
    scale_size_manual(values=c(0.75, 0.25)) +
    scale_alpha_manual(values=c(0.80, 0.5)) +
    ylab('CV2 (Allelic UMIs)') +
    xlab('Mean (allelic UMIs) ') +
    geom_point(data=topGenes, aes(x=mean, y = cv2), color='red', shape=1, inherit.aes = FALSE, size=1, stroke=0.25, alpha=0.95) +
    geom_vline(xintercept=minMean, linetype="dashed", color = "red", size=0.25) + 
    scale_x_log10() +
    scale_y_log10() +
    ggtitle(label) +
    theme_man + 
    theme(
      legend.position = c(0.8, 0.95),
      plot.title = element_text(size = 5)
    )
  return(gg)
}

### Plot bursting kinetics for sampled expession-matched genes
makePlots <- function(sampled_Genes, lncRNA_stats, nPerm, label){
  col_coding = "#39B54A"
  col_nc = "#21409A"
  
  # kon
  p = length(which(sampled_Genes$kon < lncRNA_stats$kon)) / nPerm
  ggKon <- 
    ggplot(sampled_Genes, aes(x=kon)) +
    geom_histogram(fill=col_coding, color='grey', size=0.25) +
    ggtitle(paste0(label, ', p < ', p)) +
    geom_vline(xintercept=lncRNA_stats$kon, linetype="dashed", color = col_nc, size=0.25) +
    theme_man +
    theme(
      plot.title = element_text(size = 6)
    )
  
  # Size
  p = length(which(sampled_Genes$size > lncRNA_stats$size)) / nPerm
  ggSize <- 
    ggplot(sampled_Genes, aes(x=size)) +
    geom_histogram(fill=col_coding, color='grey', size=0.25) +
    ggtitle(paste0(label, ', p < ', p)) +
    geom_vline(xintercept=lncRNA_stats$size, linetype="dashed", color = col_nc, size=0.25) +
    theme_man +
    theme(
      plot.title = element_text(size = 6)
    )
  
  # Mean
  ggMean <- 
    ggplot(sampled_Genes, aes(x=meanExp_median)) +
    geom_histogram(fill=col_coding, color='grey') +
    geom_vline(xintercept=lncRNA_stats$meanExp, linetype="dashed", color = col_nc, size=0.25) +
    xlab('Mean') +
    theme_man +
    theme(
      plot.title = element_text(size = 6)
    )
  
  # Dump plots
  dump.l <- list()
  dump.l[[1]] = ggKon
  dump.l[[2]] = ggSize
  dump.l[[3]] = ggMean
  names(dump.l) <- c('kon', 'size', 'mean')
  
  return(dump.l)
}



