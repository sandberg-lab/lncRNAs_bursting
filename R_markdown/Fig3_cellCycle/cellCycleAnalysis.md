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
  
  meta_path = file.path(prj_path, 'data/ss2_n533_meta.rds')
  biomart_path = file.path(prj_path, 'data/ss2_n533_biomart.rds')
  
  cellCycle_path <- file.path(prj_path, 'data/Whitfield.txt')

# Lib
  library('stats')
  library('ggplot2')
  library('dplyr')
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
  library('tidyr')
  library('ggpubr')
  library('ggrastr')
  library('scales')
  library('statmod')
  library('princurve') # To align cells through the cell cycle
  library('zoo')  # Running mean
```

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

``` r
  #library('Seurat')
  #library('Matrix') # Needed to read Seurat format
  
  source(file.path(prj_path, '/R_libs/Rlib_basic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_cellCycle.R'))

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
# Winzorize counts (before finding most variable genes)
  nWin = 1
  counts_filt_win <- t(apply(counts_filt, 1, winsorize, n.win = nWin))
  
# Use Searuat to get variability of genes
  # Make Seurat object and normalise data
  fib = Seurat::CreateSeuratObject(counts = counts_filt_win, project = "fib", min.cells = 0, min.features = 0)
```

    ## Warning: The following arguments are not used: row.names

``` r
  fib <- Seurat::NormalizeData(fib, normalization.method = "LogNormalize", scale.factor = 10000)
  fib <- Seurat::FindVariableFeatures(fib, selection.method = "vst", nfeatures = nrow(counts_filt_win))

  # Variability for cell cycle genes (ranked)
  Var_rank <-
    Seurat::VariableFeatures(fib) %>%
    as.data.frame() %>%
    filter(Seurat::VariableFeatures(fib) %in% cellCycle_filt$ensembl) %>%
    mutate(Rank=c(1:nrow(cellCycle_filt)))
  
  colnames(Var_rank) <- c('Gene', 'Rank')

  # Normalized data
  counts_filt_win_norm <- as.matrix(fib$RNA@data)
  

# Subset of cell cycle genes for downstream analysis
  n.topVar = 50
  topVar_rank <- Var_rank[1:n.topVar, ]

  
# Non-winsorize normalized data for Anova test
  fib_noWin = Seurat::CreateSeuratObject(counts = counts_filt, project = "fib_noWin", min.cells = 0, min.features = 0)
```

    ## Warning: The following arguments are not used: row.names

``` r
  fib_noWin <- Seurat::NormalizeData(fib_noWin, normalization.method = "LogNormalize", scale.factor = 10000)

  counts_filt_norm = as.matrix(fib_noWin$RNA@data)
```

``` r
# make df for ggplot
  df <-
    data.frame(
      'GeneStableID'=rownames(fib$RNA[[1]]),
      'mean'=fib$RNA[[1]]$vst.mean,
      'variance.standardized'=fib$RNA[[4]]$vst.variance.standardized
      ) %>%
    mutate(
      genes.sel = ifelse(GeneStableID %in% cellCycle_filt$ensembl, 'cellCycle', 'other')
    )

    df[match(topVar_rank$Gene, df$GeneStableID), ]$genes.sel <- 'top.n'
    df <- df[order(df$genes.sel), ]
  

# plot
  lab1 = paste0('n.genes: ', nrow(df))
  lab2 = paste0('n.var.cc: ', nrow(cellCycle_filt))
  lab3 = paste0('n.var.cc.top: ', length(which(df$genes.sel=='top.n')))
  
  ExtendedDataFig5A <-
    ggplot(df, aes(x=mean, y=variance.standardized, color=genes.sel, size=genes.sel, shape=genes.sel)) +
    geom_point() +
    scale_x_log10(labels=comma) +
    scale_color_manual(values=c('darkgreen', 'grey', 'red'), labels = c(lab1, lab2, lab3)) +
    scale_alpha_manual(values=c(1, 0.5, 1)) +
    scale_size_manual(values=c(1.25, 0.25, 1.25)) +
    scale_shape_manual(values=c(16, 16, 16)) +
    annotation_logticks(sides = "b") +
    ggtitle('Most variable cell cycle genes') +
    xlab('Mean') + 
    ylab('Standardized variance') +
    theme_man +
    theme(legend.position = c(0.2, 0.8),
          legend.title = element_blank(),
          legend.text=element_text(size=5)
          )

ExtendedDataFig5A
```

![](cellCycleAnalysis_files/figure-markdown_github/Plot%20most%20variable%20cell%20cycle%20genes-1.png)

``` r
# Scale data
  fib <- Seurat::ScaleData(fib, features = rownames(fib))
```

    ## Centering and scaling data matrix

``` r
  counts_filt_win_norm_scale <- fib$RNA@scale.data

# PCA using top var cell cycle genes
  counts_sel <- counts_filt_win_norm_scale[topVar_rank$Gene, ]

### PCA using selected genes
  pca = prcomp(t(counts_sel), scale=F, center=F) # Option to scale within Seurat

  df.pca <- 
    pca$x %>%
    as.data.frame() %>%
    select(PC1, PC2, PC3)

  ggPCA <-
    ggplot(df.pca, aes(x=PC1, y=PC2)) +
    geom_point() +
    theme_man

  ggPCA
```

![](cellCycleAnalysis_files/figure-markdown_github/PCA%20using%20most%20varible%20cell%20cycle%20genes-1.png)

``` r
Gas1 <- plot_pca(counts_filt_win_norm , df.pca, 'Gas1', biomart_filt)
Ccnd1 <- plot_pca(counts_filt_win_norm , df.pca, 'Ccnd1', biomart_filt)
Ccne2 <- plot_pca(counts_filt_win_norm , df.pca, 'Ccne2', biomart_filt)
Ccna2 <- plot_pca(counts_filt_win_norm , df.pca, 'Ccna2', biomart_filt)

ggarrange(Gas1, Ccnd1, Ccne2, Ccna2, ncol = 4, nrow = 1)
```

![](cellCycleAnalysis_files/figure-markdown_github/Plot%20PCA%20with%20cell%20cycle%20expression-1.png)

``` r
# Add princurve analysis to meta
  meta <-
    meta %>%
    mutate(
      lambda = principal_curve(as.matrix(df.pca[.$Sample, ]))$lambda
    ) %>%
    arrange(lambda)

# Select some genes to plot and add rolling mean
  genes_sel <- c('Ccne2', 'Ccna2', 'Gas1', 'Ccnd1')
  biomart_sel <- biomart_filt %>% filter(GeneName %in% genes_sel)

  counts_filt_win_norm_sel <- counts_filt_win_norm[biomart_sel$GeneStableID, meta$Sample]
  rownames(counts_filt_win_norm_sel) <- biomart_sel$GeneName
  
  df.cc <- rolling_mean(counts_filt_win_norm_sel,  meta$Sample, genes_sel, 15)

# Set borders for cell cycle (based on maunal inspection of rolling mean)
  cc.breaks <- c(75, 310, 430)
  cc.col <- c(col_sel[1], col_sel[2], col_sel[3], col_sel[4])

  # G0
  meta$cc.phase <- 'G0'
  meta$cc.phase.col <- cc.col[1]

  # G1
  meta[c(cc.breaks[1]:cc.breaks[2]), ]$cc.phase <- 'G1'
  meta[c(cc.breaks[1]:cc.breaks[2]), ]$cc.phase.col <- cc.col[2]

  # G1S
  meta[c(cc.breaks[2]:cc.breaks[3]), ]$cc.phase <- 'G1S'
  meta[c(cc.breaks[2]:cc.breaks[3]), ]$cc.phase.col <- cc.col[3]

  # G2M
  meta[c(cc.breaks[3]:nrow(meta)), ]$cc.phase <- 'G2M'
  meta[c(cc.breaks[3]:nrow(meta)), ]$cc.phase.col <- cc.col[4]

### Plot individual candidates (without cells)  
  gg1 <- plot_gene_cellCycle(df.cc, 'Gas1', cc.breaks, cc.col, cc.col[1])
  gg2 <- plot_gene_cellCycle(df.cc, 'Ccnd1', cc.breaks, cc.col, cc.col[3])
  gg3 <- plot_gene_cellCycle(df.cc, 'Ccne2', cc.breaks, cc.col, cc.col[2])
  gg4 <- plot_gene_cellCycle(df.cc, 'Ccna2', cc.breaks, cc.col, cc.col[4])
  
  # ExtendedDataFig 5a
  ggarrange(gg1, gg2, gg3, gg4, ncol = 2, nrow = 2)
```

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

![](cellCycleAnalysis_files/figure-markdown_github/Use%20princurve%20and%20PCA%20coordinates%20to%20align%20cells%20to%20cell%20cycle%20progression-1.png)

``` r
# Color PCA accodring to cell cycle phase
  df.pca <- 
    df.pca[meta$Sample, ] %>%
    mutate(
      cc.phase = meta$cc.phase
    )
  

  col.tmp <- c(col_sel[2], col_sel[3], col_sel[1], col_sel[4])

  ggplot(df.pca, aes(x=PC1, y=PC2, col=cc.phase)) +
    geom_point() +
    scale_color_manual(values=col.tmp) +
    theme_man
```

![](cellCycleAnalysis_files/figure-markdown_github/Plot%20PC1%20vs%20PC2%20colored%20by%20cell%20cycle%20phase-1.png)

``` r
# Make df with gene expression to plot
  df.tmp <- 
    data.frame(
      cc.phase = meta$cc.phase,
      Gas1 = counts_filt_win_norm[(biomart_filt %>% filter(GeneName %in% 'Gas1'))$GeneStableID ,meta$Sample],
      Ccnd1 = counts_filt_win_norm[(biomart_filt %>% filter(GeneName %in% 'Ccnd1'))$GeneStableID ,meta$Sample],
      Ccne2 = counts_filt_win_norm[(biomart_filt %>% filter(GeneName %in% 'Ccne2'))$GeneStableID ,meta$Sample],
      Ccna2 = counts_filt_win_norm[(biomart_filt %>% filter(GeneName %in% 'Ccna2'))$GeneStableID ,meta$Sample]
    )
  
plot_tmp <- function(df, gene, colors){
  return(   
    ggplot(df.tmp, aes(x=cc.phase, y=get(gene), colour=cc.phase)) + 
    geom_boxplot() + 
    scale_color_manual(values=col.tmp, labels = c('G0', 'G1', 'G1S', 'G2M')) +
    xlab('Cell cycle phase') +
    ylab('Normalized expression') +
    ggtitle(gene) +
    theme_man +
    theme(
      legend.position = 'none')
  )
}
  
  ggGas1 <- plot_tmp(df.tmp, gene='Gas1', colors=c(col_sel[1], col_sel[3], col_sel[2], col_sel[4]))
  ggCcnd1 <- plot_tmp(df.tmp, gene='Ccnd1', colors=c(col_sel[1], col_sel[3], col_sel[2], col_sel[4]))
  ggCcne2 <- plot_tmp(df.tmp, gene='Ccne2', colors=c(col_sel[1], col_sel[3], col_sel[2], col_sel[4]))
  ggCcna2 <- plot_tmp(df.tmp, gene='Ccna2', colors=c(col_sel[1], col_sel[3], col_sel[2], col_sel[4]))
  
  ggarrange(ggGas1, ggCcnd1, ggCcne2, ggCcna2, ncol = 2, nrow = 2)
```

![](cellCycleAnalysis_files/figure-markdown_github/Cell%20cycle%20genes%20colored%20by%20cell%20cycle%20phase-1.png)

``` r
# Level of significance
  pAnova.sig = 0.01

# Expressed lncRNAs
  lncRNA_filt <- 
    biomart_filt %>%
    filter(lncRNA %in% TRUE) %>%
    filter(GeneStableID %in% rownames(counts_filt_norm)) %>% 
    select(GeneStableID, GeneName)
  
    genes_anova <-
      anova_test(
        expression = counts_filt_norm[ ,meta$Sample],
        genes = lncRNA_filt$GeneStableID,
        cells = meta$Sample,
        group = meta$cc.phase,
        adj = 'BH'
        ) %>%
      mutate(
        GeneName = lncRNA_filt$GeneName,
        sig = ifelse(.$p.adj < pAnova.sig, TRUE, FALSE)
      ) %>%
      arrange(p.adj)

#head(genes_anova, 15)
```

``` r
# Fnc to add mean expression for a specific cell cycle phase
  mean_cellCycle_phase <- function(counts, meta, cc.sel) {
    return(
      counts[ ,(meta %>% filter(cc.phase %in% cc.sel) %>% select(Sample))[,1]] %>%
        apply(1, mean)
    )    
  }

 
# Set some threshold
  thr.cutoff = 1
  thr.fold = 1.5
  thr.p = 0.01
  
# Fold change over in cell cycle phase
  genes_anova <-
    genes_anova %>%
    mutate(
      mean.G0 = mean_cellCycle_phase(counts_filt_norm[genes_anova$GeneID, ], meta, 'G0'),
      mean.G1 = mean_cellCycle_phase(counts_filt_norm[genes_anova$GeneID, ], meta, 'G1'),
      mean.G1S = mean_cellCycle_phase(counts_filt_norm[genes_anova$GeneID, ], meta, 'G1S'),
      mean.G2M = mean_cellCycle_phase(counts_filt_norm[genes_anova$GeneID, ], meta, 'G2M') 
    ) %>%
    mutate(
      fold.G0 = mean.G0 / apply(counts_filt_norm[genes_anova$GeneID, (meta %>% filter(cc.phase %notin% 'G0') %>% select(Sample))$Sample], 1, mean),
      fold.G1 = mean.G1 / apply(counts_filt_norm[genes_anova$GeneID, (meta %>% filter(cc.phase %notin% 'G1') %>% select(Sample))$Sample], 1, mean),
      fold.G1S = mean.G1S / apply(counts_filt_norm[genes_anova$GeneID, (meta %>% filter(cc.phase %notin% 'G1S') %>% select(Sample))$Sample], 1, mean),
      fold.G2M = mean.G2M / apply(counts_filt_norm[genes_anova$GeneID, (meta %>% filter(cc.phase %notin% 'G2M') %>% select(Sample))$Sample], 1, mean)
      )

# Plot and highlight candidates (Fig 3b)
  # A function to make plots: fold changes (X-axis) vs p-values (y-axis)
  fnc_toPlot <- function(anovaStats, genesSelected, ccPhase='G0', foldSig = 1.5, pSig = 0.01){
    
    df.tmp <-
      anovaStats %>%
      select(GeneName, p.adj, starts_with('fold')) %>%
      select(GeneName, p.adj, ends_with(ccPhase))
    colnames(df.tmp) <- c('GeneName', 'p.adj', 'fold')
    
    df.tmp <-
      df.tmp %>%
      mutate(
        sig = ifelse(p.adj < pSig & fold > foldSig, TRUE, FALSE),
        genesSel = ifelse(GeneName %in% genesSelected, TRUE, FALSE)
      ) %>%
      filter(
        fold > 1
      )
    
  gg <-
    ggplot(df.tmp, aes(x=log2(fold), y=-log10(p.adj), colour=genesSel, size=genesSel)) +
    geom_point(shape=16) + 
    scale_size_manual(values = c(0.5, 1)) +
    scale_colour_manual(values = c('lightgrey', 'red')) +
    geom_vline(xintercept=log2(thr.fold), linetype="dashed", size=0.25, color='red') +
    geom_hline(yintercept=-log10(thr.p), linetype="dashed", size=0.25, color='red') +
    xlab('Fold change (log2)') +
    ylab('p-value (BH adjusted)')+
    ggtitle(paste0(ccPhase, ', nSig: ', table(df.tmp$sig)[2])) +
    geom_text(
      data = df.tmp %>% filter(GeneName %in% genesSelected),
      aes(x=log2(fold), y=-log10(p.adj), label=GeneName),
      size=4, color='red', vjust=-0.5, inherit.aes = FALSE) +
    theme_man +
    theme(legend.position = 'none')
  return(gg)
  }


# Dump plots
  ggG0 <-
    fnc_toPlot(
      genes_anova,
      c('A730056A06Rik', 'Gm12963', 'Mir22hg'),
      ccPhase='G0')
  
  ggG1S <-
    fnc_toPlot(
      genes_anova,
      c('Wincr1', '1600019K03Rik'),
      ccPhase='G1S')

  ggG2M <-
  fnc_toPlot(
    genes_anova,
    c('Lockd', '2010110K18Rik'),
    ccPhase='G2M')
  
ggarrange(ggG0, ggG1S, ggG2M, ncol = 3, nrow = 1)
```

![](cellCycleAnalysis_files/figure-markdown_github/Plot%20anova%20stats%20vs%20fold%20changes-1.png)

``` r
# Select some genes to plot and add rolling mean
  genes_sel <- c('A730056A06Rik', 'Gm12963', 'Mir22hg', 'Wincr1', '1600019K03Rik', 'Lockd', '2010110K18Rik')
  biomart_sel <- biomart_filt %>% filter(GeneName %in% genes_sel)

  counts_filt_win_norm_sel <- counts_filt_win_norm[biomart_sel$GeneStableID, meta$Sample]
  rownames(counts_filt_win_norm_sel) <- biomart_sel$GeneName
  
  df.cc <- rolling_mean(counts_filt_win_norm_sel,  meta$Sample, genes_sel, 15)

  ### Plot individual candidates (without cells)  
  gg1 <- plot_gene_cellCycle(df.cc, genes_sel[1], cc.breaks, cc.col, cc.col[3])
  gg2 <- plot_gene_cellCycle(df.cc, genes_sel[2], cc.breaks, cc.col, cc.col[3])
  gg3 <- plot_gene_cellCycle(df.cc, genes_sel[3], cc.breaks, cc.col, cc.col[3])
  gg4 <- plot_gene_cellCycle(df.cc, genes_sel[4], cc.breaks, cc.col, cc.col[3])
  gg5 <- plot_gene_cellCycle(df.cc, genes_sel[5], cc.breaks, cc.col, cc.col[3])
  gg6 <- plot_gene_cellCycle(df.cc, genes_sel[6], cc.breaks, cc.col, cc.col[3])
  gg7 <- plot_gene_cellCycle(df.cc, genes_sel[7], cc.breaks, cc.col, cc.col[3])
    
ggarrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, ncol=3, nrow=3)
```

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

    ## Warning: Removed 14 row(s) containing missing values (geom_path).

![](cellCycleAnalysis_files/figure-markdown_github/Expression%20of%20lncRNA%20candidates%20during%20cell%20cycle%20progression%20(rolling%20mean)-1.png)
