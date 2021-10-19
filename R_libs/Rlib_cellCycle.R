### Winsorize
  # use: t(apply(ed, 1, winsorize, 2))
  winsorize <- function(x, n.win = 1){
    n.vals = length(x)
    fraction = n.win / n.vals
    if(length(fraction) != 1 || fraction < 0 || fraction > 0.5){
      stop("bad value for 'fraction'")
    }
  
    win.sorted.ind = order(x)
    x[win.sorted.ind[1:n.win]] = x[win.sorted.ind[n.win + 1]]
    x[win.sorted.ind[(n.vals - n.win + 1):n.vals]] = x[win.sorted.ind[n.vals - n.win]]
  
    ##fraction-based
    if(0){
      lim = quantile(x, probs = c(fraction, 1 - fraction))
    
      ##set extreme values
      x[ x < lim[1] ] = lim[1]
      x[ x > lim[2] ] = lim[2]
    }
    return(x)
  }

# Plot PCA 
  plot_pca <- function(exp, df, gene, biomart) {
    biomart_sel <- biomart %>% filter(GeneName %in% gene | GeneStableID %in% gene)
  
    gene_sel <- exp[biomart_sel$GeneStableID, ]
    gg <-
      ggplot(df, aes(x=PC1, y=PC2, colour=gene_sel)) + 
      geom_point(shape=16, size=1) + 
      scale_colour_gradient(low="yellow", high="red", guide = guide_colourbar(title='')) +
      xlab('PC1') +
      ylab('PC2') +
      ggtitle(biomart_sel$GeneName) +
      theme_man +
      theme(
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.width = unit(0.1,"cm"))
    return(gg)
  }

  
### Anova test with adjustement
  anova_test <- function(expression, genes, cells, group, adj) {
    stats <- vector()
    
    for(i in 1:length(genes)) {
      temp.df <- data.frame('exp'=expression[genes[i], cells], group)
      res.aov <- aov(exp ~ group, data=temp.df)
      stats[i] <- summary(res.aov)[[1]][["Pr(>F)"]][1]
    }
    
    stats.adj <- p.adjust(stats, method=adj)
    
    df <- data.frame('GeneID'=genes, p=stats, p.adj=stats.adj)
    
    return(df)
  }
  
  
### Add rolling mean to df
  rolling_mean <- function(counts, cellOrder, genes, k){
    # k: size of window
    # k must be even
    
    countsOrdered <- counts[ ,cellOrder]
    tmp.x <- 1:length(cellOrder)
  
    df <-
      data.frame(
        'CellID' = cellOrder,
        'cell.order' = tmp.x
        )
  
    for(i in 1:length(genes)) {
      tmp.lab = colnames(df)
      df$tmp <- countsOrdered[genes[i], ]
      df$tmp.roll <- zoo::rollmean(countsOrdered[genes[i], ], k=k, na.pad=TRUE)
      tmp.lab <- c(tmp.lab, c(genes[i], paste0(genes[i], '.roll')))
      colnames(df) <- tmp.lab
    }
    return(df)
  }


# Plot gene expression (rolling mean) throughout cell cycle progression
  plot_gene_cellCycle <- function(df, geneSel, cellCycleBorders, cellCycleBorders_col, geneSel_col) {
    # Use input (df) from function rolling_mean()
    tmp <- which(colnames(df) == geneSel)
    
    df.tmp <-
      cbind(
        df$cell.order,
        df[ ,tmp:(tmp+1)]
      )
    colnames(df.tmp) <- c('cellOrder', 'gene', 'roll')
    
    y.tmp <- max(df.tmp$roll[!is.na(df.tmp$roll)])*1.05
    
    gg <-
      ggplot(df.tmp) +
      geom_line(data=df.tmp, aes(cellOrder, roll), colour=geneSel_col) + 
      geom_vline(xintercept=cellCycleBorders, linetype="dashed", color = "grey", size=0.25) +
      geom_segment(aes(x = 0, y = y.tmp, xend = cellCycleBorders[1], yend = y.tmp), colour=cellCycleBorders_col[1]) +  
      geom_segment(aes(x = cellCycleBorders[1], y = y.tmp, xend = cellCycleBorders[2], yend = y.tmp), colour=cellCycleBorders_col[3]) +  
      geom_segment(aes(x = cellCycleBorders[2], y = y.tmp, xend = cellCycleBorders[3], yend = y.tmp), colour=cellCycleBorders_col[2]) +  
      geom_segment(aes(x = cellCycleBorders[3] , y = y.tmp, xend = nrow(df.tmp), yend = y.tmp), colour=cellCycleBorders_col[4]) +  
      xlab('Cells (sorted by princurve lambda)') +
      ylab('Normalized exp (Seurat)') +
      theme_man +
      ggtitle(paste0(geneSel)) +
      theme(plot.title=element_text(face='italic', size=8, hjust=0.5))
    
    return(gg)
    
    }
    
    
    
    
  
    
    
    
    


