# Match mean expression
  match_mean <- function(gene, geneStats, subset_genes, n_matched) {
    # gene: selected gene (to be matched)
    # geneStats: output from basic_summary_stats (mean, var, cv2, GeneStableID)
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
  
  # Dump gene(s) with most equal mean expression
  return(geneStats_sel$GeneStableID[1:n_matched])
}


# Rank cv2 (of i.e. lncRNAs) to genes of most similar mean expression (i.e. 100 mRNAs)
  rank_var <- function(gene, geneStats, exprMatchedGenes){
    # gene: selected gene
    # geneStats: output from basic_summary_stats (mean, var, cv2, GeneStableID)
    # exprMatchedGenes: output from match_mean (list with matched genes)
    
    tmp <-
      length(
        which(
          (geneStats %>% filter(GeneStableID %in% gene))$cv2 >
          (geneStats %>% filter(GeneStableID %in% exprMatchedGenes[[gene]]))$cv2
        )
      )
    
    df <- data.frame(
      GeneStableID = gene,
      rank = tmp/length(exprMatchedGenes[[gene]])
    )
  }

  
# Rank variability of randomly sampled genes (as abone, i.e. random mRNAs to other mRNAs)
  rank_var_perm <- function(i, geneStats, exprMatchedGenes, nGenes_sample) {
    genes.sampled <- sample(names(exprMatchedGenes), nGenes_sample)
    exprMatchedGenes.sampled <- exprMatchedGenes[genes.sampled]
    
    df <- 
      do.call(
        rbind.data.frame,
        lapply(genes.sampled, rank_var, geneStats=geneStats, exprMatchedGenes=exprMatchedGenes.sampled)
      ) %>%
      mutate(i = factor(i))
    return(df)
  }
  
  # Rank variability of randomly sampled genes (as above, i.e. random mRNAs to other mRNAs)
  sample_matchedMean <- function(nPerm, matched.df, expressionStats, wReplacement = TRUE){
    # nPerm: Number of permutations
    # matched.df: matched genes (lncRNA), output from match_mean in long df format
    # expressionStats: output from basic_summary_stats (mean, var, cv2, GeneStableID)
    `%notin%` <- Negate(`%in%`)
    
    ### Allow for the same gene to be selected/matched several times
    if(wReplacement){
      matched_genes_sampled <- 
        (matched.df %>%
           group_by(geneToMatch) %>%
           sample_n(1))$geneMatched
      
      stats_sampled <- expressionStats[match(matched_genes_sampled, expressionStats$GeneStableID), ]
      
      df.out <- 
        data.frame(
          cv2 = median(stats_sampled$cv2),
          meanExp_median = median(stats_sampled$mean),
          meanExp_mean = mean(stats_sampled$mean)
        )
    }

    ### Do not allow for the same gene to be selected/matched several times
    # NOTE: This is slow!
    if(!wReplacement){
      # Note; this may fail if no unique genes can be matched
      df.out = data.frame()
      
      # Randomise order of genes to match
      genes_sel <- unique(matched.df$geneToMatch)
      genes_sel <- genes_sel[sample(length(genes_sel))]
      
      # Sampled genes
      genes_sampled = vector()
      
      # For each gene; exclude previously matched genes and draw one of the remaining
      # NOTE: high frequency of non-matched genes
      for(i in 1:length(genes_sel)){
          
        # Genes to match (that have not been matched to other genes already)
          tmp <-
            matched.df %>%
            filter(geneToMatch %in% genes_sel[i]) %>%
            filter(geneMatched %notin% genes_sampled)
        
        # Add gene to vector if sampling was successful
          if(nrow(tmp)>0){
            genes_sampled[i] <- (tmp %>% sample_n(1))$geneMatched
            }
      
          }
      
      # Add to df.out If all genes have been matched (or dump empty df if not all genes were matched)
      if(length(which(is.na(genes_sampled))) == 0){
        expressionStats_sampled <- expressionStats %>% filter(GeneStableID %in% genes_sampled)

        df.out <- 
          data.frame(
            cv2 = median(expressionStats_sampled$cv2),
            meanExp_median = median(expressionStats_sampled$mean),
            meanExp_mean = mean(expressionStats_sampled$mean)
          )
      }
    }
    
    return(df.out)
  }
  
  
  
    