# R code to annotate gene types ()

find_minDist2TSS <- function(geneSel, exp, geneAnnotations) {
  geneAnnotation_geneSel <- geneAnnotations %>% filter(GeneStableID %in% geneSel)

  geneAnnotations_chr <-
    geneAnnotations %>%
    filter(
      GeneStableID %in% rownames(exp) &
      Chromosome %in% geneAnnotation_geneSel$Chromosome &
      GeneStableID != geneAnnotation_geneSel$GeneStableID
    ) %>%
    mutate(
      Dist2TSS = abs(.$TSS - geneAnnotation_geneSel$TSS)
    ) %>%
    arrange(Dist2TSS)

  geneInProx <- geneAnnotations_chr[1, ] %>% select(Dist2TSS, GeneStableID)
  colnames(geneInProx) <- c('minDist2TSS', 'minDist2TSSGene')

  return(geneInProx)
}

find_Divergent <- function(geneSel, geneAnnotions, dist){
  geneAnnotions_geneSel <- geneAnnotions %>% filter(GeneStableID %in% geneSel)
  geneAnnotions_geneInProx <- geneAnnotions %>% filter(GeneStableID %in% geneAnnotions_geneSel$minDist2TSSGene)

  # tmp output
  tmp = FALSE

    # Divergent on positive strand
    if(geneAnnotions_geneSel$Strand > 0 & geneAnnotions_geneInProx$Strand < 0 & geneAnnotions_geneSel$minDist2TSS < dist){
      if((geneAnnotions_geneSel$TSS - geneAnnotions_geneInProx$TSS) > 0){
        tmp = TRUE
      }
    }

    # Divergent on Negative strand
    if(geneAnnotions_geneSel$Strand < 0 & geneAnnotions_geneInProx$Strand > 0 &  geneAnnotions_geneSel$minDist2TSS < dist){
      if((geneAnnotions_geneSel$TSS - geneAnnotions_geneInProx$TSS) < 0){
        tmp = TRUE
      }
    }
  return(tmp)
  }

  find_Convergent <- function(geneSel, geneAnnotions, dist){
    geneAnnotions_geneSel <- geneAnnotions %>% filter(GeneStableID %in% geneSel)
    geneAnnotions_geneInProx <- geneAnnotions %>% filter(GeneStableID %in% geneAnnotions_geneSel$minDist2TSSGene)

    # tmp output
    tmp = FALSE

      # Divergent on positive strand
      if(geneAnnotions_geneSel$Strand > 0 & geneAnnotions_geneInProx$Strand < 0 & geneAnnotions_geneSel$minDist2TSS < dist){
        if((geneAnnotions_geneSel$TSS - geneAnnotions_geneInProx$TSS) < 0){
          tmp = TRUE
        }
      }

      # Divergent on Negative strand
      if(geneAnnotions_geneSel$Strand < 0 & geneAnnotions_geneInProx$Strand > 0 &  geneAnnotions_geneSel$minDist2TSS < dist){
        if((geneAnnotions_geneSel$TSS - geneAnnotions_geneInProx$TSS) > 0){
          tmp = TRUE
        }
      }
    return(tmp)
    }
