### Binomial test for allelic imbalance
  binomTest_AllelicImbalance <- function(counts_allele1, counts_allele2){
    sum_allele1 <- rowSums(counts_allele1)
    sumTot <- rowSums(counts_allele1 + counts_allele2)
  
    test = vector()
    for(i in 1:length(sum_allele1)){
      test[i] <- binom.test(sum_allele1[i], sumTot[i], alternative=c('two.sided'))$p.value
    }
  
    tmp.l <- list()
    tmp.l[[1]] <- test
    tmp.l[[2]] <- p.adjust(test, method = 'BH', n = length(test))
    names(tmp.l) <- c('pVal', 'pVal_adj')

    return(tmp.l)
  }  

  
### Get genes in prox with stats (allelicImbalance)
  get_genesInProx_wStats <- function(geneSel, biomart, dist){
    
    # Genes on same chromosome
    biomartSel <- biomart %>% filter(GeneStableID %in% geneSel)
    
    # Genes in proximity
    biomart.tmp <-
      biomart %>%
      filter(
        Chromosome %in% biomartSel$Chromosome
      ) %>%
      mutate(
        distToTSS = TSS - biomartSel$TSS
      ) %>%
      filter(
        coding & abs(.$distToTSS) < dist
      )
    
    # Stats for genes in proximity (if any genes)
    if(nrow(biomart.tmp) > 0){
      biomart.tmp <- 
        biomart.tmp %>%
        mutate(
          lncRNASel = biomartSel$GeneStableID,
          allelicImbalance_lncRNA = biomartSel$allelicImbalance,
          allelicImbalanceTest =
            abs(.$allelicImbalance) + abs(biomartSel$allelicImbalance) -
            abs(abs(.$allelicImbalance) - abs(biomartSel$allelicImbalance))
        ) %>%
        mutate(
          allelicDirection = ifelse(.$allelicImbalance_lncRNA > 0 & .$allelicImbalance > 0 | .$allelicImbalance_lncRNA < 0 & .$allelicImbalance < 0, 1, -1)
        ) %>%
        filter(
          GeneStableID %notin% biomartSel$GeneStableID
        )
      
    }
    return(biomart.tmp)
  }
  

### Get genes in prox
  get_genesInProx <- function(geneSel, biomart, dist){
  
  # Genes on same chromosome
    biomartSel <- biomart %>% filter(GeneStableID %in% geneSel | GeneName %in% geneSel)
  
  # Genes in proximity
    biomart.tmp <-
      biomart %>%
      filter(
        Chromosome %in% biomartSel$Chromosome
      ) %>%
      mutate(
        distToTSS = TSS - biomartSel$TSS
      ) %>%
      filter(
        coding & abs(.$distToTSS) < dist
      )
    return(biomart.tmp)
  }



### Random sampling/translocation of genes in proximity
  randomSampling_allelicImbalance <- function(gene, biomart, dist, randomGenes){
  # Get gene of interest
    biomartGene <- biomart %>% filter(GeneStableID %in% gene)
  
  # Exlclude genes in proximity for downstream permutation test
    genesToExclude <-
      biomart %>%
      filter(Chromosome %in% biomartGene$Chromosome) %>%
      filter(
        TSS < biomartGene$TSS + 2*dist &
          TSS > biomartGene$TSS - 2*dist
      )
  
    biomartPass <-
      biomart %>%
      filter(
        coding & GeneStableID %notin% genesToExclude$GeneStableID
      )
  
  # Random sampling
    FNC <- function(randomGenes, biomartGene, biomartPass, dist){
      # Sample a gene
      biomartRandom <- 
        biomartPass %>%
        filter(GeneStableID %in% randomGenes)
    
    # Get genes in proximity
      genesInProx <- 
        biomartPass %>%
        filter(Chromosome %in% biomartRandom$Chromosome) %>%
        filter(
          TSS < biomartRandom$TSS + dist &
            TSS > biomartRandom$TSS - dist
        ) %>%
        filter(
          GeneStableID %notin% biomartRandom$GeneStableID
        )
    
    # Get stats for allelic imbalance
      if(nrow(genesInProx)>0){
        test.df <- 
          data.frame(
            'distToTSS' = genesInProx$TSS - biomartRandom$TSS,
            'allelicImbalanceTest' =
              (abs(genesInProx$allelicImbalance) + abs(biomartGene$allelicImbalance)) -
              abs(abs(genesInProx$allelicImbalance) - abs(biomartGene$allelicImbalance)),
            'permutation' = i,
            'randomGene' = biomartRandom$GeneStableID,
            'selectedGene' = biomartGene$GeneStableID
          )
      }
      if(nrow(genesInProx)==0){
        test.df <- 
          data.frame(
            'distToTSS' = NA,
            'allelicImbalanceTest' = NA,
            'permutation' = i,
            'randomGene' = NA,
            'selectedGene' = NA
          )
      }
      return(test.df)
    }
    test.l <- parallel::mclapply(randomGenes, FNC, biomartGene, biomartPass, dist, mc.cores = 10)  
    test.df <- do.call(rbind.data.frame, test.l)
    return(test.df)
  }
