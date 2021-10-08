``` r
# Load data, packages and set output
  options(stringsAsFactors = FALSE)
  `%notin%` <- Negate(`%in%`)

# Basic folders
  home.dir <- '/home/perj/'

# Project path
  prj_name <- 'prj_lncRNAs'
  prj_path = file.path(home.dir, 'projects', prj_name)
  analysis = '1_make_biomart_geneAnnotation'

# Data
  data_path = file.path(prj_path, 'data')
  counts_path = file.path(data_path, 'ss2_n751_readCounts.rds')
    
# Resources
  biomart_path <- file.path(home.dir, 'resources/Biomart/Mouse/GRCm38.p6.csv')

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
  source(file.path(prj_path, '/R_libs/Rlib_basic.R'))
  source(file.path(prj_path, '/R_libs/Rlib_genetypes.R'))
```

``` r
# QC genes (note that data input is already QC)
  min.counts = 5
  min.cells = 2

# Filter counts and biomart  
  counts_qc <- counts[filter_genes(counts, min.counts, min.cells), ]
  
  biomart_qc <- 
    biomart %>%
    filter(GeneStableID %in% rownames(counts_qc)) %>%
    distinct(GeneStableID, .keep_all = TRUE) %>%
    select(GeneStableID, GeneName, GeneType, Strand, Chromosome, TSS, GeneStart, GeneEnd, TranscriptLength)
```

``` r
# Define transcript types
  dist_divergent = 500
  dist_convergent = 2000
  dist_independent = 4000
  dist_intergenic = 4000

# Add Transcript type (as logical annotations)
  biomart_qc <-
    biomart_qc %>%
    mutate(
      coding = ifelse(.$GeneType == 'protein_coding', TRUE, FALSE),
      lncRNA =  ifelse(.$GeneType == 'lncRNA', TRUE, FALSE)
    )
  
# Identify distance to gene in proximity (limited to QC expressed genes)
  genesInProx.l <- parallel::mclapply(biomart_qc$GeneStableID, find_minDist2TSS, counts, biomart_qc, mc.cores = 20)
  genesInProx.df <- as.data.frame(do.call(rbind, genesInProx.l))

  biomart_qc <- cbind(biomart_qc, genesInProx.df)

# Annotate transcript type
  
  # Independent
  biomart_qc <- 
    biomart_qc %>%
    mutate(
      independent = ifelse(minDist2TSS > dist_independent, TRUE, FALSE)
    )
  
  # Divergent
  biomart_qc <- 
    biomart_qc %>%
    mutate(
      divergent =
        as.vector(
          do.call(rbind,
                  parallel::mclapply(biomart_qc$GeneStableID, find_Divergent, biomart_qc, dist_divergent, mc.cores = 20)))
      )
  
  # Convergent  
  biomart_qc <- 
    biomart_qc %>%
    mutate(
      convergent = 
        as.vector(
          do.call(rbind,
                  parallel::mclapply(biomart_qc$GeneStableID, find_Convergent, biomart_qc, dist_convergent, mc.cores = 20)))
    )
  
# Intergenic genes
  GR <-
    GRanges(
    seqnames = biomart_qc$Chromosome,
    ranges = IRanges(start = biomart_qc$GeneStart,
                     end = biomart_qc$GeneEnd),
    strand = biomart_qc$GeneStrand
  )
  names(GR) <- biomart_qc$GeneName

  ol1 <- findOverlaps(GR, select='first', maxgap=dist_intergenic)
  ol2 <- findOverlaps(GR, select='last', maxgap=dist_intergenic)

  tmp <- data.frame(GeneName = biomart_qc$GeneName, ol1, ol2)
  biomart_qc$GeneOverlapEnsemble <- biomart_qc[tmp$ol2, ]$GeneStableID
  biomart_qc$GeneOverlapRefSeq <- biomart_qc[tmp$ol2, ]$GeneName

  biomart_qc$intergenic <-
    ifelse(
      biomart_qc$GeneStableID == biomart_qc$GeneOverlapEnsemble &
        biomart_qc$minDist2TSS > dist_intergenic, TRUE, FALSE)
```

``` r
dir.create('./data_dump', recursive=TRUE)
```

    ## Warning in dir.create("./data_dump", recursive = TRUE): './data_dump' already
    ## exists

``` r
  saveRDS(biomart_qc, file='./data_dump/ss2_fibs_biomart_wGeneAnnotations.rds')
```
