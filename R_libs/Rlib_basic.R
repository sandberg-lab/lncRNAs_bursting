### Nice label for plot log axis
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

### Add biotype/GeneType
addBiotype <- function(df, biomart, addBiomart, lab){
  tmp <- colnames(df)
  biomart <- biomart[match(df$GeneStableID, biomart$GeneStableID), ]
  df <- df %>% mutate(tmp = addBiomart)
  colnames(df) <- c(tmp, lab)
  return(df)
}

### Filter genes
filter_genes <- function(counts, min_expr, min_cells){
  passed_genes <- row.names(counts)[which(rowSums(counts >= min_expr) >= min_cells)]
  return(passed_genes)
}

### Filter samples on counts
filter_samples_4counts <- function(counts, min.counts){
  sample.counts <- which(colSums(counts) >= min.counts)
  pass.samples = names(sample.counts)
  return(pass.samples)
}

### Filter samples on genes
filter_samples_4genes <- function(counts, min.counts, min.genes){
	passed_samples <- col.names(counts)[which(colSums(counts >= min.counts) >= min.genes)]
	return(passed_samples)
}

### Stats summary
basic_summary_stats <- function(exp, genes) {
  exp.sel <- exp[genes, ]
  gene.means = rowMeans(exp.sel)
  gene.vars = apply(exp.sel, 1, var)
  gene.cv2 = gene.vars / gene.means^2
  gene.stats = cbind(gene.means, gene.vars, gene.cv2)
  colnames(gene.stats) = c('mean', 'var', 'cv2')
  gene.stats <- gene.stats %>% as.data.frame() %>% mutate(GeneStableID = genes)
  rownames(gene.stats) = c()
  return(gene.stats)
}

### Get genes in proximity
add_genes_in_prox <- function(gene, biomart, dist_lim){
  biomart_gene <- biomart %>% filter(GeneStableID %in% gene)
    biomart_tmp <-
    biomart %>%
    filter(Chromosome %in% biomart_gene$Chromosome) %>%
    filter(GeneStart > (biomart_gene$GeneStart - dist_lim) & GeneStart < (biomart_gene$GeneStart + dist_lim)) %>%
    mutate(dist_to_tss = GeneStart - biomart_gene$GeneStart)
    return(biomart_tmp)
  }

### Split genes into autosomal, X, imprinter-autosmal
splitGenes_autosomeXimprint <- function(genes_v, biomart_df, imprint_df) {
   biomart_tmp <- biomart_df %>% filter(GeneStableID %in% genes_v)

   genes.auto <- biomart_tmp[which(biomart_tmp$Chromosome != 'X' & biomart_tmp$Chromosome != 'Y'), ]$GeneStableID
   genes.X <-  biomart_tmp[which(biomart_tmp$Chromosome == 'X'), ]$GeneStableID

   genes.imprint <-
     biomart_tmp %>%
     filter(GeneName %in% imprint_df$Gene) %>%
     select(GeneStableID)

   genes.auto <- setdiff(genes.auto,  genes.imprint$GeneStableID)
   genes.imprint <- setdiff(genes.imprint$GeneStableID, genes.X)

   genes.l <- list()
   genes.l[[1]] <-  genes.auto
   genes.l[[2]] <- genes.X
   genes.l[[3]] <- genes.imprint
   names(genes.l) <- c('autosomes', 'X', 'imprint')
   return(genes.l)
}
