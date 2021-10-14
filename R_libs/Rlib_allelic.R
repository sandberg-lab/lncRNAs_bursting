### Filter genes (for allelic counts)
filter_genes_4allelicCounts <- function(counts1, counts2, min_expr, min_cells){
  counts1[is.na(counts1)] <- 0
  counts2[is.na(counts1)] <- 0
  counts1_genes <- row.names(counts1)[which(rowSums(counts1 >= min_expr) >= min_cells)]
  counts2_genes <- row.names(counts2)[which(rowSums(counts2 >= min_expr) >= min_cells)]

  passed_genes <- union(counts1_genes, counts2_genes)
  return(passed_genes)
}


### Number of genes with allelic counts (per cell)
genes_per_cell_wAllelicCounts <- function(counts1, counts2, min_counts) {
  counts1_genes <- apply(counts1, 2, function(x) {which(x >= min_counts)})
  counts2_genes <- apply(counts2, 2, function(x) {which(x >= min_counts)})

  tmp = vector()
  for(i in 1:ncol(counts1)) {
    tmp[i] <- length(union(counts1_genes[[i]], counts2_genes[[i]]))
  }
  return(tmp)
}

### Allelic distribution of read counts
allelic_balance_perCell <- function(counts1, counts2, genes) {
	counts1_tmp <- counts1[genes, ]
	counts2_tmp <- counts2[genes, ]
	counts_sum <- colSums(counts1_tmp) + colSums(counts2_tmp)
	dScore <- colSums(counts1_tmp) / counts_sum - 0.5
	return(dScore)
	}

dscore_gene <- function(counts1, counts2, genes) {
	counts1_tmp <- counts1[genes, ]
	counts2_tmp <- counts2[genes, ]
	counts_sum <- rowSums(counts1_tmp) + rowSums(counts2_tmp)
	dScore <- rowSums(counts1_tmp) / counts_sum - 0.5
	return(dScore)
	}
