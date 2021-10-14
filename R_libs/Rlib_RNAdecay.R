# Calculate k (decay rate) from half-life
  getDecayRate <- function(tHalfLife){
    return(-log(0.5/1) * 1/tHalfLife)
  }

# Calculates y: y = a exp(-k * x) when k, x and a are known
  getExpectedExpression <- function(k, t, expression.t0){
    #y = a*exp(-k*x)
    return(expression.t0 * exp(-k*t))
  }

# Normalizes  data to some genes with known k / half-life
  normData_bySubSetGene <- function(expr, meta, biomart, sampleID, selectedGenes_wExpExpr){
    #selectedGenes_wExpExpr (with Expexted Expression): output from getExpectedExpression

  # Extract data for input parameters
    expr.sel <- expr[ ,(meta %>% filter(ID %in% sampleID))$XC]
    meta.sel <- meta %>% filter(ID %in% sampleID)
    biomart.sel <- biomart %>% filter(GeneName %in% selectedGenes_wExpExpr$GeneName) %>% select(GeneStableID, GeneName)

    selectedGenes_wExpExpr_long <- selectedGenes_wExpExpr %>% gather('timePoint', 'frcOf_t0', 5:10)

  # Calculate normalization factor for each time point
    normFactor = vector()
    expr.sel.norm <- expr.sel

    sample.t0 <- (meta.sel %>% filter(treatment %in% 't0'))$XC
    expr.t0 = expr.sel[biomart.sel$GeneStableID, sample.t0]

    for(i in c('t1', 't2', 't3', 't4', 't5')){
      sample.i <- (meta.sel %>% filter(treatment %in% i))$XC
      selectedGenes_wExpExpr.i <- selectedGenes_wExpExpr_long %>% filter(timePoint %in% i)
      expr.t0 <- expr.t0[selectedGenes_wExpExpr.i$GeneStableID]

      if(i == 't5' | i == 't4'){
        selectedGenes_wExpExpr.i <- selectedGenes_wExpExpr.i %>% filter(HalfLife_h > 2)
        expr.t0 <- expr.t0[selectedGenes_wExpExpr.i$GeneStableID]
      }

    tmp <-
      data.frame(
        GeneName = selectedGenes_wExpExpr.i$GeneName,
        expr.t0 = expr.t0,
        expr.exp = expr.t0 * selectedGenes_wExpExpr.i$frcOf_t0,
        expr.obs = expr.sel[names(expr.t0), sample.i]
      ) %>%
      mutate(
        normFactor =  .$expr.obs / .$expr.exp
      )

    normFactor[i] = median(tmp$expr.obs / tmp$expr.exp)
    expr.sel.norm[ ,sample.i] <- expr.sel.norm[ ,sample.i] / median(tmp$normFactor)
    }

  return(expr.sel.norm)
  }

# Calculate decay parameters
  getDecayParameters <- function(gene, rpkm.list, meta){
    df <-
      data.frame(
        rep1 = rpkm.list[[1]][gene, ],
        rep2 = rpkm.list[[2]][gene, ],
        rep3 = rpkm.list[[3]][gene, ],
        rep4 = rpkm.list[[4]][gene, ],
        time = unique(meta$t_min/60)
      )
    df.long <- df %>% gather(tTreatment, expr, 1:4)
    model <- drm(expr ~ time, fct = DRC.expoDecay(), data = df.long)
    #plot(model, log="")
    a = model$coefficients[1]
    k = model$coefficients[2]
    x = -log(0.5) / k
    return(data.frame(GeneStableID = gene, a = a , k = k, HalfLife_h = x))
  }

plotModelToRealData <- function(geneSel, expr, meta, biomart, decayParameters){
  # Select gene
  geneSel_ensemble <- biomart %>% filter(GeneName %in% geneSel)
  exprSel <- expr[geneSel_ensemble$GeneStableID, meta$XC]
  decayParametersSel <- decayParameters %>% filter(GeneName %in% geneSel)

  # Real data to plot
  t = c(0:960)/60

  df.obs <- rbind(
    data.frame(
      expr = exprSel,
      t = meta$t_min/60,
      sample = meta$ID
    ),
    data.frame(
      expr = decayParametersSel$a * exp(-decayParametersSel$k * t),
      t = t,
      sample = 'Fit'
    )
  )

  col_sel <-  c("#39B54A", "#27AAE1", "#21409A", "#EF4136", 'black')
  # Plot
  gg <-
    ggplot(df.obs, aes(x=t, y=expr, color=sample, size=sample)) +
    geom_point(shape = 16, alpha=0.95) +
    scale_color_manual(values=col_sel) +
    scale_size_manual(values=c(1.25, 1.25, 1.25, 1.25, 0.1)) +
    scale_x_continuous(limits=c(0, 16), breaks=c(0, 2, 4, 6, 8, 10, 12, 14, 16)) +
    xlab('Time (hrs)') +
    ylab('RPKM (norm)') +
    geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5) +
    geom_vline(xintercept=decayParametersSel$HalfLife_h, linetype="dashed", color = "red", size=0.5) +
    ggtitle(geneSel) +
    theme_man +
    theme(legend.position = c(0.8, 0.8))

    return(gg)
}
