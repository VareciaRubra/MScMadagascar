DriftTestInterpretation <- function (tree, mean.list, cov.matrix.list, W.fixed = NULL, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(mean.list))) 
    stop("All tip labels must be in names(mean.list).")
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  
    if(is.null(W.fixed)) W.fixed <- cov.matrices[[length(cov.matrices)]] else W.fixed <-W.fixed
  
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  if (!any(node.mask)) 
    stop("For the regression test:/nAt least one node must have more than 4 descendents in mean.list")
  
  ind.cont<- ldply(mean.list, function(x) x) 
  rownames(ind.cont) <- ind.cont[,1]
  ind.cont <- ind.cont[,-1]
  ind.cont <- apply(ind.cont, 2, FUN = function (x) ape::pic(x, tree) ) 
  
  BW.compare <- llply(nodes[node.mask], function(node) DumBW.compare(means = getMeans(mean.list, tree, node), 
                                                                     contrasts = getContrasts(ind.cont, tree, node), 
                                                                     W.mx = cov.matrices[[node]])
  )
  names(BW.compare) <- nodes[node.mask]
  test.list.cor <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getMeans(mean.list, tree, node), 
                                                                      taxons = names(getMeans(mean.list, tree, node)), 
                                                                      cov.matrix =  W.fixed, 
                                                                      show.plots = FALSE,
                                                                      title.plot =  node) )
  names(test.list.cor) <- nodes[node.mask]
  
  test.list.reg <- llply(nodes[node.mask], function(node) DriftTest0(means = getMeans(mean.list, tree, node), 
                                                                               cov.matrix = W.fixed, 
                                                                               show.plot = FALSE))
  names(test.list.reg) <- nodes[node.mask]
  
  return(list ("Correlation.W.fixed" = test.list.cor,
               "Regression.W.fixed" = test.list.reg,
               "BW.compare" = BW.compare
) )
}

Drift.results$fixed <- DriftTestInterpretation(tree = Trees$all.with.ed, 
                        mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                        cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                        W.fixed = NULL,
                        #W.fixed = ancestral.mx$`42`,
                        sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )

Drift.results$fixed <- DriftTestInterpretation(tree = Trees$extant.sp.tree, 
                                               mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                                               cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                                               W.fixed = NULL,
                                               #W.fixed = ancestral.mx$Daubentonia_madagascariensis,
                                               sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree] )

Drift.results.Toplot$Fixed$Plots$Corr.Contrasts  <- Drift.results$fixed$Correlation.W.fixed %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 
Drift.results.Toplot$FixedW.Extants$Results$Node.ref <- Drift.results$fixed$Correlation.W.fixed %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,1]
Drift.results.Toplot$FixedW.Extants$Results$Corr.Ed.1 <- Drift.results$fixed$Correlation.W.fixed %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,2]
Drift.results.Toplot$FixedW.Extants$Results$Regression.Ed <- Drift.results$fixed$Regression.W.fixed %>% ldply(function(x) x$drift_rejected )  %>% .[,2]

Plot.Drift.Results(tree = Trees$all.with.ed, 
                   info = Drift.results.Toplot$All.sp$Results, 
                   font.size = 0.7,
                   contrasts = F,
                   paleta.deriva = c("#ECA400", "#1BAAB7", "#CBC9E0"))

Plot.Drift.Results(tree = Trees$extant.sp.tree, 
                   info = Drift.results.Toplot$FixedW.Extants$Results, 
                   font.size = 0.7,
                   contrasts = F,
                   paleta.deriva = c("#ECA400", "#1BAAB7", "#CBC9E0"))



Drift.results.Toplot$Fixed$Plots$Corr.Ed.1  <- Drift.results$fixed$Correlation.W.fixed %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 

Drift.results$fixed$Correlation.W.fixed %>% llply(Plot.Drift.regression)

Drift.results$fixed$Correlation.W.fixed$`135`$P.value.plot

Drift.results.Toplot$Extants$Results$Regression <- Drift.results$extant.sp$Regression.test %>% llply(Plot.Drift.regression)

Correlation.Tree.plot.Wfixed <- plot_grid(
  Drift.results$fixed$Correlation.W.fixed$`131`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lorisiformes\n Galagidae x Lorisidae") ,
  Drift.results$fixed$Correlation.W.fixed$`132`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lorisidae") ,
  Drift.results$fixed$Correlation.W.fixed$`134`$P.value.plot + theme(legend.position = "none") + ggtitle ("Nyc-Lor"), 
  
  Drift.results$fixed$Correlation.W.fixed$`99`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuridae x Indridae") ,
  Drift.results$fixed$Correlation.W.fixed$`100`$P.value.plot + theme(legend.position = "none") + ggtitle ("Indridae") ,
  Drift.results$fixed$Correlation.W.fixed$`112`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuridae") ,
  
  
  Drift.results$fixed$Correlation.W.fixed$`76`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lepilemuridae x Cheirogaleidae") ,
  Drift.results$fixed$Correlation.W.fixed$`77`$P.value.plot + theme(legend.position = "none") + ggtitle ("Cheirogaleidae") ,
  Drift.results$fixed$Correlation.W.fixed$`88`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lepilemuridae") ,
  
  Drift.results$fixed$Correlation.W.fixed$`73`$P.value.plot + theme(legend.position = "none") + ggtitle ("Strepsirrhini") ,
  Drift.results$fixed$Correlation.W.fixed$`74`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuriformes x Daubentonidae") ,
  Drift.results$fixed$Correlation.W.fixed$`75`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuriformes") ,
  Drift.results$fixed$Correlation.W.fixed$`71`$P.value.plot + theme(legend.position = "none") + ggtitle ("Prosimian\n Strepsirrhini + Tarsiidae") ,
  ncol = 3)

Correlation.Tree.plot.Wfixed
save_plot(filename = "Figures/Correlation_Tree_plot_Wfixed.pdf", plot = Correlation.Tree.plot.Wfixed, 
          base_aspect_ratio = 0.3, base_height = 15, base_width = 9.5)












GetBW<- function (tree, mean.list, cov.matrix.list, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(mean.list))) 
    stop("All tip labels must be in names(mean.list).")
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  if (!any(node.mask)) 
    stop("For the regression test:/nAt least one node must have more than 4 descendents in mean.list")
  
  ind.cont<- ldply(mean.list, function(x) x) 
  rownames(ind.cont) <- ind.cont[,1]
  ind.cont <- ind.cont[,-1]
  ind.cont <- apply(ind.cont, 2, FUN = function (x) ape::pic(x, tree) ) 
  
  BW.compare <- llply(nodes[node.mask], function(node) DumBW.compare(means = getMeans(mean.list, tree, node), 
                                                                     contrasts = getContrasts(ind.cont, tree, node), 
                                                                     W.mx = cov.matrices[[node]])  )
  names(BW.compare) <- nodes[node.mask]
  
  test.list.reg <- llply(nodes[node.mask], function(node) DriftTest0(means = getMeans(mean.list, tree, node), 
                                                                     cov.matrix = cov.matrices[[node]], 
                                                                     show.plot = FALSE))
  names(test.list.reg) <- nodes[node.mask]
  
  return(list ("Observed.Regression.Test" = test.list.reg,
               "BW.compare" = BW.compare
  ) )
}

BWtoSimulate <- vector("list")

BWtoSimulate$All.sp <- GetBW(tree = Trees$all.with.ed, 
                             mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                             cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                             sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )

BWtoSimulate$Extants <- GetBW(tree = Trees$extant.sp.tree, 
                              mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                              cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                              sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree] )


P <- BWtoSimulate$Extants$BW.compare$`71`$B.ed
G <- BWtoSimulate$Extants$BW.compare$`71`$W
pop <- BWtoSimulate$Extants$BW.compare$`71`$B.sample.size


