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
                                               #W.fixed = NULL,
                                               W.fixed = ancestral.mx$Daubentonia_madagascariensis,
                                               sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree] )

Drift.results.Toplot$Fixed$Plots$Corr.Contrasts  <- Drift.results$fixed$Correlation.W.fixed %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 
Drift.results.Toplot$FixedW.Extants$Results$Node.ref <- Drift.results$fixed$Correlation.W.fixed %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,1]
Drift.results.Toplot$FixedW.Extants$Results$Corr.Ed.1 <- Drift.results$fixed$Correlation.W.fixed %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,2]
Drift.results.Toplot$FixedW.Extants$Results$Regression.Ed <- Drift.results$fixed$Regression.W.fixed %>% ldply(function(x) x$drift_rejected )  %>% .[,2]

Plot.Drift.Results(tree = Trees$extant.sp.tree, 
                   info = Drift.results.Toplot$FixedW.Extants$Results, 
                   font.size = 0.4,
                   contrasts = F)

Drift.results.Toplot$Fixed$Plots$Corr.Ed.1  <- Drift.results$fixed$Correlation.W.fixed %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 

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


