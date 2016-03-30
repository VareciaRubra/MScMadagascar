PCScoreCorr <- 
function (means, cov.matrix, taxons = names(means), show.plots = FALSE) 
{
  if (is.data.frame(means) | (!is.array(means) & !is.list(means))) stop("means must be in a list or an array.")
  if (!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if (is.list(means)) mean.array <- laply(means, identity) else  mean.array <- means
  if (is.null(taxons)) taxons <- names(means) 
  n.taxon <- length(taxons) - 1
  #if (n.taxon >= dim(cov.matrix)[1]) n.taxon <- dim(cov.matrix)[1]
  
  n.taxon <-dim(cov.matrix)[1]
  W.pc <- eigen(cov.matrix)
  proj.med <- as.matrix(mean.array) %*% W.pc$vectors
  colnames(proj.med) <- paste0("PC", 1:ncol(proj.med))
  mat.pcs.deriva <- matrix(NA, ncol = n.taxon, nrow = n.taxon)
  if (show.plots) {
    plots <- vector("list", (n.taxon^2 - n.taxon)/2)
    current.plot <- 1
  }
  for (i in 1:n.taxon) {
    for (j in 1:i) {
      if (j != i) {
        test = cor.test(proj.med[, i], proj.med[, j])
        mat.pcs.deriva[i, j] <- test$estimate # Lower triangle of outputput 
        mat.pcs.deriva[j, i] <- test$p.value # Upper triangle are p.values
        if (show.plots == TRUE) {
          plots[[current.plot]] <- 
            ggplot(data.frame(x = proj.med[, i], y = proj.med[, j], taxons = taxons), aes_string("x", "y")) + 
            geom_text(aes_string(label = "taxons")) + 
            geom_smooth(method = "lm", color = "black") + 
            labs(x = paste0("PC", i), y = paste0("PC", j)) + 
            theme_bw() + 
            ggtitle(paste (taxons, "\ncorr=",round(mat.pcs.deriva[i, j], 2), "p=",round(mat.pcs.deriva[j, i], 4), sep =  " " ))
          current.plot <- current.plot + 1
        }
      }
    }
  }
  mx.to.bartlett <- mat.pcs.deriva
  mx.to.bartlett[upper.tri(mx.to.bartlett)] <- t(mx.to.bartlett)[upper.tri(mx.to.bartlett)]
  diag(mx.to.bartlett) <- 1
  Bartlett.t <- psych::cortest.bartlett (R = mx.to.bartlett, n= length(taxons))
  mx.to.bonferroni <- mat.pcs.deriva
  mx.to.bonferroni[lower.tri(mx.to.bonferroni)] <- t(mx.to.bonferroni)[lower.tri(mx.to.bonferroni)]
  colnames(mx.to.bonferroni) <- paste0("PC", 1:ncol(mx.to.bonferroni))
  rownames(mx.to.bonferroni) <- paste0("PC", 1:nrow(mx.to.bonferroni))
  n.veiz <- ( nrow(mx.to.bonferroni)^2 - nrow(mx.to.bonferroni) )/2
  
  #Correction.p.bonferroni <- corr.p(r = mx.to.bonferroni, n = n.veiz, adjust = "bonferroni") 

  mx.bonferroni <- mx.to.bonferroni
  mx.bonferroni[lower.tri(mx.bonferroni)] <- mx.bonferroni[lower.tri(mx.bonferroni)] > 0.05
  mx.bonferroni[upper.tri(mx.bonferroni)] <- mx.bonferroni[upper.tri(mx.bonferroni)] > 0.05/n.veiz
  Correction.p.bonferroni <- mx.bonferroni 
  Correction.p.bonferroni <- as.logical(Correction.p.bonferroni)
  Correction.p.bonferroni <- matrix(Correction.p.bonferroni, nrow = nrow(mx.to.bonferroni), ncol = ncol(mx.to.bonferroni))
  rownames(Correction.p.bonferroni) <- paste0('PC', 1:nrow(mx.to.bonferroni))
  colnames(Correction.p.bonferroni) <- paste0('PC', 1:ncol(mx.to.bonferroni))
  Correction.p.bonferroni[upper.tri(Correction.p.bonferroni)] <- NA
  
  rejected.drift<- which(!Correction.p.bonferroni, arr.ind = T)
  mx.bonferroni[mx.bonferroni == TRUE] <- "Not significant"
  mx.bonferroni[mx.bonferroni == 0] <- "Significative"
  mx.bonferroni %<>% melt 
  
  p.value.plot <- mx.bonferroni %>% 
    ggplot () +
      geom_tile(aes(x = Var2, y = Var1, fill = as.factor(value) ), alpha = 0.6, color = "darkgrey") +
      scale_y_discrete(limits = rev(levels(mx.bonferroni$Var1))) +
      #geom_text(aes(x = Var2, y = Var1, label = value), size = 4) +
      scale_fill_manual(values = c("lightgrey", "red")) + labs(fill = "p.value") +
      ylab ('') + xlab ('') + labs(title = "Original and Bonferroni correction" ) + 
      theme_minimal() +  
      theme(plot.title = element_text(face = "bold", size = 20),
            axis.text.x = element_text(angle = 270, hjust = 0, size =9),
            axis.ticks = element_line(size = 0),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10) ,
            rect = element_blank(), 
            line = element_blank())
  if (show.plots == TRUE) 
    return(list("Correlation.p.value" = mat.pcs.deriva, 
                "Bartlett" = Bartlett.t, 
                "Bonferroni" = Correction.p.bonferroni,
                "P.value.plot" = p.value.plot,
                "Resume.table" = rejected.drift,
                "plots" = plots))
  
  else return(list("Correlation.p.value" = mat.pcs.deriva, 
              "Bartlett" = Bartlett.t, 
              "Bonferroni" = Correction.p.bonferroni,
              "Resume.table" = rejected.drift,
              "P.value.plot" = p.value.plot
              ) )
}

TreeDriftTestPCScoresCorr <- function (tree, mean.list, cov.matrix.list, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(mean.list))) 
    stop("All tip labels must be in names(mean.list).")
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  
  if (!any(node.mask)) 
   stop("At least one node must have more than 4 descendents in mean.list")
  test.list <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getMeans(mean.list, tree, node), 
                                                                  taxons = names(nodes), 
                                                                          cov.matrix =  cov.matrices[[node]], 
                                                                  show.plots = TRUE))
  names(test.list) <- nodes[node.mask]
  return(test.list)
}

getMeans <- function(mean.list, tree, node){
  mean.list[na.omit(tree$tip.label[getDescendants(tree, node)])]
}

getContrasts <- function(contr, tree, node){
  contr[unlist(dimnames(contr)[1]) %in% as.character(getDescendants(tree = tree, node = node, curr = node) ),  ]
}


TreeDriftTestAll <- function (tree, mean.list, cov.matrix.list, sample.sizes = NULL) 
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
                                                         W.mx = cov.matrices[[node]])
                                                          )
    names(BW.compare) <- nodes[node.mask]
    test.list.cor <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getMeans(mean.list, tree, node), 
                                                              taxons = names(nodes), 
                                                              cov.matrix =  cov.matrices[[node]], 
                                                              show.plots = FALSE))
    names(test.list.cor) <- nodes[node.mask]
    
    test.list.cor.contrasts <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getContrasts(ind.cont, tree, node), 
                                                             taxons = rownames(getContrasts(ind.cont, tree, node)), 
                                                             cov.matrix =  cov.matrices[[node]], 
                                                             show.plots = FALSE))
    
    names(test.list.cor.contrasts) <-  nodes[node.mask]
    test.list.reg <- llply(nodes[node.mask], function(node) DriftTest(means = getMeans(mean.list, tree, node), 
                                                                      cov.matrix = cov.matrices[[node]], 
                                                                      show.plot = FALSE))
    names(test.list.reg) <- nodes[node.mask]
    
    test.list.reg.contrasts <- llply(nodes[node.mask], function(node) DriftTest(means = getContrasts(ind.cont, tree, node), 
                                                                      cov.matrix = cov.matrices[[node]], 
                                                                      show.plot = FALSE))
    names(test.list.reg.contrasts) <- nodes[node.mask]
    
    return(list ("Correlation.test.Regular" = test.list.cor,
                 "Correlation.test.Contrasts" = test.list.cor.contrasts, 
                 "Regression.test" = test.list.reg,
                 "Regression.test.Contrasts" = test.list.reg.contrasts,
                 "BW.compare" = BW.compare) )
    }

  Drift.results <- vector("list", 3)
  
  Drift.results$with.mx <- TreeDriftTestAll (tree = pruned.tree.with.mx, 
                                             mean.list = ed.means[mask][-41], 
                                             cov.matrix.list = cov.list[-41], 
                                             sample.sizes = sample.size[-c(41, 43, 44)])
  
  

drift.vai.porra <- vector("list", 5)

drift.vai.porra$drift.T.F.vectors <- 
  cbind(Drift.alltests.tree$Correlation.test.Regular %>% ldply(function(x) dim(x$Resume.table)[[1]] >1 ),
        Drift.alltests.tree$Correlation.test.Contrasts %>% ldply(function(x) dim(x$Resume.table)[[1]] >1 ) %>% .[,2],
        Drift.alltests.tree$Regression.test %>% ldply(function(x) x$drift_rejected ) %>% .[,2],
        Drift.alltests.tree$Regression.test.Contrasts %>% ldply(function(x) x$drift_rejected )  %>% .[,2])

colnames(drift.vai.porra$drift.T.F.vectors) <- c("node", "cor", "cor.ci", "reg", "reg.ci")
str(drift.vai.porra$drift.T.F.vectors)
drift.vai.porra$drift.T.F.vectors$node <- as.numeric(drift.vai.porra$drift.T.F.vectors$node)

par(mfrow = c(1,1))
plot.phylo(pruned.tree.with.mx, font = 3, no.margin = T)
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 8, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$cor)+1), frame = "n")
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 17, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$reg)+3), frame = "n")

plot.phylo(pruned.tree.with.mx, font = 3, no.margin = T)
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 8, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$cor.ci)+1), frame = "n")
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 17, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$reg.ci)+3), frame = "n")

drift.vai.porra$sum.abs.values <- Drift.alltests.tree$Correlation.test.Regular %>% llply(function (x) x$Correlation.p.value[1:39,1:39]) %>% llply(abs) %>% laply( function (x) x[lower.tri(x)]) %>% colSums 
temp <- matrix(NA, 39, 39, byrow = F)
temp[lower.tri(temp)] <- drift.vai.porra$sum.abs.values
drift.vai.porra$mean.abs.values <- temp / length(Drift.alltests.tree$Correlation.test.Regular)



contrasts <- apply(contrasts, 2, FUN = function (x) ape::pic(x, pruned.tree.with.mx) ) 
B.var <- sp.main.data[mask][-41] %>% ldply(function(x) x$ed.means) %>% .[, -1] %>% var
MatrixCompare(var(contrasts), Ancestral.Matrices$`42`)


DumBW.compare <- function (means, contrasts, W.mx) {
  B.ed <- var(laply(means, identity) ) 
  B.ic <- var(as.matrix(contrasts)  ) 
  BW.compare <- data.frame (BedBic = MatrixCompare(B.ed, B.ic)$correlation,
                           BedW = MatrixCompare(B.ed, W.mx)$correlation,
                           BicW = MatrixCompare(B.ic, W.mx)$correlation , 
                           method = c("RS", "Mantel", "KRZ", "PCA"), 
                           row.names = c("RS", "Mantel", "KRZ", "PCA") )
  
  return (list ("BW.compare" = BW.compare,
                "B.ed" = B.ed,
                "B.ic" = B.ic,
                "W" = W.mx) )
  
}

DumBW.compare(means = ed.means.with.mx[-41], contrasts = contrasts, W.mx = ancestral.mx$`45`)
