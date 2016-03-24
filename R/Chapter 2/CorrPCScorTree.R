PCScoreCorr <- 
function (means, cov.matrix, taxons = names(means), show.plots = FALSE) 
{
  if (is.data.frame(means) | (!is.array(means) & !is.list(means))) stop("means must be in a list or an array.")
  if (!isSymmetric(cov.matrix)) stop("covariance matrix must be symmetric.")
  if (is.list(means)) mean.array <- laply(means, identity) else  mean.array <- means
  if (is.null(taxons)) taxons <- names(means) 
  n.taxon <- length(taxons) - 1
  if (n.taxon >= dim(cov.matrix)[1]) n.taxon <- dim(cov.matrix)[1]
  
  W.pc <- eigen(cov.matrix)
  proj.med <- as.matrix(mean.array) %*% W.pc$vectors
  colnames(proj.med) <- paste0("PC", 1:ncol(proj.med))
  mat.pcs.deriva <- matrix(1, ncol = n.taxon, nrow = n.taxon)
  if (show.plots) {
    plots <- vector("list", (n.taxon^2 - n.taxon)/2)
    current.plot <- 1
  }
  for (i in 1:n.taxon) {
    for (j in 1:i) {
      if (j != i) {
        test = cor.test(proj.med[, i], proj.med[, j])
        mat.pcs.deriva[i, j] <- test$estimate
        mat.pcs.deriva[j, i] <- test$p.value
        if (show.plots == TRUE) {
          plots[[current.plot]] <- 
            ggplot(data.frame(x = proj.med[, i], y = proj.med[, j], taxons = taxons), aes_string("x", "y")) + 
            geom_text(aes_string(label = "taxons")) + 
            geom_smooth(method = "lm", color = "black") + 
            labs(x = paste0("PC", i), y = paste0("PC", j)) + theme_bw() + 
            ggtitle(paste (taxons, "\ncorr=",round(mat.pcs.deriva[i, j], 2), "p=",round(mat.pcs.deriva[j, i], 4), sep =  " " ))
          current.plot <- current.plot + 1
        }
      }
    }
  }
  mx.to.bartlett <- mat.pcs.deriva
  mx.to.bartlett[upper.tri(mx.to.bartlett)] <- t(mx.to.bartlett)[upper.tri(mx.to.bartlett)]
  Bartlett.t <- psych::cortest.bartlett (R = mx.to.bartlett, n= n.taxon)
  
  if (show.plots == TRUE) 
    return(list(correlation_p.value = mat.pcs.deriva, Bartlett.test = Bartlett.t, plots = plots))
  else return(mat.pcs.deriva)
}


TreeDriftTestPCScoresCorr <- function (tree, mean.list, cov.matrix.list, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(mean.list))) 
    stop("All tip labels must be in names(mean.list).")
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, 
                                                        tree, x))) > 3
  
  if (!any(node.mask)) 
    stop("At least one node must have more than 4 descendents in mean.list")
  test.list <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getMeans(mean.list, 
                                                                         tree, node), 
                                                                  taxons = names(nodes), 
                                                                  cov.matrix =  cov.matrices[[node]], 
                                                                  show.plots = TRUE))
  names(test.list) <- nodes[node.mask]
  return(test.list)
}

getMeans <- function(mean.list, tree, node){
  means <- mean.list[na.omit(tree$tip.label[getDescendants(tree, node)])]
}

plot(pruned.tree.with.mx)
node

corr.drift.test <- PCScoreCorrelation(means = ed.means[mask], cov.matrix = W.matrix, taxons = names(ed.means[mask]), show.plots = T)

corr.drift.test.tree <- TreeDriftTestPCScoresCorr (tree = pruned.tree.with.mx, mean.list = ed.means[mask][-41], cov.matrix.list = cov.list[-41], sample.sizes = sample.size[-c(41, 43, 44)])
corr.drift.test.tree$`42`$plots[4]

TreeDriftTest(tree = pruned.tree.with.mx, mean.list = ed.means[mask][-41], cov.matrix.list = cov.list[-41], sample.sizes = sample.size[-c(41, 43, 44)])



corr.drift.test.tree %>% llply(function(x) x$Bartlett.test$p.value)


