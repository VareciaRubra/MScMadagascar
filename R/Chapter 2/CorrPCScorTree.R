PCScoreCorr <- 
function (means, cov.matrix, taxons = names(means), show.plots = FALSE, title.plot = "Correlation Dtift Test") 
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
        mat.pcs.deriva[i, j] <- test$estimate # Lower triangle of outputput : vai receber os valores de correlaçao
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
  mx.to.bartlett <- mat.pcs.deriva # matriz que tem no triangulo superior os p.values das correlaçoes e no de baixo os valores observados de correlaçao
  mx.to.bartlett[upper.tri(mx.to.bartlett)] <- t(mx.to.bartlett)[upper.tri(mx.to.bartlett)] # transformando numa matriz toda de correla,coes
  signal.cor <- mx.to.bartlett # pegando o sinal da correlaçao
  signal.cor[lower.tri(signal.cor)] <- signal.cor[lower.tri(signal.cor)] > 0 # perguntando se a correla,cão é negativa
  signal.cor[upper.tri(signal.cor)] <- signal.cor[upper.tri(signal.cor)] > 0 # perguntando se a correla,cão é negativa
  diag(mx.to.bartlett) <- 1
  Bartlett.t <- psych::cortest.bartlett (R = mx.to.bartlett, n= length(taxons))
  mx.to.bonferroni <- mat.pcs.deriva
  mx.to.bonferroni[lower.tri(mx.to.bonferroni)] <- t(mx.to.bonferroni)[lower.tri(mx.to.bonferroni)]
  colnames(mx.to.bonferroni) <- paste0("PC", 1:ncol(mx.to.bonferroni))
  rownames(mx.to.bonferroni) <- paste0("PC", 1:nrow(mx.to.bonferroni))
  n.veiz <- ( nrow(mx.to.bonferroni)^2 - nrow(mx.to.bonferroni) )/2
  
  #Correction.p.bonferroni <- corr.p(r = mx.to.bonferroni, n = n.veiz, adjust = "bonferroni") 

  mx.bonferroni <- mx.to.bonferroni[1:10,1:10]
  mx.bonferroni[lower.tri(mx.bonferroni)] <- mx.bonferroni[lower.tri(mx.bonferroni)] < 0.05 # perguntando se o p.value é menor que 0.05 
  mx.bonferroni[upper.tri(mx.bonferroni)] <- mx.bonferroni[upper.tri(mx.bonferroni)] < 0.05/n.veiz # correçao de bonferroni 
  Correction.p.bonferroni <- mx.bonferroni 
  Correction.p.bonferroni <- as.logical(Correction.p.bonferroni)
  Correction.p.bonferroni <- matrix(Correction.p.bonferroni, nrow = nrow(mx.bonferroni), ncol = ncol(mx.bonferroni))
  rownames(Correction.p.bonferroni) <- paste0('PC', 1:nrow(mx.bonferroni))
  colnames(Correction.p.bonferroni) <- paste0('PC', 1:ncol(mx.bonferroni))
  Correction.p.bonferroni[lower.tri(Correction.p.bonferroni)] <- NA
  signal.cor <- signal.cor[1:10,1:10]
  rejected.drift<- which(!Correction.p.bonferroni, arr.ind = T)
  mx.bonferroni[mx.bonferroni == 1] <- "Not Significative"
  mx.bonferroni[mx.bonferroni == 0 & signal.cor == 0] <- "Significative (+)"
  mx.bonferroni[mx.bonferroni == 0 & signal.cor == 1] <- "Significative (-)"
  mx.bonferroni %<>% melt 
  
  p.value.plot <- mx.bonferroni %>% 
    ggplot () +
      geom_tile(aes(x = Var2, y = Var1, fill = value ), alpha = 0.6, color = "darkgrey") +
      scale_y_discrete(limits = rev(levels(mx.bonferroni$Var1))) +
      #geom_text(aes(x = Var2, y = Var1, label = value), size = 4) +
      scale_fill_manual(values = c("#f5f5f5", "#5ab4ac", "#f1a340")) + labs(fill = "p.value") +
      ylab ('') + xlab ('') + labs(title = paste("Node:", title.plot) ) + 
      theme_minimal() +  
      theme(plot.title = element_text(face = "bold", size = 17),
            axis.text.x = element_text(angle = 270, hjust = 0),
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

PCScoreCorr(means = All.sp.data$means, cov.matrix = ancestral.mx$`42`, show.plots = F)

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

DriftTest0 <- 
function (means, cov.matrix, show.plot = TRUE) 
{
  if (is.data.frame(means) | (!is.array(means) & !is.list(means))) 
    stop("means must be in a list or an array.")
  if (!isSymmetric(cov.matrix)) 
    stop("covariance matrix must be symmetric.")
  if (is.list(means)) {
    mean.array <- laply(means, identity)
  }
  else {
    mean.array <- means
  }
  W.pc <- eigen(cov.matrix)
  projection.Wpc <- as.matrix(mean.array) %*% W.pc$vectors
  log.B_variance <- log(apply(projection.Wpc, 2, var))
  log.W_eVals <- log(abs(W.pc$values))
  regression <- lm(log.B_variance ~ log.W_eVals)
  reg.plot <- ggplot(data.frame(log.B_variance, log.W_eVals, 
                                names = 1:(dim(mean.array)[2])), aes_string("log.W_eVals", 
                                                                            "log.B_variance")) + geom_text(aes_string(label = "names")) + 
    geom_smooth(method = "lm", color = "black") + labs(x = "log(W Eigenvalues)", 
                                                       y = "log(B variances)") + theme_bw()
  if (show.plot) 
    print(reg.plot)
  containsOne <- function(x) ifelse(x[1] < 1 & x[2] > 1, TRUE, 
                                    FALSE)
  test <- !containsOne(confint(regression)[2, ])
  names(test) <- "5 %"
  objeto <- list(regression = regression, 
                 coefficient_CI_95 = confint(regression), 
                 log.between_group_variance = log.B_variance, 
                 log.W_eVals = log.W_eVals, 
                 drift_rejected = test, 
                 plot = reg.plot,
                 means = means,
                 W.mx = cov.matrix)
  return(objeto)
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
                                                              taxons = names(getMeans(mean.list, tree, node)), 
                                                              cov.matrix =  cov.matrices[[node]], 
                                                              show.plots = FALSE,
                                                              title.plot =  node) )
    names(test.list.cor) <- nodes[node.mask]
    
    test.list.cor.contrasts <- llply(nodes[node.mask], function(node) PCScoreCorr(means = getContrasts(ind.cont, tree, node), 
                                                             taxons = rownames(getContrasts(ind.cont, tree, node)), 
                                                             cov.matrix =  cov.matrices[[node]], 
                                                             show.plots = FALSE,
                                                             title.plot =  node) )
    
    names(test.list.cor.contrasts) <-  nodes[node.mask]
    test.list.reg <- llply(nodes[node.mask], function(node) DriftTest0(means = getMeans(mean.list, tree, node), 
                                                                      cov.matrix = cov.matrices[[node]], 
                                                                      show.plot = FALSE))
    names(test.list.reg) <- nodes[node.mask]
    
    test.list.reg.contrasts <- llply(nodes[node.mask], function(node) DriftTest0(means = getContrasts(ind.cont, tree, node), 
                                                                      cov.matrix = cov.matrices[[node]], 
                                                                      show.plot = FALSE))
    names(test.list.reg.contrasts) <- nodes[node.mask]
    
    test.list.reg.contrasts.W.fixed <- llply(nodes[node.mask], function(node) DriftTest0(means = getContrasts(ind.cont, tree, node), 
                                                                                 cov.matrix = cov.matrices[[node]], 
                                                                                 show.plot = FALSE))
    names(test.list.reg.contrasts) <- nodes[node.mask]
    
    
    return(list ("Correlation.test.Regular" = test.list.cor,
                 "Correlation.test.Contrasts" = test.list.cor.contrasts, 
                 "Regression.test" = test.list.reg,
                 "Regression.test.Contrasts" = test.list.reg.contrasts,
                 "BW.compare" = BW.compare) )
}


  Drift.results <- vector("list")
  
  Drift.results$with.mx <- TreeDriftTestAll (tree = pruned.tree.with.mx, 
                                             mean.list = ed.means[mask][-41], 
                                             cov.matrix.list = cov.list[-41], 
                                             sample.sizes = sample.size[-c(41, 43, 44)])
  Drift.results$with.mx$Correlation.test.Regular$`42`$P.value.plot
  

drift.vai.porra$sum.abs.values <- 
  Drift.results$with.mx$Correlation.test.Regular %>% 
  llply(function (x) x$Correlation.p.value[1:39,1:39]) %>% 
  llply(abs) %>% laply( function (x) x[lower.tri(x)]) %>% colSums 
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
