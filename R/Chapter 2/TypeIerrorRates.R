DriftTestInterpretation <- function (tree, mean.list, cov.matrix.list, W.fixed, sample.sizes = NULL) 
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
                                                                      cov.matrix =  W.fixed, 
                                                                      show.plots = FALSE,
                                                                      title.plot =  node) )
  names(test.list.cor) <- nodes[node.mask]
  
  test.list.reg <- llply(nodes[node.mask], function(node) DriftTest0(means = getMeans(mean.list, tree, node), 
                                                                               cov.matrix = W.fixed, 
                                                                               show.plot = FALSE))
  names(test.list.reg) <- nodes[node.mask]
  
  beta <- vector()
  prob <- vector()
  #generate samples from pop populations with n obs each. the ancestral vector is composed of m zeros.
  M <-mvrnorm(pop, rep(0, ncol(G)), tNe*G)
  group <- factor( rep( seq(1:pop), each = n))
  data <- matrix(0, pop*n, ncol(P) )
  #for (j in 1: (n*pop) )
  #{
  #  data[j,] <- mvrnorm(1, M[group[j],], P)
  #}
  
  data<-adply(1: (n*pop), 1, function(j) mvrnorm(1, M[group[j],], P))[,-1]
  #calculate matrix of mean vectors from simulations
  
  zmeans <- matrix(unlist(by(data, group, colMeans)),  nrow=pop, ncol = ncol(G), byrow=T)
  
  #calculate within-group phenotypic covariance and extract eigenvalues
  
  eigW <- eigen(cov( data - zmeans[rep(1:nrow(zmeans), each = n),]))
  
  #Project mean vectors for the pop pops on within group eigenvectors
  
  zm.proj <- zmeans %*% eigW$vectors
  
  #calculate among-group variance
  
  v <- diag(cov(zm.proj))
  
  #perform Ackermann and Cheverud test with t-test for unity slope
  
  model <- summary(lm(log(v) ~ log(eigW$values)))
  prob[i] <- pt(abs((model$coefficients[2,1]-1) / model$coefficients[2,2]), df = model$df[2], lower.tail=F)*2
  beta[i] <- model$coefficients[2,1]
}
#Calculate 95% confidence limits
Beta <- data.frame(Lower = NA, Upper = NA, Mean = NA)
Beta$Upper<-mean(beta)+1.96*sd(beta)
Beta$Lower<-mean(beta)-1.96*sd(beta)
Beta$Mean <- mean(beta)

tIe <- mean(prob < 0.05)
  
  return(list ("Correlation.W.fixed" = test.list.cor,
               "Regression.W.fixed" = test.list.reg,
               "BW.compare" = BW.compare,
               "TypeIerror" = tIe,
               "Probabilities" = prob,
               "Beta.ic" = Beta,
               "betas.dist" = beta) )
}

DriftTestInterpretation(tree = Trees$all.with.ed, 
                        mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                        cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                        W.fixed = ancestral.mx$`42`,
                        sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )
