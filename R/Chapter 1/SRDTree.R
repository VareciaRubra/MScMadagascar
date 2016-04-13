


SRD.Tree <- function (tree, cov.matrix.list, sample.sizes = NULL) 
{
  if (!all(tree$tip.label %in% names(cov.matrix.list))) 
    stop("All tip labels must be in names(cov.matrix.list).")
  cov.matrices <- PhyloW(tree, cov.matrix.list, sample.sizes)
  
  nodes <- names(cov.matrices)
  node.mask <- laply(nodes, function(x) length(getMeans(mean.list, tree, x))) > 3
  
  newvec <- names(cov.matrices[1:70])
  for (i in unique(tree$tip.label)) newvec[names(cov.matrices[1:70]) == i ] <- c(1:70)[tree$tip.label == i ]
  names(cov.matrices)[1:70] <- newvec
  
  SRD.by.node <- llply(nodes[node.mask], 
                       function(node) 
                        {
                          indices <- as.character(tree$edge[,2] [tree$edge[,1] == node])
                          tryCatch (expr = SRD(cov.matrices [[ indices [1] ]], cov.matrices [[ indices [2] ]]) , error = function(cond) NA )
                        })
  names(SRD.by.node) <- nodes[node.mask]

  return(list ("SRD.by.node" = SRD.by.node) )
}


plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(cex = 0.7)

SRD.results <- SRD.Tree(tree = Trees$extant.sp.tree, 
                        cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree],
                        sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree])


plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(cex = 0.7)
plot(SRD.results$SRD.by.node$'100')
plot(SRD(mx.list.taxonomy$W.Indridae, mx.list.taxonomy$W.Lemuridae))

plot(SRD.results$SRD.by.node$'76')
plot(SRD(mx.list.taxonomy$W.Cheirogaleidae, mx.list.taxonomy$Lepilemur))

plot(SRD(cov.matrices[as.character(tree$edge[,2][tree$edge[,1] == 100])])[1,2][[1]] )


cov.matrices[as.character(tree$edge[,2][tree$edge[,1] == as.character(100)]) ]



as.character(tree$edge[,2][tree$edge[,1] == as.character(100)]) 

tryCatch(expr = SRD(cov.matrices$'1', cov.matrices$'2'), error = function(cond) NA)
