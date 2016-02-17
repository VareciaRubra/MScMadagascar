############################# Vendo efeito de filogenia ##########################################
notat.tree<- is.na(match(dimnames(mx.compare$RS$correlations)[[1]], treefile$tip.label)) 
names.at.tree <- dimnames(mx.compare$RS$correlations)[[1]][!notat.tree] 
mx.all.at.tree<- mx.compare$RS$correlations[!notat.tree,!notat.tree]
pruned.tree.all<- drop.tip(treefile,treefile$tip.label[-match( names.at.tree, treefile$tip.label)])
phylo.dist.all.at.tree<- cophenetic.phylo(pruned.tree.all)
phylo.dist.all.at.tree<- phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)]
MatrixCor(phylo.dist.all.at.tree, mx.all.at.tree)
### bom indicativo de que quanto mais próximo na filogenia mais semelhante sao as matrizes
#### vale a pena fazer um intervalo de confiança dessa porqueira?

##################################################################################################

################## Ponderar os valores de comparação pelo tamanho amostral : média armonica