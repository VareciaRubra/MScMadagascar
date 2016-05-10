##### Partition of variation: 
#### referencias: Yves Desdevises et al 2003: Evolution: QUANTIFYING PHYLOGENETICALLY STRUCTURED ENVIRONMENTAL VARIATION
####              Legendre, 1998: Livro: Numerical ecology

cov.mx <- sp.main.data %>% llply(function(x) x$matrix$cov)
cor.mx <- sp.main.data %>% llply(function(x) x$matrix$cor)
mx.rep <- sp.main.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
mx.rep.mc <- sp.main.data.t%>% ldply(function(x) x$Mx.Rep$MCRep) 
mx.rep.all <- cbind(mx.rep, mx.rep.mc[,-1])
n.size <- sp.main.data %>% ldply(function(x) x$sample.size)

mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
names (cov.mx[mask])

#1# Verificar se as matrizes de resultados de comparação sao euclidianas: ade4::is.euclid()
mx.compare %>% llply(function(x) x$correlations) %>% llply(., .fun = dist) %>% llply(., .fun = is.euclid)
Steppan[1:2] %>% llply(., .fun = dist) %>% llply(., .fun = is.euclid)
#ok, todas é!
# vendo se as matrizes sao simetricas
mx.compare %>% llply(function(x) x$correlations) %>%  llply(., .fun = isSymmetric ) 
### a matriz nao é simetrica pois sao os outputs das funções de comparação do pacote.
#mas o defaut da função dist é pegar a lower tri (onde tao os valores nao corrigidos das comparações, entao ta tudo bem)


#2# Transformar os valores de similaridade de matrizes em valores de diferença usando sqrt(1-semelhança)
mx.dissimilarity <- mx.compare %>% llply(function(x) x$correlations) %>% llply( function (x) sqrt(1-x) )
mx.dissimilarity$Steppan.RS <- sqrt(abs(1-Steppan$RS) )
mx.dissimilarity$Steppan.KRZ <- sqrt(abs(1-Steppan$KRZ) )

names(mx.dissimilarity)
dimnames(mx.dissimilarity$BS.RS)

# tirando a P e a G de Saguinus da lista de comparações:
mx.dissimilarity %<>% llply(function(x) x[-c(41,43:44), -c(41,43:44)] )
dimnames(mx.dissimilarity$BS.RS)

#3# Arrumando a árvore de quem entra na análise: pruned.tree.with.mx
# treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
plot(treefile, cex = 0.5)
nodelabels(bg = NULL, frame = "none")

notat.tree <- is.na(match(dimnames(mx.dissimilarity$BS.RS)[[1]], treefile$tip.label)) 
names.at.tree <- dimnames(mx.dissimilarity$BS.RS)[[1]][!notat.tree] 
mx.all.at.tree<- mx.compare$RS$correlations[!notat.tree,!notat.tree]
pruned.tree.with.mx<- drop.tip(treefile,treefile$tip.label[-match( names.at.tree, treefile$tip.label)])
phylo.dist<- cophenetic.phylo(pruned.tree.with.mx)
phylo.dist<- phylo.dist[rownames( mx.all.at.tree), rownames( mx.all.at.tree)]

dimnames(mx.all.at.tree)[[1]] %in% dimnames(mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])[[1]]

#4# Gerando a matriz ancestral
mx.at.tree <- cov.mx[mask][-41]
ancestral.mx <- PhyloW(tree = pruned.tree.with.mx, tip.data = mx.at.tree, tip.sample.size = n.size[mask,2][-41])
plot(pruned.tree.with.mx, cex = 0.5)
nodelabels()
W.matrix <- ancestral.mx$'42'
#5# evolqg::MultiMahalanobis distancia de Mahalanobis rntre as médias multivariadas das especies
multi.means <- sp.main.data %>% llply(function(x) x$ed.means)
mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
multi.means<- multi.means[mask]

MahalaLemur<- MultiMahalanobis(means =multi.means, cov.matrix = W.matrix)
MahalaLemur %>%  dist %>% is.euclid

#6# pegando os valores médios das distribuiçoes de R2 geradas pela evolqg::BootstrapR2()
r2.means <- sp.main.data %>% llply(function (x) mean(x$BootsR2$ed.means) )
mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
r2.means <- r2.means[mask][-41]
r2.means %>% dist %>% is.euclid
r2.means <- as.matrix(dist(r2.means, diag = T))
isSymmetric(r2.means)


#7# montando as PCOA de tudo essas coisa ae: ape:: pcoa
mx.pcoa <- mx.dissimilarity %>% llply(., .fun = dist) %>% llply(., rn = dimnames(mx.all.at.tree)[[1]], .fun = pcoa, .progress = "text") 
mx.pcoa$BS.RS$values
mx.pcoa$Steppan.RS$values
mx.pcoa$BS.KRZ$values

phylomorphospace(tree = pruned.tree.with.mx, X = mx.pcoa$MC.RS$vectors[,c(1,2)], label = "horizontal")

phylo.pcoa <- phylo.dist %>% dist %>% pcoa
phylo.pcoa$values
mahala.pcoa <- MahalaLemur %>% dist %>% pcoa
mahala.pcoa$values

r2.pcoa <- r2.means %>% dist %>% pcoa

phylomorphospace(tree = pruned.tree.with.mx, X = mahala.pcoa$vectors[,c(1,2)], label = "horizontal")
phylomorphospace(tree = pruned.tree.with.mx, X = phylo.pcoa$vectors[,c(1,2)], label = "horizontal")
phylomorphospace(tree = pruned.tree.with.mx, X = r2.pcoa$vectors[,c(1,2)], label = "horizontal")


# tentando jogar essas coisas no vegan
varpart(mx.pcoa$BS.RS$vectors[1:10], r2.pcoa$vectors[1:10], phylo.pcoa$vectors[1:10], mahala.pcoa$vectors[1:10]  )
rda(mx.pcoa$BS.RS$vectors[1:10], phylo.pcoa$vectors[1:10])
rda(mx.pcoa$BS.RS$vectors[1:10], mahala.pcoa$vectors[1:10])


# multi.mantel {phytools} Multiple matrix regression (partial Mantel test) ####
#Y	== single "dependent" square matrix. Can be either a symmetric matrix of class "matrix" or a distance matrix of class "dist".
#X	== a single independent matrix or multiple independent matrices in a list. As with Y can be a object of class "matrix" or class "dist".
dist.dissim <-  mx.dissimilarity$MC.KRZ
dist.dissim[upper.tri(dist.dissim)] <- NA
diag(dist.dissim) <- 0
dist.dissim <- dist(dist.dissim)

dist.sim <- mx.compare$BS.RS$correlations[-c(41,43:44), -c(41, 43:44)]
dist.sim[upper.tri(dist.sim)] <- NA
diag(dist.sim) <- 0
dist.sim <- dist(dist.sim)


Multi.Mantel <- list()
Multi.Mantel$diss.phylo <- multi.mantel(Y = dist.dissim, X = cophenetic.phylo(pruned.tree.with.mx), nperm=1000)
Multi.Mantel$diss.phylo$probF
Multi.Mantel$sim.phylo <- multi.mantel(Y = dist.sim, X = cophenetic.phylo(pruned.tree.with.mx), nperm=1000)
Multi.Mantel$sim.phylo$probF

