##### Partition of variation: 
#### referencias: Yves Desdevises et al 2003: Evolution: QUANTIFYING PHYLOGENETICALLY STRUCTURED ENVIRONMENTAL VARIATION
####              Legendre, 1998: Livro: Numerical ecology

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

#6# montando as PCOA de tudo essas coisa ae: ape:: pcoa
mx.pcoa <- mx.dissimilarity %>% llply(., .fun = dist) %>% llply(., rn = dimnames(mx.all.at.tree)[[1]], .fun = pcoa, .progress = "text") 
mx.pcoa$BS.RS$values
mx.pcoa$Steppan.RS$values
mx.pcoa$BS.KRZ$values

phylo.pcoa <- phylo.dist %>% dist %>% pcoa
phylo.pcoa$values
mahala.pcoa <- MahalaLemur %>% dist %>% pcoa
mahala.pcoa$values
# tentando jogar essas coisas no vegan
varpart(mx.pcoa$BS.RS$vectors[1:10], phylo.pcoa$vectors[1:10], mahala.pcoa$vectors[1:10]  )
rda(mx.pcoa$BS.RS$vectors[1:10], phylo.pcoa$vectors[1:10])
rda(mx.pcoa$BS.RS$vectors[1:10], mahala.pcoa$vectors[1:10])


