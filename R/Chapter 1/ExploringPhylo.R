mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
eds.mx <- sp.main.data[mask] %>% ldply(., function (x) x$ed.means) 
rownames(eds.mx) <- eds.mx$.id
eds.mx.info <- data.frame (".id"= eds.mx$.id)
eds.mx.info$info <- eds.mx.info$.id
eds.mx.info <- eds.mx.info %>% separate(info, c('Genus', 'Species'), sep = "\\_")
eds.mx <- eds.mx[,-1]
eds.mx <- t(eds.mx) 

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


gm.means <- sp.main.data[mask] %>% ldply(., function (x) x$gm.mean) 
dimnames(gm.means)[[1]] <- gm.means$.id
at.tree <- !is.na(match(dimnames(gm.means)[[1]], pruned.tree.with.mx$tip.label)) 
gm.means <- gm.means[at.tree,]
gm.means$.id <- factor(gm.means$.id, levels = unique(gm.means$.id) ) 
gm.means$.id <- reorder.factor(gm.means$.id, new.order=pruned.tree.with.mx$tip.label)

gm.means <- gm.means%>% arrange(.id)

pruned.tree.with.mx$tip.label == as.character(gm.means$.id)
rownames(gm.means) <- 
Phylo.Sig <- phylosig(tree = pruned.tree.with.mx , x = gm.means[,-1],  method="K", test=T, nsim=1000, se=NULL, start=NULL)

tree.t<-pbtree(n=25)
X<-fastBM(tree.t,nsim=2)
phylomorphospace(tree.t,X,xlab="trait 1",ylab="trait 2")
