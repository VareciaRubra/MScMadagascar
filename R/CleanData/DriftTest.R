
#Teste na raiz apenas, usando a matriz ancestral:
### nao roda pq? Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 0 (non-NA) cases
DriftTest(means = means.All.Fuckers, cov.matrix = Wmat.All.Fuckers, show.plot = TRUE) 
#rodam, e rejeita em todos:
DriftTest(means = means.madagascar.Fuckers, cov.matrix = Wmat.madagascar.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.madagascar.Fuckers, cov.matrix = Wmat.extant.madagascar.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.madagascar.Fuckers, cov.matrix = Wmat.madagascar.Fuckers, show.plot = TRUE)

current.data<- sp.main.data

# selecting parts of the main.data.list to perform the matrices comparrison and drift test 
cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) 
n.size <- current.data %>% ldply(function(x) x$sample.size) 
ed.means <- current.data %>% llply(function(x) x$ed.means) 
gm.mean <- current.data %>% ldply(function(x) x$gm.mean) 

# creating the masks to index the current set of matrices, means and repetabilities
mask.n.size <- n.size[,2]>40
mask.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
mask.rep <- !is.na(mx.rep)[,2]
#Selecting accordingly to masks
###Sample sizes #####
sample.no.na <- n.size[mask.na.cov,2]
names(sample.no.na) <- n.size[,1][mask.na.cov]
##### Cov mx ########
cov.no.na <- cov.mx[mask.na.cov]
##### Mx repetabilities via RS ######
rep.no.na <- mx.rep[mask.na.cov,2]
names(rep.no.na) <- mx.rep[,1][mask.na.cov]
##### ED means ######
ed.means.no.na <- ed.means[mask.na.cov]
gm.mean.no.na <- gm.mean [mask.na.cov, 2]
# Árvore James: 
treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
species <- treefile$tip.label[treefile$tip.label %in% names(sample.no.na)]
pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)

tree.drift.test<- TreeDriftTest(tree = pruned.tree, mean.list = ed.means.no.na , cov.matrix.list = cov.no.na, sample.sizes = sample.no.na)
results <- llply(tree.drift.test, function(x) x$drift_rejected)
PlotTreeDriftTest(test.list = tree.drift.test, tree = pruned.tree)
nodelabels()
tree.drift.test$`45`

######### changing plot 
tested.nodes <- as.numeric(names(tree.drift.test))
non.drift.nodes <- laply(tree.drift.test, function(x) x$drift_rejected)
i.c.5 <- ldply(tree.drift.test, function(x) x$coefficient_CI_95[2,])
i.c.5<- ldply(tree.drift.test, function(x) x$coefficient_CI_95[1,])
names(i.c.5) <- c('node', 'min', 'max')
names(intercept.ic) <- c('node', 'min', 'max')
i.c.5$min <- round(i.c.5$min, digits=2)
i.c.5$max <- round(i.c.5$max, digits=2)

# par(mfrow = c(1,1))
# par(oma=c(1,1,1,1))
# par(mar=c(0,2,1,1)+0.1)
plot(pruned.tree, font = 3, label.offset = 15, no.margin = T, edge.color = "grey80", edge.width = 5, cex= 1.1)
#nodelabels(node = tested.nodes, thermo = as.numeric(non.drift.nodes))
i.c.5 <- ldply(tree.drift.test, function(x) x$coefficient_CI_95[2,])
names(i.c.5) <- c('node', 'min', 'max')
names(intercept.ic) <- c('node', 'min', 'max')
i.c.5$min <- round(i.c.5$min, digits=2)
i.c.5$max <- round(i.c.5$max, digits=2)
nodelabels(node = tested.nodes , text = i.c.5$min, adj = 1.2, bg = "transparent", col = 'red', frame = "n", cex = i.c.5$min)
nodelabels(node = tested.nodes, i.c.5$max, adj = -0.2, bg = "transparent", col = 'blue', frame = "n", cex = i.c.5$max)
nodelabels(node = tested.nodes, pch = (as.numeric(non.drift.nodes)+17), cex=1.4, col = (as.numeric(non.drift.nodes)+5) )
compare.stuff<- cbind( n.size[n.size[,1] %in% pruned.tree$tip.label, ], gm.mean[gm.mean[,1] %in% pruned.tree$tip.label, 2])
names(compare.stuff) <- c('Especie', 'n', 'gm')
compare.stuff$Especie <- as.character(compare.stuff$Especie)
compare.stuff<- compare.stuff[with(compare.stuff, order(Especie)), ]
tiplabels(pch = 21, cex =compare.stuff$n/18, adj = 6)
tiplabels(pch = 15, cex = compare.stuff$gm/9, adj = 10)
nodelabels(node = tested.nodes, col = 'grey30', adj = c(0.5, 1.25), cex = 0.5, frame = "n")
# compare.stuff$Especie <- compare.stuff[sort(compare.stuff[,1]),]

par(mar=c(0,0,1,3)+0.1)
plot(pruned.tree, font = 3, label.offset = -3, no.margin = T, edge.color = "grey80", edge.width = 5, cex= 1.1, direction ="leftwards", show.tip.label = T, col = "white")
#nodelabels(node = tested.nodes, thermo = as.numeric(non.drift.nodes))
i.c.5<- ldply(tree.drift.test, function(x) x$coefficient_CI_95[1,])
names(i.c.5) <- c('node', 'min', 'max')
names(intercept.ic) <- c('node', 'min', 'max')
i.c.5$min <- round(i.c.5$min, digits=2)
i.c.5$max <- round(i.c.5$max, digits=2)
nodelabels(node = tested.nodes , i.c.5$min, adj = -0.2, bg = "transparent", col = 'red', frame = "n", cex = i.c.5$min*0.6)
nodelabels(node = tested.nodes, i.c.5$max, adj = 1.2, bg = "transparent", col = 'blue', frame = "n", cex = i.c.5$max*0.6)
nodelabels(node = tested.nodes, pch = (as.numeric(non.drift.nodes)+17), cex=1.4, col = (as.numeric(non.drift.nodes)+9) )
tiplabels(pch = 19, cex = gm.mean.no.na/10, adj = -2.5)



tree.drift.test$`47`
nodelabels()


mask.n.size <- n.size[,2]>40
###Sample sizes #####
sample.40 <- n.size[mask.n.size,2]
names(sample.40) <- n.size[mask.n.size,1]
##### Cov mx ########
cov.40 <- cov.mx[mask.n.size]
##### Mx repetabilities via RS ######
rep.40 <- mx.rep[mask.n.size,2]
names(rep.40) <- mx.rep[,1][mask.n.size]
##### ED means ######
ed.means.40 <- ed.means[mask.n.size]

names(ed.means.40) == names(rep.40) == names(cov.mx[mask.n.size]) == mx.rep[,1][mask.n.size]

species.40 <- treefile$tip.label[treefile$tip.label %in% names(sample.40)]
pruned.tree.40<-drop.tip(treefile,treefile$tip.label[-match(species.40, treefile$tip.label)])
plot(pruned.tree.40)

tree.drift.test.40<- TreeDriftTest(tree = pruned.tree.40, mean.list = ed.means.40 , cov.matrix.list = cov.40, sample.sizes = sample.40)
results <- llply(tree.drift.test, function(x) x$drift_rejected)
PlotTreeDriftTest(test.list = tree.drift.test.40, tree = pruned.tree.40)
nodelabels()


