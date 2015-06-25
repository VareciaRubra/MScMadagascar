
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

# Árvore James: 
treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
species <- treefile$tip.label[treefile$tip.label %in% names(sample.no.na)]
pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)

tree.drift.test<- TreeDriftTest(tree = pruned.tree, mean.list = ed.means.no.na , cov.matrix.list = cov.no.na, sample.sizes = sample.no.na)
results <- llply(tree.drift.test, function(x) x$drift_rejected)
PlotTreeDriftTest(test.list = tree.drift.test, tree = pruned.tree)
nodelabels()
tree.drift.test$`36`+
coord_fixed()
table(ldply(tree.drift.test, function(x) x$drift_rejected))

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


