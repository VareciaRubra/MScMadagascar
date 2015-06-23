

Wmat.All.Fuckers = CalculateMatrix(manova(as.matrix(all.main.data$All$ed)  ~ Especie, data = as.data.frame(all.main.data$All$info) ) )
Wmat.extant.Fuckers = CalculateMatrix(manova(as.matrix(extant.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.main.data$All$info) ) )
Wmat.madagascar.Fuckers = CalculateMatrix(manova(as.matrix(madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(madagascar.main.data$All$info) ) )
Wmat.extant.madagascar.Fuckers = CalculateMatrix(manova(as.matrix(extant.madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.madagascar.main.data$All$info) ) )

mask.extant<- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Status != "Extinct"])
mask.madagascar<- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Regiao == "Madagascar"])
mask.madagascar.extant <- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Status != "Extinct" &All.raw.main.data$All$Regiao == "Madagascar"])

means.All.Fuckers<- sp.master.main.data %>% llply(function(x) x$ed.means )
means.extant.Fuckers<- means.All.Fuckers[mask.extant] 
means.madagascar.Fuckers<- means.All.Fuckers[mask.madagascar]
means.extant.madagascar.Fuckers<- means.All.Fuckers[mask.madagascar.extant]
  
DriftTest(means = means.All.Fuckers, cov.matrix = Wmat.All.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.Fuckers, cov.matrix = Wmat.extant.Fuckers, show.plot = TRUE)
DriftTest(means = means.madagascar.Fuckers, cov.matrix = Wmat.madagascar.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.madagascar.Fuckers, cov.matrix = Wmat.extant.madagascar.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.madagascar.Fuckers, cov.matrix = Wmat.madagascar.Fuckers, show.plot = TRUE)


cov.mx <- sp.master.main.data %>% llply(function(x) x$matrix$cov)
mx.rep <- sp.master.main.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) 
n.size <- sp.master.main.data %>% ldply(function(x) x$sample.size) 
ed.means <- sp.master.main.data %>% llply(function(x) x$ed.means) 
gm.mean <- sp.master.main.data %>% ldply(function(x) x$gm.mean) 

###Sample sizes #####
sample.no.na <- n.size[mask.na.cov,2]
names(sample.no.na) <- n.size[,1][mask.na.cov]
##### Cov mx ########
cov.no.na <- cov.mx[mask.na.cov]
##### Mx repetabilities via RS ######
rep.no.na <- mx.rep[mask.na.cov,2]
names(rep.no.na) <- mx.rep[,1][mask.na.cov]
##### ED means ######
ed.means <- ed.means[mask.na.cov]

treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
plot(treefile)
species <- treefile$tip.label[treefile$tip.label %in% names(sample.no.na)]

pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)

tree.drift.test<- TreeDriftTest(tree = pruned.tree, mean.list = ed.means , cov.matrix.list = cov.no.na, sample.sizes = sample.no.na)
results <- llply(tree.drift.test, function(x) x$drift_rejected)
PlotTreeDriftTest(test.list = tree.drift.test, tree = pruned.tree)
nodelabels()
tree.drift.test$`49`
