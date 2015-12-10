
#Pegando a matriz ancestral de Indridae
Ancestral.Matrices$"56"

ed.means <- sp.main.data %>% llply(function(x) x$ed.means) 
cov.mx.sp <- sp.main.data %>% llply(function(x) x$matrix$cov)
cov.mx.genus <- gen.main.data %>% llply(function(x) x$matrix$cov)
n.size.sp <- sp.main.data %>% ldply(function(x) x$sample.size) 

#Pegando os caras da familia indridae e os extintos, que todos os que eu tenho tao dentro dela
indridae.family <- grepl("Indri", names(ed.means)) | grepl("Avahi", names(ed.means)) | grepl("Propithecus", names(ed.means))
extintos <- sp.main.data %>% laply(function(x) all(x$info$Status == "Extinct") )
extintos.names <- names(ed.means)[extintos]
viventes.names <- names(ed.means)[indridae.family]
#Pegando um vetor que vai selecionar somente esses caras
indroidea <- extintos | indridae.family
indroidea.names <- names(ed.means)[indroidea]

#indroidea só viventes
ed.means.indridae.family <- ed.means [indridae.family] 
ed.means.indroidea <- ed.means [indroidea] 

#pegando as médias de todos os de indroidea
indroidea.means <- sp.main.data[indroidea] %>% ldply(function (x) x$ed.means)
rowSums(is.na(indroidea.means[, -1]))

#arvore trimada pros indridae viventes e outra com os indroidea
treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
plot(treefile, cex = 0.4)
nodelabels()
pruned.indroidea <- extract.clade(treefile, 222)
plot(pruned.indroidea, cex = 1)

pruned.tree.indroidea<- drop.tip(treefile,treefile$tip.label[-match(indroidea.names, treefile$tip.label)])
plot(pruned.tree.indroidea)

pruned.tree.indridae<- drop.tip(treefile,treefile$tip.label[-match(viventes.names, treefile$tip.label)])
plot(pruned.tree.indridae)

cov.indridae <- cov.mx.sp[indridae.family]
mask.cov.na.indridae <- ldply(cov.indridae, function(x) !is.na(x[1]))[,2]
names(cov.indridae)[mask.cov.na.indridae == F]

cov.indridae$Avahi_meridionalis <- cov.mx.genus$Avahi
cov.indridae$Avahi_occidentalis <- cov.mx.genus$Avahi
cov.indridae$Propithecus_candidus <- cov.mx.genus$Propithecus
cov.indridae$Propithecus_coquereli <- cov.mx.genus$Propithecus
cov.indridae$Propithecus_tattersalli <- cov.mx.genus$Propithecus

n.size.indridae <- n.size.sp[indridae.family,2]
mask.cov.na.indridae <- ldply(cov.indridae, function(x) !is.na(x[1]))[,2]
n.size.indridae[mask.cov.na.indridae == F] <- 1
n.size.sp

Ancestral.Matrices.indridae<- PhyloW(tree = pruned.tree.indridae, 
                                     tip.data = cov.indridae, 
                                     tip.sample.size = n.size.indridae )
Ancestral.Matrices.indridae$"14"

drift.test.indridae <- TreeDriftTest(tree = pruned.tree.indridae, mean.list = ed.means.indridae.family , cov.matrix.list = cov.indridae, sample.sizes = n.size.indridae)
PlotTreeDriftTest(test.list = drift.test.indridae, tree = pruned.tree.indridae)
