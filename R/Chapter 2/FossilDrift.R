###################### BRINCANDO COM OS FOSSEIS E TESTE DE DERIVA ############################
#olhando pro pedaço da árvore que quero trabalhar:
#arvore trimada para Indroidea: diversidade total desse clado
treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
plot(treefile, cex = 0.4)
nodelabels()
pruned.indroidea.total <- extract.clade(treefile, 222)
plot(pruned.indroidea.total, cex = 1)

#separando os fosseis e viventes dessa familia indroidea apenas e os que tenho em minha amostra
indroidea <- names(ed.means) %in% pruned.indroidea$tip.label
indroidea.names <- names(ed.means)[indroidea]

###############
# 1. estimar matrizes pra tudo quanto é coisa que tem amostra >2
# 2. Rodar a Pholow(  ):
###### 3. PhyloW(tree = tree.file, tip.data = cov.matrices.low.sample.sizes, tip.sample.size = real.sample.sizes)
# 4. colocar nas matrizes que só tem um bicho a matriz do nó mais próximo.

############################ 1.  estimar matrizes pra tudo quanto é coisa que tem amostra >2
#Tirando quem tem algum NA na amostra

indroidea.data <- sp.main.data[indroidea]

na.na.na<- rowSums(is.na(indroidea.data %>% ldply(function (x) Matrix::colMeans(x$ed, na.rm = T) ))) !=0

cov.mx.all <- indroidea.data[!na.na.na] %>% llply(function (x) cov(x$ed))
ed.means.all <- indroidea.data[!na.na.na] %>% ldply(function (x) Matrix::colMeans(x$ed, na.rm = T) ) 
sample.size.all <- indroidea.data[!na.na.na] %>% ldply(function(x) x$sample.size) %>% .[,2]
names.all <- names(cov.mx.all)
pruned.tree.all<- drop.tip(treefile,treefile$tip.label[-match(names.all, treefile$tip.label)])
plot(pruned.tree.all)

PhyloW(tree = )




















#Olhando quem são os fosseis da minha amostra
extintos <- sp.main.data %>% laply(function(x) all(x$info$Status == "Extinct") )
extintos.names <- names(ed.means)[extintos]


#Pegando os caras da familia indridae e os extintos, que todos os que eu tenho tao dentro dela
indridae <- grepl("Indri", names(ed.means)) | grepl("Avahi", names(ed.means)) | grepl("Propithecus", names(ed.means))
indridae.names <- names(ed.means)[indridae]

#Arvore dos indroidea: viventes e extintos: apenas os que tenho em minha amostra
pruned.tree.indroidea<- drop.tip(treefile,treefile$tip.label[-match(indroidea.names, treefile$tip.label)])
plot(pruned.tree.indroidea)
#Árvore dos indridade: apenas viventes que tenho na minha amostra
pruned.tree.indridae<- drop.tip(treefile,treefile$tip.label[-match(indridae.names, treefile$tip.label)])
plot(pruned.tree.indridae)
indridae.names <- pruned.tree.indridae$tip.label

#Filtrando matrizes, médias e tamanhos amostrais pra este clado
ed.means <- sp.main.data %>% llply(function(x) x$ed.means) 
cov.mx.sp <- sp.main.data %>% llply(function(x) x$matrix$cov)
cov.mx.genus <- gen.main.data %>% llply(function(x) x$matrix$cov)
n.size.sp <- sp.main.data %>% ldply(function(x) x$sample.size) 

#indroidea só viventes
ed.means.indridae <- ed.means [indridae] 
#pegando as médias de todos os de indroidea
indroidea.means <- sp.main.data[indroidea] %>% llply(function (x) Matrix::colMeans(x$ed, na.rm = T)) 
indroidea.na<- indroidea.means %>% ldply(function(x) dim(table(is.na(x) ) ) !=1)
indroidea.names <- names(ed.means)[indroidea][!indroidea.na[,2]]
pruned.tree.indroidea<- drop.tip(treefile,treefile$tip.label[-match(indroidea.names, treefile$tip.label)])
plot(pruned.tree.indroidea)

cov.indridae <- cov.mx.sp[indridae]
mask.cov.na.indridae <- ldply(cov.indridae, function(x) !is.na(x[1]))[,2]
names(cov.indridae)[mask.cov.na.indridae == F]
n.size.indridae <- n.size.sp[indridae,2]

RandomSkewers(cov.x = list("Propithecus" = cov.mx.genus$Propithecus,
                           "P.candidus" = cov(sp.main.data$Propithecus_candidus$ed),
                           "P.coquereli" = cov(sp.main.data$Propithecus_coquereli$ed),
                           "P.verreauxi" = cov.indridae$Propithecus_verreauxi),
                            num.vectors = 1000)

cov.indridae$Propithecus_candidus <- cov(sp.main.data$Propithecus_candidus$ed)
cov.indridae$Propithecus_coquereli <- cov(sp.main.data$Propithecus_coquereli$ed)

cov.indridae$Avahi_meridionalis <- cov.mx.genus$Avahi
cov.indridae$Avahi_occidentalis <- cov.mx.genus$Avahi
cov.indridae$Propithecus_tattersalli <- cov.mx.genus$Propithecus
#cov.indridae$Propithecus_coquereli <- cov.mx.genus$Propithecus
#cov.indridae$Propithecus_candidus <- cov.mx.genus$Propithecus

n.size.indridae[mask.cov.na.indridae == F] <- 1



Ancestral.Matrices.indridae<- PhyloW(tree = pruned.tree.indridae, 
                                     tip.data = cov.indridae, 
                                     tip.sample.size = n.size.indridae )
Ancestral.Matrices.indridae$"14"

drift.test.indridae <- TreeDriftTest(tree = pruned.tree.indridae, mean.list = ed.means.indridae , cov.matrix.list = cov.indridae, sample.sizes = n.size.indridae)
PlotTreeDriftTest(test.list = drift.test.indridae, tree = pruned.tree.indridae)

cov.indroidea <- cov.mx.sp[indroidea][indroidea.no.na[,2]]
mask.cov.na.indroidea <- ldply(cov.indroidea, function(x) is.na(x[1]) )[,2]
names(cov.indroidea)[mask.cov.na.indroidea]

indroidea.means<- indroidea.means[!indroidea.na[,2]]
indroidea.sample.size <- sp.main.data[indroidea][!indroidea.na[,2]] %>% ldply(function (x) x$sample.size)
#indroidea.sample.size$V2 <- rowSums(is.na(indroidea.means[, -1]))
indroidea.sample.size<- indroidea.sample.size[,2]
indroidea.sample.size[mask.cov.na.indroidea]

cov.indroidea$Avahi_meridionalis <- cov.mx.genus$Avahi
cov.indroidea$Avahi_occidentalis <- cov.mx.genus$Avahi
cov.indroidea$Propithecus_candidus <- cov.mx.genus$Propithecus
cov.indroidea$Propithecus_coquereli <- cov.mx.genus$Propithecus
cov.indroidea$Propithecus_tattersalli <- cov.mx.genus$Propithecus

mask.cov.fossil.indroidea <- ldply(cov.indroidea, function(x) is.na(x[1]) ) [,2]
str(cov.indroidea[mask.cov.fossil.indroidea] )
cov.indroidea[mask.cov.fossil.indroidea] <- list(Ancestral.Matrices.indridae$"14")

n.size.indroidea <- n.size.sp[indroidea,2]
n.size.indroidea[!mask.cov.na.indroidea] <- 1

drift.test.indroidea <- TreeDriftTest(tree = pruned.tree.indroidea, mean.list = indroidea.means , cov.matrix.list = cov.indroidea, sample.sizes = n.size.indroidea)
PlotTreeDriftTest(test.list = drift.test.indroidea, tree = pruned.tree.indroidea)
nodelabels(bg = NULL)
Ancestral.Matrices.indridae<- PhyloW(tree = pruned.tree.indridae, 
                                     tip.data = cov.indridae, 
                                     tip.sample.size = n.size.indridae )

drift.test.indridae <- TreeDriftTest(tree = pruned.tree.indridae, mean.list = ed.means.indridae , cov.matrix.list = cov.indridae, sample.sizes = n.size.indridae)
PlotTreeDriftTest(test.list = drift.test.indridae, tree = pruned.tree.indridae)
nodelabels(bg = NULL)
Ancestral.Matrices.indridroidea<- PhyloW(tree = pruned.tree.indroidea, 
                                     tip.data = cov.indroidea, 
                                     tip.sample.size = n.size.indroidea )

RandomSkewers(cov.x = Ancestral.Matrices.indridae$"15",cov.y = Ancestral.Matrices.indridroidea$"27", num.vectors = 1000 )


par(mfrow=c(2,1))