 ###################### BRINCANDO COM OS FOSSEIS E TESTE DE DERIVA ############################
#olhando pro pedaço da árvore que quero trabalhar:
#arvore trimada para Indroidea: diversidade total desse clado
treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
plot(treefile, cex = 0.5)
nodelabels(bg = NULL, frame = "none", node = 222, pch = 21, col = "red", text = "Indroidea" )
tree.indroidea.total <- extract.clade(treefile, 222)
plot(tree.indroidea.total, cex = 1, adj = 0.1, font =3 )

#Filtrando matrizes, médias e tamanhos amostrais pra este clado
ed.means <- sp.main.data %>% llply(function(x) x$ed.means) 
cov.mx.sp <- sp.main.data %>% llply(function(x) x$matrix$cov)
cov.mx.genus <- gen.main.data %>% llply(function(x) x$matrix$cov)
n.size.sp <- sp.main.data %>% ldply(function(x) x$sample.size) 

#Olhando quem são os fosseis da minha amostra
extintos <- sp.main.data %>% laply(function(x) all(x$info$Status == "Extinct") )
extintos.names <- names(ed.means)[extintos]

#separando os fosseis e viventes dessa familia indroidea apenas e os que tenho em minha amostra
indroidea <- names(ed.means) %in% tree.indroidea.total$tip.label
indroidea.names <- names(ed.means)[indroidea]
#Pegando os caras da familia indridae e os extintos, que todos os que eu tenho tao dentro dela
indridae <- grepl("Indri", names(ed.means)) | grepl("Avahi", names(ed.means)) | grepl("Propithecus", names(ed.means))
indridae.names <- names(ed.means)[indridae]
#Arvore dos indroidea: viventes e extintos: apenas os que tenho em minha amostra
tree.indroidea<- drop.tip(treefile,treefile$tip.label[-match(indroidea.names, treefile$tip.label)])
plot(tree.indroidea, cex = 1, adj =0.1, font =3 )
tree.indroidea$tip.label %in% extintos.names
#Árvore dos indridade: apenas viventes que tenho na minha amostra
tree.indridae<- drop.tip(treefile,treefile$tip.label[-match(indridae.names, treefile$tip.label)])
plot(tree.indridae, cex = 1, adj = 0.1, font =3)

###############
# 1. estimar matrizes pra tudo quanto é coisa que tem amostra >2
# 2. Rodar a Pholow(  ):
###### 3. PhyloW(tree = tree.file, tip.data = cov.matrices.low.sample.sizes, tip.sample.size = real.sample.sizes)
# 4. colocar nas matrizes que só tem um bicho a matriz do nó mais próximo.

############################ 1.  estimar matrizes pra tudo quanto é coisa que tem amostra >2
#Tirando quem tem algum NA na amostra
indroidea.data <- sp.main.data[indroidea]
#Mascara para aqueles que nao tem medidas de algumas distancias at all.
na.na.na<- rowSums(is.na(indroidea.data %>% ldply(function (x) Matrix::colMeans(x$ed, na.rm = T) ))) !=0
#selecionando e montando as listas em que vou trabalhar
cov.mx.all <- indroidea.data[!na.na.na] %>% llply(function (x) cov(x$ed))
ed.means.all <- indroidea.data[!na.na.na] %>% llply(function (x) Matrix::colMeans(x$ed, na.rm = T) ) 
sample.size.all <- indroidea.data[!na.na.na] %>% ldply(function(x) x$sample.size) %>% .[,2]
names.all <- names(cov.mx.all)
pruned.tree.all<- drop.tip(treefile,treefile$tip.label[-match(names.all, treefile$tip.label)])
plot(pruned.tree.all, cex = 1, adj =0.1, font =3)

cov.mx.all.backup <- cov.mx.all
# verificando se tem algum aqui que ta como NA nas matrizes: provavelmente sao os caras que tem sample.size == 1
names(cov.mx.all)[sample.size.all == 1]
############################################# Começando a meter-lhes matrizes bem estimadas
#o que fazer agora: meter-lhes matrizes bem estimadas
#Passo 1: rodar o teste usando só quem tem matriz e rodar a PhyloW()

no.mx <- cov.mx.all %>% ldply(function(x) sum(is.na(x) ) ) %>% .[,2] !=0
names.with.mx <- names(cov.mx.all)[!no.mx]

cov.mx.with.mx <- cov.mx.all[!no.mx] 
ed.means.with.mx <- ed.means.all[!no.mx] 
sample.size.with.mx <- sample.size.all[!no.mx]
names.with.mx <- names(cov.mx.with.mx)
#Arvore dos indroidea: viventes e extintos: apenas os que tem matriz completa
tree.with.mx<- drop.tip(treefile,treefile$tip.label[-match(names.with.mx, treefile$tip.label)])
plot(tree.with.mx, cex = 1, adj =0.1, font =3)

Ancestral.indroidea.with.mx <- PhyloW(tree = tree.with.mx, 
                                      tip.data = cov.mx.with.mx, 
                                      tip.sample.size = sample.size.with.mx )

drift.test.indroidea.with.mx <- TreeDriftTest(tree = tree.with.mx, 
                                     mean.list = ed.means.with.mx, 
                                     cov.matrix.list = cov.mx.with.mx, 
                                     sample.sizes = sample.size.with.mx)

PlotTreeDriftTest(test.list = drift.test.indroidea.with.mx, tree = tree.with.mx, cex = 1, adj =0.1, font =3)
nodelabels(bg = NULL, frame = "none", adj = -1 )

#################################### vendo a quem atribuir quais matrizes ancestrais ########################
names(cov.mx.all)[no.mx]
sample.size.all
#cov.mx.all <- indroidea.data[!na.na.na] %>% llply(function (x) cov(x$ed))
#ed.means.all <- indroidea.data[!na.na.na] %>% llply(function (x) Matrix::colMeans(x$ed, na.rm = T) ) 
#sample.size.all <- indroidea.data[!na.na.na] %>% ldply(function(x) x$sample.size) %>% .[,2]
#names.all <- names(cov.mx.all)
par(mfrow=c(2,1))
PlotTreeDriftTest(test.list = drift.test.indroidea.with.mx, tree = tree.with.mx, cex = 1, adj =0.1, font =3)
nodelabels(bg = NULL, frame = "none", adj = -1 )
pruned.tree.all<- drop.tip(treefile,treefile$tip.label[-match(names.all, treefile$tip.label)])
plot(pruned.tree.all, cex = 1, adj =0.1, font =3)
par(mfrow=c(1,1))
##### substituindo nos extintos
cov.mx.all$Hadropithecus_stenognathus <- Ancestral.indroidea.with.mx$`14`
cov.mx.all$Archaeolemur_majori <- Ancestral.indroidea.with.mx$`14`
cov.mx.all$Archaeolemur_edwardsi <- Ancestral.indroidea.with.mx$`14`
cov.mx.all$Mesopropithecus_pithecoides <- Ancestral.indroidea.with.mx$`15`
cov.mx.all$Mesopropithecus_dolichobrachion <- Ancestral.indroidea.with.mx$`15`
cov.mx.all$Palaeopropithecus_ingens <- Ancestral.indroidea.with.mx$`15`
cov.mx.all$Babakotia_radafolia <- Ancestral.indroidea.with.mx$`15`
###### substituindo nos viventes que só tem média
cov.mx.all$Avahi_meridionalis <- Ancestral.indroidea.with.mx$"18"
cov.mx.all$Avahi_occidentalis <- Ancestral.indroidea.with.mx$"17"



######### verificando os novos resultados dos testes dadas as substituiçoes

Ancestral.indroidea.subistitute.mx <- PhyloW(tree = pruned.tree.all, 
                                      tip.data = cov.mx.all, 
                                      tip.sample.size = sample.size.all )

drift.test.indroidea.subistitute.mx <- TreeDriftTest(tree = pruned.tree.all, 
                                              mean.list = ed.means.all, 
                                              cov.matrix.list = cov.mx.all, 
                                              sample.sizes = sample.size.all)
length(pruned.tree.all$tip.label)
length(ed.means.all)
length(cov.mx.all)
length(sample.size.all)

PlotTreeDriftTest(test.list = drift.test.indroidea.subistitute.mx, tree = pruned.tree.all, cex = 1, adj =0.1, font =3)
nodelabels(bg = NULL, frame = "none", adj = -1 )

pruned.tree.all$edge[which(pruned.tree.all$edge[, 1] == "31"), 2]

library(phytools)
all.equal(ed.means.all[na.omit(pruned.tree.all$tip.label[getDescendants(pruned.tree.all, 31)])],
ed.means.with.mx[na.omit(tree.with.mx$tip.label[getDescendants(tree.with.mx, 19)])])


###################### Porque diabos o nó do genero Propithecus estaria mudando?
RandomSkewers(cov.x = Ancestral.Matrices.indridae$"15",cov.y = Ancestral.Matrices.indridroidea$"27", num.vectors = 1000 )
par(mfrow=c(2,1))
plot_grid(labels = c("node 19: drift rejected IC 95%= 1.02  1.36 ", "node 31: drift accepted IC 95%= 0.82 1.26"), drift.test.indroidea.with.mx$`19`$plot, drift.test.indroidea.subistitute.mx$`31`$plot)
Propithecus.list<- list("Propithecus Genus" = cov.mx.genus$Propithecus,
                        "P.edwardsi" = cov.mx.all$Propithecus_edwardsi,
                        "P.diadema" = cov.mx.all$Propithecus_diadema,
                        "P.tattersalli" = cov.mx.all$Propithecus_tattersalli,
                        "P.coquereli" = cov.mx.all$Propithecus_coquereli,
                        "P.verreauxi" = cov.mx.all$Propithecus_verreauxi,
                        "P.deckenii" = cov.mx.all$Propithecus_deckenii,
                        "P.coronatus" = cov.mx.all$Propithecus_coronatus,
                        "P.candidus" = cov.mx.all$Propithecus_candidus,
                        "node 19" = Ancestral.indroidea.with.mx$"19",
                        "node 20" = Ancestral.indroidea.with.mx$"20",
                        "node 21" = Ancestral.indroidea.with.mx$"21",
                        "node 22" = Ancestral.indroidea.with.mx$"22",
                        "node 23" = Ancestral.indroidea.with.mx$"23",
                        "node 24" = Ancestral.indroidea.with.mx$"24",
                        "node 31" = Ancestral.indroidea.subistitute.mx$"31")
Propithecus <- grepl("Propithecus", names(sp.main.data)) 
Sample.Propithecus <- c(201, 28, 43, 7, 14, 58, 21, 17, 13, 201, 130, 180, 167, 109, 27, 201)

RS.Propithecus <- RandomSkewers(cov.x = Propithecus.list, num.vectors = 1000)$correlations
RS.Propithecus <- t(RS.Propithecus)
KRZ.Propithecus <- KrzCor(cov.x = Propithecus.list, ret.dim = 16)
KRZ.Propithecus <- t(KRZ.Propithecus)

plot.propithecus.compare.all.mx<- Combine.Mx.Plot(Mx1 = RS.Propithecus, Mx2 = KRZ.Propithecus, diag.info = Sample.Propithecus, titulo = "Propithecus genus comparison")

#conferindo o nome dos terminais filhos dos respectivos nós de Propithecus.
names(evolqg:::getMeans(ed.means.all, pruned.tree.all, "31"))
names(evolqg:::getMeans(ed.means.with.mx, tree.with.mx, "19"))

# tava errada a evolqg:::getMeans mas agora o Ograo ajeitou, ta ok.

cov.indridae$Propithecus_candidus <- Ancestral.indroidea.with.mx$"21"
cov.indridae$Propithecus_coquereli <- Ancestral.indroidea.with.mx$"24"
cov.indridae$Avahi_meridionalis <- Ancestral.indroidea.with.mx$"18"
cov.indridae$Avahi_occidentalis <- Ancestral.indroidea.with.mx$"17"
cov.indridae$Propithecus_tattersalli <-  Ancestral.indroidea.with.mx$"24"
#cov.indridae$Propithecus_coquereli <- cov.mx.genus$Propithecus
#cov.indridae$Propithecus_candidus <- cov.mx.genus$Propithecus

#n.size.indridae[mask.cov.na.indridae == F] <- 1

Ancestral.Matrices.indridae<- PhyloW(tree = tree.indridae, 
                                     tip.data = cov.indridae, 
                                     tip.sample.size = n.size.indridae )
######## Ancestral de indridae viventes apenas
Ancestral.Matrices.indridae$"14"
# bem parecido com o considerando o unico fossil que da pra ter matriz 9mas tbm, essa porra tem só 2 bichos)
RandomSkewers(Ancestral.indroidea.with.mx$"14", Ancestral.Matrices.indridae$"14", num.vectors = 1000) # = 0.998861314 

drift.test.indridae <- TreeDriftTest(tree = tree.indridae, mean.list = ed.means.indridae , cov.matrix.list = cov.indridae, sample.sizes = n.size.indridae[,2])
PlotTreeDriftTest(test.list = drift.test.indridae, tree = tree.indridae, cex = 1, adj =0.1, font =3)
nodelabels(bg = NULL, frame = "none", adj = -0.4 )

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
PlotTreeDriftTest(test.list = drift.test.indroidea, tree = pruned.tree.indroidea, cex = 1, adj =0.1, font =3)
nodelabels(bg = NULL, frame = "none", adj = -0.4 )

Ancestral.Matrices.indridae<- PhyloW(tree = pruned.tree.indridae, 
                                     tip.data = cov.indridae, 
                                     tip.sample.size = n.size.indridae )

drift.test.indridae <- TreeDriftTest(tree = pruned.tree.indridae, mean.list = ed.means.indridae , cov.matrix.list = cov.indridae, sample.sizes = n.size.indridae)
PlotTreeDriftTest(test.list = drift.test.indridae, tree = tree.indridae, cex = 1, adj =0.1, font =3)
nodelabels(bg = NULL, frame = "none", adj = -0.4 )

Ancestral.Matrices.indridroidea<- PhyloW(tree = pruned.tree.indroidea, 
                                     tip.data = cov.indroidea, 
                                     tip.sample.size = n.size.indroidea )






RandomSkewers(Ancestral.indroidea.with.mx$"14", Ancestral.Matrices.indridae$"14", num.vectors = 1000)
