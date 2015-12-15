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
cov.mx.all <- indroidea.data[!na.na.na] %>% llply(function (x) cov(x$ed, use =))
ed.means.all <- indroidea.data[!na.na.na] %>% llply(function (x) Matrix::colMeans(x$ed, na.rm = T) ) 
sample.size.all <- indroidea.data[!na.na.na] %>% ldply(function(x) x$sample.size) %>% .[,2]
names.all <- names(cov.mx.all)
pruned.tree.all<- drop.tip(treefile,treefile$tip.label[-match(names.all, treefile$tip.label)])
plot(pruned.tree.all, cex = 1, adj =0.1, font =3)

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
nodelabels(bg = NULL, frame = "none", adj = -0.4 )

########## comparando os resultados apenas com matrizes bem estimadas
cov.indridae <- cov.mx.sp[indridae]
mask.cov.na.indridae <- ldply(cov.indridae, function(x) !is.na(x[1]))[,2]
names(cov.indridae)[mask.cov.na.indridae == F]
n.size.indridae <- n.size.sp[indridae,2]

RandomSkewers(cov.x = list("Propithecus Genus" = cov.mx.genus$Propithecus,
                           "P.edwardsi" = cov.mx.all$Propithecus_edwardsi,
                           "P.diadema" = cov.mx.all$Propithecus_diadema,
                           "P.tattersalli" = cov.mx.all$Propithecus_tattersalli,
                           "P.coquereli" = cov.mx.all$Propithecus_coquereli,
                           "P.verreauxi" = cov.mx.all$Propithecus_verreauxi,
                           "P.candidus" = cov.mx.all$Propithecus_candidus,
                           "Ancestral P.19" = Ancestral.indroidea.with.mx$"19",
                           "Ancestral P.20" = Ancestral.indroidea.with.mx$"20",
                           "Ancestral P.21" = Ancestral.indroidea.with.mx$"21",
                           "Ancestral P.22" = Ancestral.indroidea.with.mx$"22",
                           "Ancestral P.23" = Ancestral.indroidea.with.mx$"23",
                           "Ancestral P.24" = Ancestral.indroidea.with.mx$"24"),
              num.vectors = 1000)$correlations

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

drift.test.indridae <- TreeDriftTest(tree = tree.indridae, mean.list = ed.means.indridae , cov.matrix.list = cov.indridae, sample.sizes = n.size.indridae)
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

RandomSkewers(cov.x = Ancestral.Matrices.indridae$"15",cov.y = Ancestral.Matrices.indridroidea$"27", num.vectors = 1000 )
par(mfrow=c(2,1))




RandomSkewers(Ancestral.indroidea.with.mx$"14", Ancestral.Matrices.indridae$"14", num.vectors = 1000)
