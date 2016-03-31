#Drift test dos dois jeitos, com 3 conjudos de dados  #vai porra!

cov.mx <- sp.main.data %>% llply(function(x) x) # todas as matrizes de cov
mask.extant <- ldply(sp.main.data, function(x) unique(x$info$Status) == "Extant")[,2] #mascara para os vventes
mask.is.na.cov <- ldply(cov.mx, function(x) is.na(x[1]))[,2] # mascara dos sem matriz
mask.no.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2] #mascara dos com matriz, mesmo que seja tosca
n.size <- sp.main.data %>% ldply(function(x) x$sample.size) %>% .[,2] #tamanho amostral de todos
mask.n.size <- n.size> 30 # mascara de tamanho amostral >30
ed.means <- sp.main.data %>% llply(function(x) x$ed.means) #médias de todos

B.var <- sp.main.data[mask.no.na.cov] %>% ldply(function(x) x$ed.means) %>% .[, -1] %>% var
B.all <- sp.main.data[mask.no.na.cov] %>% ldply(function(x) x$ed)
B.all.lm <- lm(as.matrix(B.all[,-1]) ~ B.all[, 1])
B.all.lm.1 <- lm(as.matrix(B.all[,-1]) ~ 1)

sumsqr.W <- t(B.all.lm$residuals) %*% B.all.lm$residuals 
W.lm <- CalculateMatrix(B.all.lm)
W.lm - Ancestral.Matrices$`42`
sumsqr.W.t <- t(B.all.lm.1$residuals) %*% B.all.lm.1$residuals 

B.sumsqr <- (sumsqr.W.t - sumsqr.W)

MatrixCompare(B.sumsqr, Ancestral.Matrices$`45`)

#Matrizes de genero a serem atribuidas por espécie
names(cov.mx)[mask.is.na.cov]

#treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
TREEFILE = read.nexus(file = "attaches/fbd421agerange_edited.tre")
plot(TREEFILE, cex = 0.5)
nodelabels()
Species <- vector("list", 3)
Species$all.my.species <- treefile$tip.label[treefile$tip.label %in% names(ed.means)]
Species$species.extants <- treefile$tip.label[treefile$tip.label %in% names(ed.means[mask.extant])]
Species$species.with.mx <- treefile$tip.label[treefile$tip.label %in% names(ed.means[mask.no.na.cov])]

Trees <- vector("list")
Trees$all.my.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(Species$all.my.species, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$all.my.sp.tree, no.margin = T, cex = 0.5)
nodelabels()
Trees$extant.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(Species$species.extants, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$extant.sp.tree, no.margin = T)
nodelabels()
Trees$with.mx.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(Species$species.with.mx, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$with.mx.sp.tree, no.margin = T)
nodelabels()

ancestral.mx <- PhyloW(tree = pruned.tree.with.mx, tip.data = mx.at.tree, tip.sample.size = n.size[mask][-41])
plot.phylo(Trees$with.mx.sp.tree, no.margin = T) # os numeros dos nós desta arvore é que são os nomes das matrizes ancestrais calculadas desse jeito
nodelabels()

#Passo 1: verificar quais sao as especies que tem tamanho amostral muito pequeno para ter matriz
All.sp.data <- vector("list")

All.sp.data$cov.mx <- sp.main.data %>% llply(function(x) x$matrix$cov) # todas as matrizes de cov
All.sp.data$ed <- sp.main.data %>% llply(function(x) x$ed) # todas as distancias médias de cada indivíduo
All.sp.data$means <- sp.main.data %>% llply(function(x) colMeans(x$ed, na.rm = T) ) #todos os vetores de médias
All.sp.data$gm <- sp.main.data %>% llply(function(x) x$gm.mean) #todos os vetores de médias
All.sp.data$n.sizes <- sp.main.data %>% ldply(function(x) x$sample.size )  %>% .[,2] 

# Cheirogaleidae family ####
All.sp.data$cov.mx$Allocebus_trichotis <- ancestral.mx$`50` # ancestral de microcebus e mirza

All.sp.data$cov.mx$Microcebus_berthae <- Gen.cov.list$Microcebus
All.sp.data$cov.mx$Microcebus_myoxinus <- Gen.cov.list$Microcebus
All.sp.data$cov.mx$Microcebus_rufus <- Gen.cov.list$Microcebus
All.sp.data$cov.mx$Microcebus_tanosi <- Gen.cov.list$Microcebus

All.sp.data$cov.mx$Cheirogaleus_crossleyi <- Gen.cov.list$Cheirogaleus

All.sp.data$cov.mx$Phaner_pallescens <- Gen.cov.list$Phaner
All.sp.data$cov.mx$Phaner_furcifer<- Gen.cov.list$Phaner


# Lepilemuridae family ####
All.sp.data$cov.mx$Lepilemur_tymerlachsonorum <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_dorsalis <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_ahmansonorum <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_edwardsi <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_randrianasoloi <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_petteri <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_betsileo <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_jamesorum <- Gen.cov.list$Lepilemur
All.sp.data$cov.mx$Lepilemur_microdon <- Gen.cov.list$Lepilemur

# Indridae family ####
All.sp.data$cov.mx$Avahi_meridionalis <- Gen.cov.list$Avahi
All.sp.data$cov.mx$Avahi_occidentalis <- Gen.cov.list$Avahi

All.sp.data$cov.mx$Propithecus_candidus <- Gen.cov.list$Propithecus
All.sp.data$cov.mx$Propithecus_coquereli <- Gen.cov.list$Propithecus
All.sp.data$cov.mx$Propithecus_tattersalli <- Gen.cov.list$Propithecus

# Lemuridae family ####
All.sp.data$cov.mx$Eulemur_sanfordi <- Gen.cov.list$Eulemur
All.sp.data$cov.mx$Eulemur_flavifrons <- Gen.cov.list$Eulemur

All.sp.data$cov.mx$Hapalemur_aureus <- Gen.cov.list$Hapalemur
All.sp.data$cov.mx$Hapalemur_occidentalis <- Gen.cov.list$Hapalemur
All.sp.data$cov.mx$Hapalemur_alaotrensis <- Gen.cov.list$Hapalemur

All.sp.data$cov.mx$Perodicticus_edwardsi <- Gen.cov.list$Perodicticus
All.sp.data$cov.mx$Perodicticus_ibeanus <- Gen.cov.list$Perodicticus


All.sp.data$cov.mx$Tarsius_bancanus<- Gen.cov.list$Tarsius
All.sp.data$cov.mx$Loris_tardigradus <- Gen.cov.list$Loris
All.sp.data$cov.mx$Galago_senegalensis <- Gen.cov.list$Galago
All.sp.data$cov.mx$Nycticebus_bengalensis <- Gen.cov.list$Nycticebus
All.sp.data$cov.mx$Nycticebus_menagensis<- Gen.cov.list$Nycticebus
All.sp.data$cov.mx$Nycticebus_javanicus<- Gen.cov.list$Nycticebus
All.sp.data$cov.mx$Nycticebus_coucang <- Gen.cov.list$Nycticebus
All.sp.data$cov.mx$Arctocebus_calabarensis <- Gen.cov.list$Perodicticus
All.sp.data$cov.mx$Otolemur_crassicaudatus <- Gen.cov.list$Otolemur
All.sp.data$cov.mx$Perodicticus_potto<- Gen.cov.list$Perodicticus

 # Extintos ####
All.sp.data$cov.mx$Hadropithecus_stenognathus <- ancestral.mx$`55` # recebe o ancestral da família  Indroideia sem Indri
All.sp.data$cov.mx$Archaeolemur_majori <- ancestral.mx$`55` # recebe o ancestral da família  Indroideia sem Indri
All.sp.data$cov.mx$Archaeolemur_edwardsi <- ancestral.mx$`55` # recebe o ancestral da família  Indroideia sem Indri
All.sp.data$cov.mx$Mesopropithecus_pithecoides <- Gen.cov.list$W.Indrida
All.sp.data$cov.mx$Mesopropithecus_dolichobrachion <- Gen.cov.list$W.Indrida
All.sp.data$cov.mx$Palaeopropithecus_ingens <- Gen.cov.list$W.Indrida
All.sp.data$cov.mx$Palaeopropithecus_maximus <- Gen.cov.list$W.Indrida
All.sp.data$cov.mx$Babakotia_radafolia <- Gen.cov.list$W.Indrida

All.sp.data$cov.mx$Pachylemur_jullyi <- Gen.cov.list$Varecia


All.sp.data$cov.mx$Megaladapis_edwardsi <- ancestral.mx$`46` # recebe o ancestral de madagascar sem Daubentonia
All.sp.data$cov.mx$Megaladapis_madagascariensis <- ancestral.mx$`46` # recebe o ancestral de madagascar sem Daubentonia
All.sp.data$cov.mx$Megaladapis_edwardsi <- ancestral.mx$`46` # recebe o ancestral de madagascar sem Daubentonia
All.sp.data$cov.mx$Megaladapis_grandidieri <- ancestral.mx$`46`
All.sp.data$cov.mx[mask.at.tree] %>% ldply( function(x) is.na(x)[1]) %>% .[,2] %>% table # conferindo se todos receberam uma matriz

mask.at.tree <- names(sp.main.data) %in% Trees$all.my.sp.tree$tip.label
Trees$extant.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(Species$species.extants, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.6)
nodelabels(cex = 0.4)
Drift.results$extant.sp <- TreeDriftTestAll (tree = Trees$extant.sp.tree  , 
                                          mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                                          cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                                          sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree])

# alguns dos fósseis nao tem informação de algumas medidas: quem sao eles?
All.sp.data$means[mask.at.tree]%>% llply( function(x) !is.na(x) ) %>% ldply( function(x) sum(x) <39 ) # conferindo se todos receberam uma matriz
missing.ed <- All.sp.data$means[mask.at.tree]%>% ldply( function(x) !is.na(x) ) 
row.names(missing.ed) <- missing.ed$.id
missing.ed <- missing.ed[,-1]
table(colSums(missing.ed) == 82 )
table(rowSums(missing.ed) )
which(rowSums(missing.ed) < 39)


Trees$all.my.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(Species$all.my.species, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$all.my.sp.tree, no.margin = T, cex = 0.5)
nodelabels()

Trees$all.with.ed<- drop.tip(treefile,treefile$tip.label[-match(names(All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39]), treefile$tip.label)])

table(names(tree.no.megaladapis$tip.label %in% All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39]) )

Drift.results$all.sp <- TreeDriftTestAll (tree = Trees$all.with.ed, 
                                          mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                                          cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                                          sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )
drift.vai.porra <- vector("list", 5)

drift.vai.porra$drift.T.F.vectors <- 
  cbind(Drift.alltests.tree$Correlation.test.Regular %>% ldply(function(x) dim(x$Resume.table)[[1]] >1 ),
        Drift.alltests.tree$Correlation.test.Contrasts %>% ldply(function(x) dim(x$Resume.table)[[1]] >1 ) %>% .[,2],
        Drift.alltests.tree$Regression.test %>% ldply(function(x) x$drift_rejected ) %>% .[,2],
        Drift.alltests.tree$Regression.test.Contrasts %>% ldply(function(x) x$drift_rejected )  %>% .[,2])

colnames(drift.vai.porra$drift.T.F.vectors) <- c("node", "cor", "cor.ci", "reg", "reg.ci")
str(drift.vai.porra$drift.T.F.vectors)
drift.vai.porra$drift.T.F.vectors$node <- as.numeric(drift.vai.porra$drift.T.F.vectors$node)

par(mfrow = c(1,2))
plot.phylo(pruned.tree.with.mx, font = 3, no.margin = T)
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 8, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$cor)+1), frame = "n")
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 17, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$reg)+3), frame = "n")

plot.phylo(pruned.tree.with.mx, font = 3, no.margin = T)
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 8, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$cor.ci)+1), frame = "n")
nodelabels(node = drift.vai.porra$drift.T.F.vectors$node , pch = 17, bg = "transparent", col = (as.numeric(drift.vai.porra$drift.T.F.vectors$reg.ci)+3), frame = "n")

drift.vai.porra$sum.abs.values <- Drift.alltests.tree$Correlation.test.Regular %>% llply(function (x) x$Correlation.p.value[1:39,1:39]) %>% llply(abs) %>% laply( function (x) x[lower.tri(x)]) %>% colSums 
temp <- matrix(NA, 39, 39, byrow = F)
temp[lower.tri(temp)] <- drift.vai.porra$sum.abs.values
drift.vai.porra$mean.abs.values <- temp / length(Drift.alltests.tree$Correlation.test.Regular)



######### changing plot 
tested.nodes <- as.numeric(names(tree.drift.test))
non.drift.nodes <- laply(tree.drift.test, function(x) x$drift_rejected)
i.c.5 <- ldply(tree.drift.test, function(x) x$coefficient_CI_95[2,])
i.c.5<- ldply(tree.drift.test, function(x) x$coefficient_CI_95[1,])
names(i.c.5) <- c('node', 'min', 'max')
i.c.5$min <- round(i.c.5$min, digits=2)
i.c.5$max <- round(i.c.5$max, digits=2)

# par(mfrow = c(1,1))
par(mar=c(1,1,1,1))

plot(pruned.tree.edited, font = 3, label.offset = 11, no.margin = T, edge.color = "grey80", edge.width = 5, cex= 0.9, x.lim = c(-8, 100) )

#nodelabels(node = tested.nodes, thermo = as.numeric(non.drift.nodes))
i.c.5 <- ldply(tree.drift.test, function(x) x$coefficient_CI_95[2,])
names(i.c.5) <- c('node', 'min', 'max')
i.c.5$min <- round(i.c.5$min, digits=2)
i.c.5$max <- round(i.c.5$max, digits=2)
nodelabels(node = tested.nodes , text = i.c.5$min, adj = 1.2, bg = "transparent", col = 'red', frame = "n", cex = i.c.5$min)
nodelabels(node = tested.nodes, i.c.5$max, adj = -0.2, bg = "transparent", col = 'blue', frame = "n", cex = i.c.5$max)
nodelabels(node = tested.nodes, pch = (as.numeric(non.drift.nodes)+17), cex=1.4, col = (as.numeric(non.drift.nodes)+5) )

compare.stuff<- cbind( n.size[n.size[,1] %in% pruned.tree.edited$tip.label, ], gm.mean[gm.mean[,1] %in% pruned.tree.edited$tip.label, 2])
names(compare.stuff) <- c('Especie', 'n', 'gm')
compare.stuff<- compare.stuff[with(compare.stuff, order(Especie)), ]
rownames(compare.stuff)<- compare.stuff$Especie
compare.stuff<- picante::match.phylo.data(phy = pruned.tree.edited, data = compare.stuff)$data
compare.stuff$n <- as.numeric(as.character(compare.stuff$n))
compare.stuff$gm <- as.numeric(as.character(compare.stuff$gm))
tiplabels(pch = 21, cex =compare.stuff$n/18, adj = 6)
tiplabels(pch = 15, cex = compare.stuff$gm/9, adj = 9)
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

PCScoreCorrelation(means = ed.means, cov.matrix = cov.mx, show.plots = T)
