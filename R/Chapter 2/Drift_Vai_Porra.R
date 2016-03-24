#Drift test dos dois jeitos, com 3 conjudos de dados  #vai porra!

#pegando só os viventes:
cov.mx<- sp.main.data %>% llply(function(x) x$matrix$cov) # todas as matrizes de cov
mask.extant <- ldply(sp.main.data, function(x) unique(x$info$Status) == "Extant")[,2] #mascara para os vventes
mask.is.na.cov <- ldply(cov.mx, function(x) is.na(x[1]))[,2] # mascara dos sem matriz
mask.no.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2] #mascara dos com matriz, mesmo que seja tosca
n.size <- sp.main.data %>% ldply(function(x) x$sample.size) %>% .[,2] #tamanho amostral de todos
mask.n.size <- table(n.size[,2] > 30) # mascara de tamanho amostral >30
ed.means <- sp.main.data %>% llply(function(x) x$ed.means) #médias de todos

#Matrizes de genero a serem atribuidas por espécie
names(cov.mx)[mask.is.na.cov]

#treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
treefile = read.nexus(file = "attaches/fbd421agerange_edited.tre")
plot(treefile, cex = 0.5)
all.my.species <- treefile$tip.label[treefile$tip.label %in% names(ed.means)]
species.extants <- treefile$tip.label[treefile$tip.label %in% names(ed.means[mask.extant])]
species.with.mx <- treefile$tip.label[treefile$tip.label %in% names(ed.means[mask.no.na.cov])]


Trees <- vector("list", 3)

Trees$all.my.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(all.my.species, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$all.my.sp.tree, no.margin = T)
Trees$extant.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(species.extants, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$extant.sp.tree, no.margin = T)
Trees$with.mx.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(species.with.mx, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$with.mx.sp.tree, no.margin = T)


W.geral <- vector("list", 3)

Ancestral.Matrices <- PhyloW(tree = pruned.tree, tip.data = cov.mx[mask.no.na.cov], tip.sample.size = n.size[mask.no.na.cov,2])


#Passo 1: verificar quais sao as especies que tem tamanho amostral muito pequeno para ter matriz
trimed.sp.main.data <- sp.main.data
trimed.sp.main.data$Microcebus_berthae$matrix$cov <- gen.main.data$Microcebus$matrix$cov
trimed.sp.main.data$Microcebus_myoxinus$matrix$cov <- gen.main.data$Microcebus$matrix$cov
trimed.sp.main.data$Microcebus_rufus$matrix$cov <- gen.main.data$Microcebus$matrix$cov
trimed.sp.main.data$Microcebus_tanosi$matrix$cov <- gen.main.data$Microcebus$matrix$cov

trimed.sp.main.data$Cheirogaleus_crossleyi$matrix$cov <- gen.main.data$Cheirogaleus$matrix$cov

trimed.sp.main.data$Phaner_pallescens$matrix$cov <- gen.main.data$Phaner$matrix$cov

trimed.sp.main.data$Lepilemur_tymerlachsonorum$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_dorsalis$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_ahmansonorum$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_edwardsi$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_randrianasoloi$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_petteri$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_betsileo$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_jamesorum$matrix$cov <- gen.main.data$Lepilemur$matrix$cov
trimed.sp.main.data$Lepilemur_microdon$matrix$cov <- gen.main.data$Lepilemur$matrix$cov


trimed.sp.main.data$Avahi_meridionalis$matrix$cov <- gen.main.data$Avahi$matrix$cov
trimed.sp.main.data$Avahi_occidentalis$matrix$cov <- gen.main.data$Avahi$matrix$cov

trimed.sp.main.data$Propithecus_candidus$matrix$cov <- gen.main.data$Propithecus$matrix$cov
trimed.sp.main.data$Propithecus_coquereli$matrix$cov <- gen.main.data$Propithecus$matrix$cov
trimed.sp.main.data$Propithecus_tattersalli$matrix$cov <- gen.main.data$Propithecus$matrix$cov


trimed.sp.main.data$Eulemur_sanfordi$matrix$cov <- gen.main.data$Eulemur$matrix$cov
trimed.sp.main.data$Eulemur_flavifrons$matrix$cov <- gen.main.data$Eulemur$matrix$cov

trimed.sp.main.data$Hapalemur_aureus$matrix$cov <- gen.main.data$Hapalemur$matrix$cov
trimed.sp.main.data$Hapalemur_occidentalis$matrix$cov <- gen.main.data$Hapalemur$matrix$cov
trimed.sp.main.data$Hapalemur_alaotrensis$matrix$cov <- gen.main.data$Hapalemur$matrix$cov


trimed.sp.main.data$Perodicticus_edwardsi$matrix$cov <- gen.main.data$Perodicticus$matrix$cov
trimed.sp.main.data$Perodicticus_ibeanus$matrix$cov <- gen.main.data$Perodicticus$matrix$cov

trimed.sp.main.data$Allocebus_trichotis$matrix$cov <- ((gen.main.data$Mirza$matrix$cov/gen.main.data$Mirza$sample.size) + (gen.main.data$Microcebus$matrix$cov/gen.main.data$Microcebus$sample.size) )

trimed.sp.main.data$Tarsius_bancanus<- gen.main.data$Tarsius
trimed.sp.main.data$Loris_tardigradus <- gen.main.data$Loris
trimed.sp.main.data$Galago_senegalensis <- gen.main.data$Galago
trimed.sp.main.data$Nycticebus_coucang <- gen.main.data$Nycticebus
trimed.sp.main.data$Arctocebus_calabarensis$matrix$cov <- gen.main.data$Perodicticus$matrix$cov
trimed.sp.main.data$Otolemur_crassicaudatus <- gen.main.data$Otolemur
trimed.sp.main.data$Perodicticus_potto<- gen.main.data$Perodicticus
trimed.sp.main.data$Phaner_furcifer<- gen.main.data$Phaner

trimed.main.data<- trimed.sp.main.data[mask.extant]
trimed.main.data <- trimed.main.data[-c(2, 15, 18, 64, 65, 66, 68, 69, 70, 72, 73, 76)]


current.data <- trimed.main.data
#Separar:
#média geométrica por espécie:
gm.mean <- current.data  %>% ldply(function(x) x$gm.mean) 
#médias dos caracteres por especie
ed.means <- current.data  %>% llply(function(x) x$ed.means) 
#Tamanho amostral:
n.size <- current.data  %>% ldply(function(x) x$sample.size) 



# PC1 ao PC4 por espécie
PCs1<- current.data %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,1]) 
PCs2<- current.data %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,2]) 
PCs3<- current.data %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,3]) 
PCs4<- current.data %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,4]) 

cov.mx <- current.data %>% llply(function(x) x$matrix$cov)

write.csv(gm.mean, "gm_mean_Strepsirrhini.csv")
write.csv(ed.means, "ed_means_Strepsirrhini.csv")
write.csv(n.size, "n_size_Strepsirrhini.csv")
write.csv(PCs1, "PCs1_Strepsirrhini.csv")
write.csv(PCs2, "PCs2_Strepsirrhini.csv")
write.csv(PCs3, "PCs3_Strepsirrhini.csv")
write.csv(PCs4, "PCs4_Strepsirrhini.csv")


tree.best.edited = read.nexus(file = "~/ataches/fbd421agerange_edited.tre")
pruned.tree.edited<-drop.tip(tree.best.edited,tree.best.edited$tip.label[-match(names(trimed.main.data), tree.best.edited$tip.label)])
pruned.tree.edited$tip.label
writeNexus(pruned.tree.edited, "Strep_Tree_trim.nex")



tree.drift.test<- TreeDriftTest(tree = pruned.tree.edited, mean.list = ed.means , cov.matrix.list = cov.mx, sample.sizes = n.size[,2])
results <- llply(tree.drift.test, function(x) x$drift_rejected)
PlotTreeDriftTest(test.list = tree.drift.test, tree = pruned.tree.edited)
nodelabels()
tree.drift.test$`45`

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
