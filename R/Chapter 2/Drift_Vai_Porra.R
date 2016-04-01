#Drift test dos dois jeitos, com 3 conjudos de dados  #vai porra!

cov.mx <- sp.main.data %>% llply(function(x) x$matrix$cov) # todas as matrizes de cov
mask.extant <- sp.main.data %>% ldply( function(x) unique(x$info$Status) == "Extant") %>% .[,2] #mascara para os vventes
mask.is.na.cov <- cov.mx %>% ldply( function(x) is.na(x[1])) %>% .[,2] # mascara dos sem matriz
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
names(cov.mx[mask.is.na.cov])

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
Trees$File <-  quote(read.nexus(file = "attaches/fbd421agerange_edited.tre"))
Trees$TREEFILE <- TREEFILE

# calculando as matrizes ancestrais usando as espécies que tenho matriz
ancestral.mx <- PhyloW(tree = pruned.tree.with.mx, tip.data = mx.at.tree, tip.sample.size = n.size[mask][-41])
plot.phylo(Trees$with.mx.sp.tree, no.margin = T) # os numeros dos nós desta arvore é que são os nomes das matrizes ancestrais calculadas desse jeito
nodelabels()
negative.eigen <- ancestral.mx %>% llply(function (x) eigen(x)$values <=0) %>% llply(table) %>% ldply(dim)
W.eigen<- ancestral.mx %>% llply(function (x) eigen(x)$values[eigen(x)$values <=0]) 

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

Drift.results <- vector("list")

Drift.results$with.mx <- TreeDriftTestAll (tree = pruned.tree.with.mx, 
                                           mean.list = ed.means[mask][-41], 
                                           cov.matrix.list = cov.list[-41], 
                                           sample.sizes = sample.size[-c(41, 43, 44)])

mask.at.tree <- names(sp.main.data) %in% Trees$all.my.sp.tree$tip.label
Trees$extant.sp.tree <- drop.tip(treefile,treefile$tip.label[-match(Species$species.extants, treefile$tip.label)]) # árvore com todo mundo
plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.6)
nodelabels(cex = 0.4)
Drift.results$extant.sp <- TreeDriftTestAll (tree = Trees$extant.sp.tree  , 
                                          mean.list = All.sp.data$means[mask.extant & mask.at.tree], 
                                          cov.matrix.list = All.sp.data$cov.mx[mask.extant & mask.at.tree], 
                                          sample.sizes = All.sp.data$n.sizes[mask.extant & mask.at.tree])

Drift.results$extant.sp$Correlation.test.Regular$`71`$P.value.plot
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

par(mfrow= c(1,1))
plot.phylo(Trees$all.with.ed, no.margin = T, cex = 0.5)
nodelabels()
Drift.results$all.sp <- TreeDriftTestAll (tree = Trees$all.with.ed, 
                                          mean.list = All.sp.data$means[mask.at.tree][rowSums(missing.ed) == 39], 
                                          cov.matrix.list = All.sp.data$cov.mx[mask.at.tree][rowSums(missing.ed) == 39], 
                                          sample.sizes = All.sp.data$n.sizes[mask.at.tree][rowSums(missing.ed) == 39] )

Drift.results$with.mx$BW.compare %>% llply(function (x) x$BW.compare) 
Drift.results$all.sp$BW.compare %>% llply(function (x) x$BW.compare) 
Drift.results$extant.sp$BW.compare %>% llply(function (x) x$BW.compare) 



Drift.results.Toplot <- vector("list")
Drift.results.Toplot$All.sp <- vector("list")
Drift.results.Toplot$All.sp$Results <- vector("list")
Drift.results.Toplot$Extants <- vector("list")
Drift.results.Toplot$Extants$Results <- vector("list")

# graficos dos valores de correlaçao
Drift.results.Toplot$All.sp$Plots$Corr.Contrasts  <- Drift.results$all.sp$Correlation.test.Contrasts %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none") ) %>% cowplot::plot_grid(plotlist = .) 
Drift.results.Toplot$All.sp$Plots$Corr.Ed         <- Drift.results$all.sp$Correlation.test.Regular %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none")) %>% cowplot::plot_grid(plotlist = .)
Drift.results.Toplot$Extants$Plots$Corr.Contrasts <- Drift.results$extant.sp$Correlation.test.Contrasts %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none")) %>% cowplot::plot_grid(plotlist = .)
Drift.results.Toplot$Extants$Plots$Corr.Ed        <- Drift.results$extant.sp$Correlation.test.Regular %>% llply(function (x) x$P.value.plot) %>% llply(function(x) x + theme(legend.position = "none")) %>% cowplot::plot_grid(plotlist = .)

# Pegando em quais testes de correlaçao apareceu uma correlaçao significativa no teste, corrigindo por bonferroni
# se tiver alguma correlação que ainda com a correçao de vonferroni é significativa, entao rejeitamos deriva.
# nessa tabela a pergunta era se p.value < 0.05 / fator de correçao. (se sim seria deriva)
# quando o valor da correlaçao é menor que isso ele plota FALSE no par de correlaçao
# entao se tiver FALSE é pq rejeitamos deriva.
Drift.results$all.sp$Correlation.test.Regular$`81`$Bonferroni 
Drift.results$all.sp$Correlation.test.Regular$`81`$P.value.plot
# para facilitar vou mandar um table em quais dessas correlaçoes deu como TRUE, mandar ele excluir quem tem TRUE.
# nesse caso sobrarão apenas os que tem NA e aqueles casos em que rejeitamos Deriva. Nesses casos o resultado do table terá 2 elementos.
# se eu perguntar quem é ==2 pego aqueles em que rejeitamos deriva
Drift.results.Toplot$All.sp$Results$Node.ref <- Drift.results$all.sp$Correlation.test.Regular %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,1]
Drift.results.Toplot$All.sp$Results$Corr.Ed.1 <- Drift.results$all.sp$Correlation.test.Regular %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,2]
Drift.results.Toplot$All.sp$Results$Corr.Contrasts.1 <- Drift.results$all.sp$Correlation.test.Contrasts %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) !=1) %>% .[,2]
Drift.results.Toplot$Extants$Results$Node.ref <- Drift.results$extant.sp$Correlation.test.Regular %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,1]
Drift.results.Toplot$Extants$Results$Corr.Ed.1 <- Drift.results$extant.sp$Correlation.test.Regular %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,2]
Drift.results.Toplot$Extants$Results$Corr.Contrasts.1 <- Drift.results$extant.sp$Correlation.test.Contrasts %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,2]

# como pegar o resultado de deriva no teste de regressao:
# 1 estiver contido no Intervalo de confiança da regressao : Deriva nao é rejeitada
# quando $drift_rejected = TRUE : rejeitamos a H0 de Deriva
Drift.results$all.sp$Regression.test$`81`$coefficient_CI_95
Drift.results$all.sp$Regression.test$`81`$drift_rejected #TRUE
Drift.results$extant.sp$Regression.test$`71`$coefficient_CI_95
Drift.results$extant.sp$Regression.test$`71`$drift_rejected #TRUE
Drift.results$with.mx$Regression.test$`42`$coefficient_CI_95
Drift.results$extant.sp$Regression.test$`71`$drift_rejected #TRUE

Drift.results.Toplot$All.sp$Results$Regression.Ed <- Drift.results$all.sp$Regression.test %>% ldply(function(x) x$drift_rejected )  %>% .[,2]
Drift.results.Toplot$All.sp$Results$Regression.Contrasts <- Drift.results$all.sp$Regression.test.Contrasts %>% ldply(function(x) x$drift_rejected )  %>% .[,2]
Drift.results.Toplot$Extants$Results$Regression.Ed <- Drift.results$extant.sp$Regression.test %>% ldply(function(x) x$drift_rejected )  %>% .[,2]
Drift.results.Toplot$Extants$Results$Regression.Contrasts <- Drift.results$extant.sp$Regression.test.Contrasts %>% ldply(function(x) x$drift_rejected )  %>% .[,2]

Plot.Drift.Results <- function (tree, node.ref, info, font.size = 0.4) {
  node.ref <- as.numeric(info$Node.ref)
  results.ed <- data.frame(V1 = as.numeric(info$Corr.Ed.1) * 0.5, V2 = as.numeric(info$Regression.Ed) * 0.5  )
  results.ci <- data.frame(V1 = as.numeric(info$Corr.Contrasts.1) * 0.5, V2 = as.numeric(info$Regression.Contrasts) * 0.5) 
  
  results.ed$V3 <- abs(results.ed$V1 + results.ed$V2 -1)
  results.ci$V3 <- abs(results.ci$V1 + results.ci$V2 -1)
  results.ed<- as.matrix(results.ed)
  results.ci<- as.matrix(results.ci)
  par(mfrow = c(1,2))
  plot.phylo(tree, font = 3, no.margin = TRUE,  cex = font.size, edge.color = "darkgrey", edge.width = 3, label.offset = 1 )
  nodelabels(node = node.ref, 
             pie = results.ed[,1:3], cex=0.6, 
             piecol = c("#FFAF02", "#0912C2", "#C3BBBF")  )
  #title(main = "B.matrices obtained from ED")
  #nodelabels(node = node.ref, bg = "transparent", frame = "n",cex = 0.6, col = "black", adj = 0.9)
  legend("bottomleft", inset = .08,
         title ="Drift test \nB.matrices obtained from ED",
         text.col = "grey10",
         legend = c("Rejected - Correlation", "Rejected - Regression", "Not rejected"), 
         fill = c("#FFAF02", "#0912C2", "#C3BBBF"), 
         col = c("#FFAF02", "#0912C2", "#C3BBBF"), 
         border = "grey", box.lwd = "n",
         bg= "transparent",
         cex = 0.5)
  plot.phylo(tree, font = 3, no.margin = TRUE,  cex = font.size, edge.color = "darkgrey", edge.width = 3, label.offset = 1 )
  nodelabels(node = node.ref, 
             pie = results.ci[,1:3], cex=0.6, 
             piecol = c("#FFAF02", "#0912C2", "#C3BBBF")  )
  #title(main = "B.matrices obtained from IC")
  #nodelabels(node = node.ref, bg = "transparent", frame = "n",cex = 0.6, col = "black", adj = 0.9)
  legend("bottomleft", inset = 0.08,
         title ="Drift test \nB.matrices obtained from IC",
         text.col = "grey10",
         legend = c("Rejected - Correlation", "Rejected - Regression", "Not rejected"), 
         fill = c("#FFAF02", "#0912C2", "#C3BBBF"), 
         col = c("#FFAF02", "#0912C2", "#C3BBBF"), 
         border = "grey", box.lwd = "n",
         bg= "transparent",
         cex = 0.5)
 
    par(mfrow = c(1,1))
    return(data.frame("node" = node.ref,
           "ed" = results.ed,
           "ic" = results.ci) )
}
Drift.results.Toplot$Extants$Plots$Tree <- Plot.Drift.Results(tree = Trees$extant.sp.tree, info = Drift.results.Toplot$Extants$Results, font.size = 0.4)
Drift.results.Toplot$All.sp$Plots$Tree <- Plot.Drift.Results(tree = Trees$all.with.ed, info = Drift.results.Toplot$All.sp$Results, font.size = 0.4)


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
