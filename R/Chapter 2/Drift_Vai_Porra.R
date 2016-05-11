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
plot.phylo(Trees$extant.sp.tree, no.margin = T, cex = 0.8)
nodelabels(cex = 0.9)

All.sp.data$cov.mx$Otolemur_crassicaudatus <- gen.cov.mx$Otolemur

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

Drift.results$all.sp$BW.compare$`81`$W

Drift.results$all.sp$BW.compare %>% ldply(function (x) cbind(Pc1Percent(x$B.ed ), Pc1Percent(x$W ) ) )
Drift.results$extant.sp$BW.compare %>% ldply(function (x) cbind(Pc1Percent(x$B.ed ), Pc1Percent(x$W ) ) )

  
  
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
Drift.results.Toplot$All.sp$Results$Corr.Ed.1 <- Drift.results$all.sp$Correlation.test.Regular %>% ldply(function(x) (table(x$Bonferroni == T) %>% dimnames() %>% .[[1]] %>% length()) ==2)  %>% .[,2]
Drift.results.Toplot$All.sp$Results$Corr.Contrasts.1 <- Drift.results$all.sp$Correlation.test.Contrasts %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) !=1) %>% .[,2]
Drift.results.Toplot$Extants$Results$Node.ref <- Drift.results$extant.sp$Correlation.test.Regular %>% ldply(function(x) dim(table(x$Bonferroni == T, exclude = T)) ==2) %>% .[,1]
Drift.results.Toplot$Extants$Results$Corr.Ed.1 <- Drift.results$extant.sp$Correlation.test.Regular %>% ldply(function(x) (table(x$Bonferroni == T) %>% dimnames() %>% .[[1]] %>% length()) ==2)  %>% .[,2]
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


Plot.Drift.Results <- function (tree, node.ref, info, font.size = 0.4, contrasts = TRUE,  paleta.deriva = c("#ffaf02ff", "#1baab7ff", "#c3bbbfff") ) {
  node.ref <- as.numeric(info$Node.ref)
  results.ed <- data.frame(V1 = as.numeric(info$Corr.Ed.1) * 0.5, V2 = as.numeric(info$Regression.Ed) * 0.5  )
  if (contrasts == TRUE) results.ci <- data.frame(V1 = as.numeric(info$Corr.Contrasts.1) * 0.5, V2 = as.numeric(info$Regression.Contrasts) * 0.5) 
  
  results.ed$V3 <- abs(results.ed$V1 + results.ed$V2 -1)
  if (contrasts == TRUE) results.ci$V3 <- abs(results.ci$V1 + results.ci$V2 -1)
  results.ed<- as.matrix(results.ed)
  if (contrasts == TRUE) results.ci<- as.matrix(results.ci)
  if (contrasts == TRUE) par(mfrow = c(1,2)) else par(mfrow = c(1,1)) 
  plot.phylo(tree, font = 3, no.margin = TRUE,  cex = font.size, edge.color = "darkgrey", edge.width = 3, label.offset = 1 )
  nodelabels(node = node.ref, 
             pie = results.ed[,1:3], cex=0.6, 
             piecol = paleta.deriva  )
  nodelabels(node = node.ref, bg = "transparent", frame = "n",cex = 0.6, col = "black", adj = c(1.2, 1.1) )
  legend("bottomleft", inset = .08,
         title ="Drift test \nB.matrices obtained from ED",
         text.col = "grey10",
         legend = c("Rejected - Correlation", "Rejected - Regression", "Not rejected"), 
         fill =paleta.deriva, 
         col = paleta.deriva, 
         border = "grey", box.lwd = "n",
         bg= "transparent",
         cex = 0.5)
  
  if (contrasts == TRUE) {plot.phylo(tree, font = 3, no.margin = TRUE,  cex = font.size, edge.color = "darkgrey", edge.width = 3, label.offset = 1 )
  nodelabels(node = node.ref, 
             pie = results.ci[,1:3], cex=0.6, 
             piecol = paleta.deriva )
  nodelabels(node = node.ref, bg = "transparent", frame = "n",cex = 0.6, col = "black", adj = c(1.1, 1.1) )
  legend("bottomleft", inset = 0.08,
         title ="Drift test \nB.matrices obtained from IC",
         text.col = "grey10",
         legend = c("Rejected - Correlation", "Rejected - Regression", "Not rejected"), 
         fill = paleta.deriva , 
         col = paleta.deriva, 
         border = "grey", box.lwd = "n",
         bg= "transparent",
         cex = 0.5)
 
    par(mfrow = c(1,1)) 
    return(data.frame("node" = node.ref,
           "ed" = results.ed,
           "ic" = results.ci) ) }
  
  else return(data.frame("node" = node.ref,
                         "ed" = results.ed) ) 
}
Drift.results.Toplot$Extants$Plots$Tree <- Plot.Drift.Results(tree = Trees$extant.sp.tree, 
                                                              info = Drift.results.Toplot$Extants$Results, 
                                                              font.size = 0.7, contrasts = F, paleta.deriva = c("#ffaf02ff", "#1baab7ff", "#c3bbbfff") )

Drift.results.Toplot$All.sp$Plots$Tree <- Plot.Drift.Results(tree = Trees$all.with.ed, 
                                                             info = Drift.results.Toplot$All.sp$Results, 
                                                             font.size = 0.7, contrasts = F, paleta.deriva = c("#ffaf02ff", "#1baab7ff", "#c3bbbfff") )
Correlation.Tree.plot <- plot_grid(
Drift.results$extant.sp$Correlation.test.Regular$`131`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lorisiformes\n Galagidae x Lorisidae") ,
Drift.results$extant.sp$Correlation.test.Regular$`132`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lorisidae") ,
Drift.results$extant.sp$Correlation.test.Regular$`134`$P.value.plot + theme(legend.position = "none") + ggtitle ("Nyc-Lor"), 

Drift.results$extant.sp$Correlation.test.Regular$`99`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuridae x Indridae") ,
Drift.results$extant.sp$Correlation.test.Regular$`100`$P.value.plot + theme(legend.position = "none") + ggtitle ("Indridae") ,
Drift.results$extant.sp$Correlation.test.Regular$`112`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuridae") ,

Drift.results$extant.sp$Correlation.test.Regular$`76`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lepilemuridae x Cheirogaleidae") ,
Drift.results$extant.sp$Correlation.test.Regular$`77`$P.value.plot + theme(legend.position = "none") + ggtitle ("Cheirogaleidae") ,
Drift.results$extant.sp$Correlation.test.Regular$`88`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lepilemuridae") ,
  
Drift.results$extant.sp$Correlation.test.Regular$`73`$P.value.plot + theme(legend.position = "none") + ggtitle ("Strepsirrhini") ,
Drift.results$extant.sp$Correlation.test.Regular$`74`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuriformes x Daubentonidae") , 
Drift.results$extant.sp$Correlation.test.Regular$`75`$P.value.plot + theme(legend.position = "none") + ggtitle ("Lemuriformes") ,
Drift.results$extant.sp$Correlation.test.Regular$`71`$P.value.plot + theme(legend.position = "none") + ggtitle ("Prosimian\n Strepsirrhini + Tarsiidae") ,
ncol = 3)

Correlation.Tree.plot 
save_plot(filename = "Figures/Correlation_Tree_plot.pdf", plot = Correlation.Tree.plot, 
          base_aspect_ratio = 0.3, base_height = 15, base_width = 9.5)

Drift.results$extant.sp$Regression.test$`71`$log.between_group_variance

Regression.Tree.plot <- plot_grid(
  Drift.results$extant.sp$Regression.test$`134`$plot + ggtitle ("Nyc-Lor"), 
  Drift.results$extant.sp$Regression.test$`132`$plot + ggtitle ("Lorisidae") ,
  Drift.results$extant.sp$Regression.test$`131`$plot +  ggtitle ("Lorisiformes\n Galagidae x Lorisidae") ,
  Drift.results$extant.sp$Regression.test$`112`$plot +  ggtitle ("Lemuridae") ,
  Drift.results$extant.sp$Regression.test$`100`$plot +  ggtitle ("Indridae") ,
  Drift.results$extant.sp$Regression.test$`99`$plot +  ggtitle ("Lemuridae x Indridae") ,
  Drift.results$extant.sp$Regression.test$`88`$plot + ggtitle ("Lepilemuridae") ,
  Drift.results$extant.sp$Regression.test$`77`$plot + ggtitle ("Cheirogaleidae") ,
  Drift.results$extant.sp$Regression.test$`76`$plot +  ggtitle ("Lepilemuridae x Cheirogaleidae") ,
  Drift.results$extant.sp$Regression.test$`75`$plot +  ggtitle ("Lem-Ind x Lep-Che") ,
  Drift.results$extant.sp$Regression.test$`74`$plot +  ggtitle ("Lemuriformes\n Lemurs x Daubentonidae") ,
  Drift.results$extant.sp$Regression.test$`73`$plot +  ggtitle ("Strepsirrhini") ,
  Drift.results$extant.sp$Regression.test$`71`$plot +  ggtitle ("Prosimian\n Strepsirrhini + Tarsiidae") ,
  ncol = 3)

Regression.Tree.plot 

Drift.results$extant.sp$Regression.test$`101`$plot + geom_abline(slope = 1, color = "red") + ggtitle ("Indridae")

Plot.Drift.regression <- function(regress.result) {
  regressae <- data.frame( "log.B_variance" = regress.result$log.between_group_variance, "log.W_eVals" = regress.result$log.W_eVals, names = 1:(dim(mean.array)[2]) ) 
  
  if (regress.result$drift_rejected == T) .beta <- "red" 
  if (regress.result$drift_rejected == F) .beta <- "green"
  empirical.c.i <- regress.result$coefficient_CI_95[2, ]
  beta.coef <- lm(data = regressae, log.B_variance ~ log.W_eVals)$coefficients[2]
  regressae.1 <- regressae[-1,]
  IC95_.1 <- confint(lm(data = regressae[-1,], log.B_variance ~ log.W_eVals))[2,]
  containsOne <- function(x) ifelse(x[1] < 1 & x[2] > 1, TRUE, 
                                    FALSE)
  test <- !containsOne(confint(lm(data = regressae[-1,], log.B_variance ~ log.W_eVals))[2,])
  names(test) <- "5 %"
  if (test == T) .beta1 <- "red" 
  if (test == F) .beta1 <- "green"
  
  plotows <- 
    ggplot(regressae, aes (y = log.B_variance, x = log.W_eVals) ) +
    geom_abline(slope = beta.coef, color = .beta1, size = 1, alpha = 0.4, linetype = 2) +
    geom_text(aes(label = names, size = 5)) + 
    geom_smooth(method = "lm", color = .beta) + 
    labs(x = "log(W Eigenvalues)", y = "log(B variances)") + 
    theme_bw() +
    theme(legend.position = "none") +
    geom_abline(slope = 1, color = "black") 
  
  tabelado <- cbind("IC-" = empirical.c.i[1], "beta"= regress.result$regression$coefficients[2], "IC+" = empirical.c.i[2], "Ho" = regress.result$drift_rejected,
                    "IC-.-PC1" = IC95_.1[1], "beta.-PC1" = beta.coef, "IC+.-PC1" = IC95_.1[2], "Ho.-PC1" = test )
  
  
  return(list("plot" = plotows,
              "beta" = regress.result$regression$coefficients[2], 
              "beta.-1" = beta.coef, 
              "IC95" = empirical.c.i, 
              "IC95.-1" = IC95_.1, 
              "drift.rejected" = regress.result$drift_rejected,
              "drift.rejected.-1" = test,
              "intercept" = regress.result$regression$coefficients[1],
              "tabelinha" = tabelado) )
  rm(.beta, .beta.1)
}

Plot.Drift.regression (regress.result = Drift.results$extant.sp$Regression.test$`76`)

Drift.results.Toplot$Extants$Results$Regression <- Drift.results$extant.sp$Regression.test %>% llply(Plot.Drift.regression)
tabelinha <- Drift.results.Toplot$Extants$Results$Regression %>% ldply(function (x) x$tabelinha)
rownames(tabelinha) <- tabelinha[,1]

tabelinha[-1] %>% xtable()

Regression.Tree.plot <- plot_grid(
  Drift.results.Toplot$Extants$Results$Regression$`134`$plot  + ggtitle ("Nycticebus x Loris"), 
  Drift.results.Toplot$Extants$Results$Regression$`132`$plot  + ggtitle ("Lorisidae") ,
  Drift.results.Toplot$Extants$Results$Regression$`131`$plot + ggtitle ("Lorisiformes") ,
  Drift.results.Toplot$Extants$Results$Regression$`112`$plot + ggtitle ("Lemuridae") ,
  Drift.results.Toplot$Extants$Results$Regression$`100`$plot + ggtitle ("Indriidae") ,
  Drift.results.Toplot$Extants$Results$Regression$`99`$plot  + ggtitle ("Lemuridae x Indriidae") ,
  Drift.results.Toplot$Extants$Results$Regression$`88`$plot  + ggtitle ("Lepilemuridae") ,
  Drift.results.Toplot$Extants$Results$Regression$`77`$plot  + ggtitle ("Cheirogaleidae") ,
  Drift.results.Toplot$Extants$Results$Regression$`76`$plot  + ggtitle ("Lepilemuridae x Cheirogaleidae") ,
  Drift.results.Toplot$Extants$Results$Regression$`75`$plot  + ggtitle ("Lemuriformes") ,
  Drift.results.Toplot$Extants$Results$Regression$`74`$plot  + ggtitle ("Lemuriformes x Daubentonidae") ,
  Drift.results.Toplot$Extants$Results$Regression$`73`$plot  + ggtitle ("Strepsirrhini") ,
  Drift.results.Toplot$Extants$Results$Regression$`71`$plot  + ggtitle ("Prosimian\n Strepsirrhini + Tarsiidae") ,
  ncol = 3)

Regression.Tree.plot 

save_plot(filename = "R/Figures/Regression_Tree_plot.pdf", plot = CRegression.Tree.plot, 
          -          base_aspect_ratio = 0.9, base_height = 14, base_width = 4)

############### make tables
Drift.results$extant.sp$Regression.test$`71`$regression$coefficients[2]
Drift.results$extant.sp$Regression.test$`71`$coefficient_CI_95
Drift.results$extant.sp$Regression.test$`71`$drift_rejected
Drift.results$extant.sp$Regression.test$`71`$regression$coefficients[1]
Drift.results$extant.sp$Regression.test$`71`$log.W_eVals

Drift.results.Toplot$Extants$Results$Regression$`71`
Drift.results$extant.sp$Regression.test$`135`$means
  
table.drift.results <- Drift.results.Toplot$Extants$Results$Regression %>% ldply(function(x) cbind("2.5 %" = round(x$IC95[1],3), 
                                                                            "97.5 %" = round(x$IC95[2],3),
                                                                            "\beta" = round(x$beta.1,3),  
                                                                            "\beta_-PC1"= round(x$`beta.-1`,3), 
                                                                            "Drift rejected" = as.character(x$drift.rejected)
                                                                           ) 
                                                                    )
                                                  
table.drift.results$intercept <- Drift.results$extant.sp$Regression.test %>% ldply(function(x) beta = round(x$regression$coefficients[1],3 ) )%>% .[,2]
table.drift.results$B.W <- round(Bvar.Wvar$ratio, 3)
names(table.drift.results)[1] <- "Node"
rownames(table.drift.results) <- table.drift.results[,1]

table.drift.results[,-1] %>% as.data.frame %>% xtable()


