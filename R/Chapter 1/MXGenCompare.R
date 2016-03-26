current.data <- gen.main.data

#4# Gerando a matriz ancestral
mx.at.tree <- cov.mx[mask][-41]
ancestral.mx <- PhyloW(tree = pruned.tree.with.mx, tip.data = mx.at.tree, tip.sample.size = n.size[mask,2][-41])
plot(pruned.tree.with.mx, cex = 0.5)
nodelabels()
W.matrix <- ancestral.mx$'42'


gen.cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
str(gen.cov.mx)
gen.mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
gen.n.size <- current.data %>% ldply(function(x) x$sample.size) 

gen.mask <- ldply(gen.cov.mx, function(x) !is.na(x[1]))[,2]

gen.cov.list <- gen.cov.mx[gen.mask]
gen.cov.list$Saguinus_P.cov <- Saguinus_P.cov
gen.cov.list$Saguinus_G.cov <- Saguinus_G.cov
str(gen.cov.list)
names(gen.cov.list)

# Comparando as matrizes estimadas pelo Phylow com as matrizes estimadas com os residuos do efeito de especie####
RandomSkewers(gen.cov.list$Tarsius, ancestral.mx$'43')
RandomSkewers(gen.cov.list$Microcebus, ancestral.mx$'51')
RandomSkewers(gen.cov.list$Mirza, ancestral.mx$Mirza_coquereli)
RandomSkewers(gen.cov.list$Cheirogaleus, ancestral.mx$'52')
RandomSkewers(gen.cov.list$Phaner, ancestral.mx$Phaner_furcifer)
RandomSkewers(gen.cov.list$Lepilemur, ancestral.mx$'53')
RandomSkewers(gen.cov.list$Avahi, ancestral.mx$'58')
RandomSkewers(gen.cov.list$Propithecus, ancestral.mx$'59')
RandomSkewers(gen.cov.list$Indri, ancestral.mx$Indri_indri)
RandomSkewers(gen.cov.list$Eulemur, ancestral.mx$'65')
RandomSkewers(gen.cov.list$Hapalemur, ancestral.mx$Hapalemur_griseus)
RandomSkewers(gen.cov.list$Prolemur, ancestral.mx$Prolemur_simus)
RandomSkewers(gen.cov.list$Lemur, ancestral.mx$Lemur_catta)
RandomSkewers(gen.cov.list$Varecia, ancestral.mx$'76')
RandomSkewers(gen.cov.list$Daubentonia, ancestral.mx$'45')
RandomSkewers(gen.cov.list$Perodicticus, ancestral.mx$'79')
RandomSkewers(gen.cov.list$Loris, ancestral.mx$'79')
RandomSkewers(gen.cov.list$Nycticebus, ancestral.mx$'79')
RandomSkewers(gen.cov.list$Euoticus, ancestral.mx$'80')
RandomSkewers(gen.cov.list$Otolemur, ancestral.mx$'81')
RandomSkewers(gen.cov.list$Galago, ancestral.mx$'81')

Gen.cov.list <- list()
Gen.cov.list$Tarsius <- ancestral.mx$'43'
Gen.cov.list$Microcebus <- ancestral.mx$'51'
Gen.cov.list$Mirza <- ancestral.mx$Mirza_coquereli
Gen.cov.list$Cheirogaleus <- ancestral.mx$'52'
Gen.cov.list$Phaner <- ancestral.mx$Phaner_furcifer
Gen.cov.list$W.Cheirogaleidae <- ancestral.mx$'48'
Gen.cov.list$Lepilemur <- ancestral.mx$'53'
Gen.cov.list$Avahi <- ancestral.mx$'58'
Gen.cov.list$Propithecus <- ancestral.mx$'59'
Gen.cov.list$Indri <- ancestral.mx$Indri_indri
Gen.cov.list$W.Indridae <- ancestral.mx$'56'
Gen.cov.list$Eulemur <- ancestral.mx$'65'
Gen.cov.list$Hapalemur <- ancestral.mx$Hapalemur_griseus
Gen.cov.list$Prolemur <- ancestral.mx$Prolemur_simus
Gen.cov.list$Lemur <- ancestral.mx$Lemur_catta
Gen.cov.list$Varecia <- ancestral.mx$'76'
Gen.cov.list$W.Lemuridae <- ancestral.mx$'63'
Gen.cov.list$Daubentonia <- ancestral.mx$Daubentonia_madagascariensis
Gen.cov.list$W.Madagascar <- ancestral.mx$'45'
Gen.cov.list$Perodicticus <- ancestral.mx$Perodicticus_potto
Gen.cov.list$Loris <- ancestral.mx$Loris_tardigradus
Gen.cov.list$W.Lorisidae <- ancestral.mx$'78'
Gen.cov.list$Nycticebus <- ancestral.mx$Nycticebus_coucang
Gen.cov.list$Euoticus <- ancestral.mx$Euoticus_elegantulus
Gen.cov.list$Otolemur <- ancestral.mx$Otolemur_crassicaudatus
Gen.cov.list$Galago <- ancestral.mx$Galago_senegalensis
Gen.cov.list$W.Galagidae <- ancestral.mx$'80'
Gen.cov.list$W.OutMadagascar <- ancestral.mx$'77'
Gen.cov.list$W.Strepsirrhini <- ancestral.mx$'44'
Gen.cov.list$W.Prosimian <- ancestral.mx$'42'
Gen.cov.list$Saguinus.P <- Saguinus_P.cov
Gen.cov.list$Saguinus.G <- Saguinus_G.cov


ComparaSagui <- list ("G" = RandomSkewers(ancestral.mx[c(1:41, 45, 77,63, 56, 53, 48, 43)], Saguinus_G.cov),
                      "P" =  RandomSkewers(ancestral.mx[c(1:41, 45, 77,63, 56, 53, 48, 43)], Saguinus_P.cov) ) 
ComparaSagui %>% ldply(., function(x) summary(x$correlation[1:41]) ) %>% xtable(caption )

#names(gen.cov.list) %<>% gsub("_", ' ', .)

Gen.mx.compare = vector("list", 2)
Gen.mx.compare[1:2] <- NA
Gen.mx.compare[[1]] <- RandomSkewers(cov.x = Gen.cov.list, num.vectors = 1000 )
Gen.mx.compare[[2]]$correlations <- as.matrix(KrzCor(cov.x= Gen.cov.list, num.vectors = 1000))
names(Gen.mx.compare)[1:2] <-  c('RS', 'KRZ')

Steppan$Plot$GenusCompare <- 
  Combine.Mx.Plot(Mx1 = t(Gen.mx.compare$RS$correlations),
                Mx2 = t(Gen.mx.compare$KRZ$correlations),
                prob = Gen.mx.compare$RS$probabilities,
                diag.info = rep(NA, dim(Gen.mx.compare$RS$correlations)[[1]]), 
                method =  "Random Skewers",
                titulo = "Genus Matrices raw values of comparison via \nKrzanowski and")
Steppan$Plot$GenusCompare + theme(legend.position= "bottom") + theme(legend.position= c(-0.1,-0.15)) 


par (mar = c(5, 4, 4, 2) + 0.1)
plot(pruned.tree.with.mx, cex = 0.9, no.margin =T, edge.color = "darkgrey", edge.width = 2, x.lim = c(-20, 150) )
#nodelabels()
nodelabels(bg = NULL, frame = "none", node = 56, pch = NULL, col = "black", text = "Indroidea" )
nodelabels(bg = NULL, frame = "none", node = 45, pch = NULL, col = "black", text = "Malagasy \nprimates" )
nodelabels(bg = NULL, frame = "none", node = 77, pch = NULL, col = "black", text = "Out of     \nMadagascar", adj = c(0.8, 0.5) )
nodelabels(bg = NULL, frame = "none", node = 42, pch = NULL, col = "black", text = "Prosimian", adj = c(0.8, 0.5) )
nodelabels(bg = NULL, frame = "none", node = 44, pch = NULL, col = "black", text = "Strepsirrhini", adj = c(0.8, 0.5) )
nodelabels(bg = NULL, frame = "none", node = 48, pch = NULL, col = "black", text = "Cheirogaleidae", adj = c(0.4, 0.5))
nodelabels(bg = NULL, frame = "none", node = 63, pch = NULL, col = "black", text = "Lemuridae" )
nodelabels(bg = NULL, frame = "none", node = 80, pch = NULL, col = "black", text = "Galagidae", adj = c(0.1, 0.5))
nodelabels(bg = NULL, frame = "none", node = 78, pch = NULL, col = "black", text = "Lorisidae", adj = c(0.1, 0.5))


gm.means[,1] == pruned.tree.with.mx$tip.label
rownames(gm.means) <- gm.means[,1]

ace(x = gm.means[,2], phy = pruned.tree.with.mx)
ace(x = gm.means[,2], phy = pruned.tree.with.mx, method = "pic")

ace.gm <- ace(x = gm.means[,2], phy = pruned.tree.with.mx, type = "continuous", CI = T, method = "GLS", use.eigen = T, corStruct = corBrownian(1,  pruned.tree.with.mx))
rownames(gm.means) == pruned.tree.with.mx$tip.label
plot(pruned.tree.with.mx, cex = 0.9, no.margin =T, edge.color = "darkgrey", edge.width = 2, x.lim = c(-20, 150) )
co <- c("red", "blue")
tiplabels(pch = 22, bg = co[as.numeric(gm.means[,2])], cex = 2, adj = 1)
nodelabels(thermo = ace.gm$lik.anc, piecol = co, cex = 0.75)
