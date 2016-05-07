### Remover efeito do tamanho amostral da matriz de similaridade ###

#Perdoe-me pelo script portugues... Tenho certeza que voce pode fazer melhor!

#Mas da para voce ter uma ideia: mandei um for para calcular as medias harmonicas (n/(soma dos reciprocos)) para cada
#especie

harmonic.n.size <- c(n.size[mask,2], 130, 36) ## vetor com o tamanho amostral de cada especie que entrou nas comparaçao

names(harmonic.n.size) <- c(n.size[mask,1], 'Saguinus_P', 'Saguinus_G')

#matriz de medias harmonicas
harm_matrix <- matrix(NA,nrow = length(harmonic.n.size), ncol = length(harmonic.n.size) )
for (j in 1:length(harmonic.n.size)) {
  for(i in 1:length(harmonic.n.size)) { harm_matrix[i,j] <- 2/((1/harmonic.n.size[j]) + (1/harmonic.n.size[i]))} 
}

menor.n.size.matrix <- matrix(NA,nrow = length(harmonic.n.size), ncol = length(harmonic.n.size) )
for (j in 1:length(harmonic.n.size)) {
  for(i in 1:length(harmonic.n.size)) { menor.n.size.matrix[i,j] <- ifelse(harmonic.n.size[i] <= harmonic.n.size[j], log(harmonic.n.size[i]), log(harmonic.n.size[j]) ) } 
}

### Testando o efeito do n na similaridade e removendo
#correlacao entre os valores de similaridade e os valores de media harmonica
Matrix.Similarity <- mx.compare$BS.RS$correlations
Matrix.Similarity.KRZ <- mx.compare$BS.KRZ$correlations
Harm.Mx<- harm_matrix
cor.test(Matrix.Similarity[lower.tri(Matrix.Similarity, diag = F)], harm_matrix[upper.tri(harm_matrix, diag = F)])
{# #Pearson's product-moment correlation (RS)
# 
# data:  mx.compare$BS.RS$correlations[upper.tri(mx.compare$BS.RS$correlations,  and harm_matrix[upper.tri(harm_matrix, diag = F)]    diag = F)] and harm_matrix[upper.tri(harm_matrix, diag = F)]
# t = 13.375, df = 944, p-value < 0.001
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3441578 0.4513957
# sample estimates:
#       cor 
# 0.3991408 
}
cor.test(Matrix.Similarity.KRZ[lower.tri(Matrix.Similarity.KRZ, diag = F)], harm_matrix[upper.tri(harm_matrix, diag = F)])
{# Pearson's product-moment correlation (KRZ)
# 
# data:  Matrix.Similarity.KRZ[lower.tri(Matrix.Similarity.KRZ, diag = F)] and harm_matrix[upper.tri(harm_matrix, diag = F)]
# t = -0.59072, df = 944, p-value = 0.5548
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.08285989  0.04457057
# sample estimates:
#         cor 
# -0.01922273 
}

model.harmonic.mean <- lm(mx.compare$BS.RS$correlations [lower.tri(mx.compare$BS.RS$correlations, diag = F)] ~ harm_matrix[upper.tri(harm_matrix, diag = F)])
model.harmonic.mean.k <- lm(mx.compare$BS.KRZ$correlations [lower.tri(mx.compare$BS.RS$correlations, diag = F)] ~ harm_matrix[upper.tri(harm_matrix, diag = F)])

#modelo linear: similaridade por RS explicado por medias harmonicas
Harm.Mean <- harm_matrix[upper.tri(harm_matrix, diag = F)]
###################################################################################################################
######################     ESSE GRAFICO AQUI     ##################################################################
###################################################################################################################
par(mfrow= c(1,1))
Similaridade.Mx <- t(mx.compare$MC.RS$correlations) [upper.tri(mx.compare$MC.RS$correlations, diag = F)]
model.harmonic.mean <- lm(Similaridade.Mx ~  log(Harm.Mean))
plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS', xlim=c(2.5,4.8), ylim = c(0.2, 0.9)) 
abline(coef = model.harmonic.mean$coefficients, col = 'red')
title(main = 'Sample size effect in comparison values')
cor.test.values<- cor.test(Similaridade.Mx, Harm.Mean)

###### excluindo o Otolemur kirkii pq nao tem ele na arvore do James (ele é a sp. 41. 43: 44 é a P e G de Saguinus)

Harm.Mean <- harm_matrix[-c(41, 43:44), -c(41, 43:44)]
Harm.Mean <- Harm.Mean[upper.tri(Harm.Mean, diag = F)]
par(mfrow= c(1,2))
Similaridade.Mx  <- mx.compare$MC.RS$correlations[-c(41, 43:44), -c(41, 43:44)]
Similaridade.Mx <- Similaridade.Mx [lower.tri(Similaridade.Mx, diag = F)]
model.harmonic.mean <- lm(Similaridade.Mx ~  log(Harm.Mean))
plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS', xlim=c(2.5,4.8), ylim = c(0.2, 0.9)) 
abline(coef = model.harmonic.mean$coefficients, col = 'red')
title(main = 'Raw value of comparisson')

plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS' ) 

temp5 <- array(0, c(44,44) ) #array de dez linhas por dez colunas (tenho 10 especies)
temp5[upper.tri(temp5)] <- model.harmonic.mean$residuals #pegando os residuos da regressao: variacao na similaridade nao explicada
#pela media harmonica dos n
mx.compare$BS.RS$correlations.residual5 <- temp5 + t(temp5) #tornando os residuos uma matriz quadrada
mx.compare$BS.RS$correlations.residual5 <- var(mx.compare$BS.RS$correlations.residual5)
#Agora voce pode correlacionar uma matriz de similaridade sem efeito do n com qualquer fator
MatrixCor(mx.compare$BS.RS$correlations, mx.compare$BS.RS$correlations.residual5)

############################# Vendo efeito de filogenia ##########################################
notat.tree <- is.na(match(dimnames(mx.compare$MC.RS$correlations)[[1]], treefile$tip.label)) 
names.at.tree <- dimnames(mx.compare$MC.RS$correlations)[[1]][!notat.tree] 
mx.all.at.tree <- mx.compare$MC.RS$correlations[!notat.tree,!notat.tree]
pruned.tree.all <- drop.tip(treefile,treefile$tip.label[-match( names.at.tree, treefile$tip.label)])
phylo.dist.all.at.tree <- cophenetic.phylo(pruned.tree.all)
phylo.dist.all.at.tree <- phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)]

dimnames(mx.all.at.tree)[[1]] %in% dimnames(mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])
MatrixCor(phylo.dist.all.at.tree, mx.all.at.tree)
## -0.5739466
MatrixCor(phylo.dist.all.at.tree, mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])
## -0.5729755
### bom indicativo de que quanto mais próximo na filogenia mais semelhante sao as matrizes
#### vale a pena fazer um intervalo de confiança dessa porqueira?

##################################################################################################

dimnames(phylo.dist.all.at.tree)[[1]] %in% 
dimnames(mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])[[1]] 


MatrixCor(as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)], mx.compare$MC.RS$correlations[-c(41,43:44), -c(41,43:44)])

Phylo.Dist<- phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)]
Phylo.Dist <- Phylo.Dist[upper.tri(Phylo.Dist, diag=F)]
Sim.Mx.bs.rs <- mx.compare$BS.RS$correlations[-c(41, 43:44), -c(41, 43:44)]
Sim.Mx.bs.rs <- Sim.Mx.bs.rs[upper.tri(t(Sim.Mx.bs.rs), diag = F)]
Sim.Mx.bs.krz <- mx.compare$BS.KRZ$correlations[-c(41, 43:44), -c(41, 43:44)]
Sim.Mx.bs.krz <- Sim.Mx.bs.krz[upper.tri(t(Sim.Mx.bs.krz), diag = F)]
Sim.Mx.mc.rs <- mx.compare$MC.RS$correlations[-c(41, 43:44), -c(41, 43:44)]
Sim.Mx.mc.rs <- Sim.Mx.mc.rs[upper.tri(t(Sim.Mx.mc.rs), diag = F)]
Sim.Mx.mc.krz <- mx.compare$MC.KRZ$correlations[-c(41, 43:44), -c(41, 43:44)]
Sim.Mx.mc.krz <- Sim.Mx.mc.krz[upper.tri(t(Sim.Mx.mc.krz), diag = F)]
Harm.M <- as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)]
Harm.M  <- Harm.M[upper.tri(Harm.M, diag = F)]

COr.PHyHM <- cbind(Sim.Mx.bs.rs, 
                   Sim.Mx.bs.krz, 
                   Sim.Mx.mc.rs, 
                   Sim.Mx.mc.krz, 
                   Harm.M, 
                   Phylo.Dist)
COr.PHyHM <- as.data.frame(COr.PHyHM)
names(COr.PHyHM) <- c("Matrix.Similarity.RS", 
                      "Matrix.Similarity.KRZ", 
                      "Matrix.Similarity.RS.mc", 
                      "Matrix.Similarity.KRZ.mc",
                      "Harmonic.Mean", 
                      "Phylogenetic.Distance")

COr.PHyHM$Harmonic.Mean.inverso<- 1/(COr.PHyHM$Harmonic.Mean)

#######################################################################################
#Corrigindo os valores de similaridade observado do jeito do Steppan 1997 @ Evolution (I)
######### Valor corrigido = Valor observado - (slope)/ (média harmonica dos tamanhos amostrais)
############################################## slope = slope da regressao valor observado ~ 1/ (média harmonica dos tamanhos amostrais)

summary(lm(Matrix.Similarity.RS ~ COr.PHyHM$Harmonic.Mean.inverso, data = COr.PHyHM)) ### slope = -6.43707
summary(lm(Matrix.Similarity.KRZ ~ COr.PHyHM$Harmonic.Mean.inverso, data = COr.PHyHM)) ### slope = -0.389137
summary(lm(Matrix.Similarity.RS.mc ~ COr.PHyHM$Harmonic.Mean.inverso, data = COr.PHyHM)) ### slope = -6.28693 
summary(lm(Matrix.Similarity.KRZ.mc ~ COr.PHyHM$Harmonic.Mean.inverso, data = COr.PHyHM)) ### slope = -5.18420

Steppan <- list("BS.RS" = mx.compare$BS.RS$correlations - (-6.43707 /harm_matrix),
                "BS.KRZ" = mx.compare$BS.KRZ$correlations - (-0.389137 /harm_matrix),
                "MC.RS" = mx.compare$MC.RS$correlations - (-6.28693 /harm_matrix),
                "MC.KRZ" = mx.compare$MC.KRZ$correlations - (-5.18420 /harm_matrix))


COr.PHyHM$Matrix.Similarity.KRZ <- Sim.Mx.bs.krz
 Plot.Phy.SIM<- 
 COr.PHyHM %>% 
  ggplot( ., aes(x = Phylogenetic.Distance, y = Matrix.Similarity.RS), varwidth = T) +
  geom_point(color = "grey") +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_label(x = 26,  y = 0.4, label = "r.squared = 0.35 \n 818 DF,  p-value: < 0.001 \n MatrixCor = -0.58 ", size = 3) + 
  #ggtitle("Geometric mean x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
        xlab( "Phylogenetic distance") + 
        ylab("Matrix similarity by RS")

summary(lm(Matrix.Similarity.RS ~ Phylogenetic.Distance , data = COr.PHyHM))
MatrixCor(phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)], mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])

COr.PHyHM$Matrix.Similarity.KRZ <- Sim.Mx.bs.krz
Plot.Phy.SIM.k<- 
  COr.PHyHM %>% 
  ggplot( ., aes(x = sqrt(Phylogenetic.Distance), y = Matrix.Similarity.KRZ.mc), varwidth = T) +
  geom_point(color = "grey") +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme_bw() +
  geom_label(x = 2.5,  y = 0.55, label = "r.squared = 0.16 \n 818 DF,  p-value: < 0.001 \n MatrixCor = -0.44 ", size = 2) + 
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab( "Phylogenetic distance") + 
  ylab("Matrix similarity by KRZ")

summary(lm(Matrix.Similarity.KRZ.mc ~ sqrt(Phylogenetic.Distance ), data = COr.PHyHM))
MatrixCor( sqrt(phylo.dist.all.at.tree[rownames( mx.all.at.tree) , rownames( mx.all.at.tree)] ),  mx.compare$BS.KRZ$correlations[-c(41,43:44), -c(41,43:44)])

Plot.HM.SIM<- 
  COr.PHyHM %>% 
  ggplot( ., aes(x = Harmonic.Mean.inverso, y = Matrix.Similarity.RS), varwidth = T) +
  geom_point(color = "grey") +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_label(x = 0.023,  y = 0.4, label = "r.squared = 0.35 \n 818 DF,  p-value: < 0.001 \n MatrixCor = -0.65", size = 3) + 
  #ggtitle("Geometric mean x PC1.percent") +
  #theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold") ) +
        xlab( "1/ sample size harmonic mean") +
        ylab("Matrix similarity by RS")

summary(lm(Matrix.Similarity.RS ~ COr.PHyHM$Harmonic.Mean.inverso, data = COr.PHyHM))
MatrixCor(1/(as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)] ), mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])

temp <- plot_grid(Plot.HM.SIM, Plot.Phy.SIM,  labels = LETTERS[1:4], ncol = 2)

Plot.HM.PHy<- 
  COr.PHyHM %>% 
  ggplot( ., aes(x = log(Harmonic.Mean), y = Phylogenetic.Distance), varwidth = T) +
  geom_point(color = "grey") +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_text(x = 3.8,  y = 0.4, label = "r.squared = 0.05665 \n 818 DF,  p-value: < 0.001 \n MatrixCor = -0.228", color = "red", size = 2) + 
  #ggtitle("Geometric mean x PC1.percent") +
  #theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold") ) +
  xlab( "Harmonic mean") +
  ylab("Phylogenetic distance")

summary(lm(Phylogenetic.Distance ~ log(Harmonic.Mean) , data = COr.PHyHM))
MatrixCor(as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)], phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)])

<<<<<<< HEAD
=======
mantel.partial(xdis = 1-mx.compare$MC.RS$correlations[-c(41, 43:44), -c(41, 43:44)], 
               zdis = as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)], 
               ydis = phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)], 
               method = "pearson")
mantel.partial(xdis = 1-mx.compare$MC.KRZ$correlations[-c(41, 43:44), -c(41, 43:44)], 
               zdis = as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)], 
               ydis = phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)], method = "pearson")
>>>>>>> a6f0b3fa071439c4d680c4d364e9bebc4f3ed6d4
