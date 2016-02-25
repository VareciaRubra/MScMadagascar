### Remover efeito do tamanho amostral da matriz de similaridade ###

#Perdoe-me pelo script portugues... Tenho certeza que voce pode fazer melhor!

#Mas da para voce ter uma ideia: mandei um for para calcular as medias harmonicas (n/(soma dos reciprocos)) para cada
#especie

harmonic.n.size <- c(n.size[mask,2], 130, 12) ## vetor com o tamanho amostral de cada especie que entrou nas comparaçao

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
Harm.Mx<- harm_matrix

cor.test(Matrix.Similarity[upper.tri(Matrix.Similarity, diag = F)], harm_matrix[upper.tri(harm_matrix, diag = F)])

cor.test(Matrix.Similarity[lower.tri(Matrix.Similarity, diag = F)], harm_matrix[upper.tri(harm_matrix, diag = F)])



#Pearson's product-moment correlation
# 
# data:  mx.compare$BS.RS$correlations[upper.tri(mx.compare$BS.RS$correlations,  and harm_matrix[upper.tri(harm_matrix, diag = F)]    diag = F)] and harm_matrix[upper.tri(harm_matrix, diag = F)]
# t = 13.375, df = 944, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3441578 0.4513957
# sample estimates:
#       cor 
# 0.3991408 



model.harmonic.mean <- lm(mx.compare$BS.RS$correlations [upper.tri(mx.compare$BS.RS$correlations, diag = F)] ~ harm_matrix[upper.tri(harm_matrix, diag = F)])
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
###################################################################################################################
Similaridade.Mx <- mx.compare$MC.RS$correlations [upper.tri(mx.compare$MC.RS$correlations, diag = F)]
model.harmonic.mean <- lm(Similaridade.Mx ~  log(Harm.Mean))
plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS', xlim=c(2.5,4.8), ylim = c(0.2, 0.9)) 
abline(coef = model.harmonic.mean$coefficients, col = 'red')
title(main = 'Comparison value corrected by MC Repetability')


Harm.Mean <- harm_matrix[-c(43:44), -c(43:44)]
Harm.Mean <- Harm.Mean[upper.tri(Harm.Mean, diag = F)]
par(mfrow= c(1,2))
Similaridade.Mx  <- mx.compare$MC.RS$correlations[-c(43:44), -c(43:44)]
Similaridade.Mx <- Similaridade.Mx [lower.tri(Similaridade.Mx, diag = F)]
model.harmonic.mean <- lm(Similaridade.Mx ~  log(Harm.Mean))
plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS', xlim=c(2.5,4.8), ylim = c(0.2, 0.9)) 
abline(coef = model.harmonic.mean$coefficients, col = 'red')
title(main = 'Raw value of comparisson')

Similaridade.Mx  <- mx.compare$MC.RS$correlations[-c(43:44), -c(43:44)]
Similaridade.Mx <- Similaridade.Mx [upper.tri(Similaridade.Mx, diag = F)]
model.harmonic.mean <- lm(Similaridade.Mx ~  log(Harm.Mean))
plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS', xlim=c(2.5,4.8), ylim = c(0.2, 0.9)) 
abline(coef = model.harmonic.mean$coefficients, col = 'red')
title(main = 'Comparisson value corrected by MC Repetability')



plot(y = Similaridade.Mx, x = log(Harm.Mean), pch = 20, col = 'grey', xlab = "log sample sizes' harmonic mean", ylab = 'Matrices similarity via RS' ) 

temp5 <- array(0, c(44,44) ) #array de dez linhas por dez colunas (tenho 10 especies)
temp5[upper.tri(temp5)] <- model.harmonic.mean$residuals #pegando os residuos da regressao: variacao na similaridade nao explicada
#pela media harmonica dos n
CalculateMatrix(lm(mx.compare$BS.RS$correlations [upper.tri(mx.compare$BS.RS$correlations)] ~ harm_matrix[upper.tri(harm_matrix)]))

mx.compare$BS.RS$correlations.residual5 <- temp5 + t(temp5) #tornando os residuos uma matriz quadrada

#Agora voce pode correlacionar uma matriz de similaridade sem efeito do n com qualquer fator

############################# Vendo efeito de filogenia ##########################################
notat.tree<- is.na(match(dimnames(mx.compare$RS$correlations)[[1]], treefile$tip.label)) 
names.at.tree <- dimnames(mx.compare$RS$correlations)[[1]][!notat.tree] 
mx.all.at.tree<- mx.compare$RS$correlations[!notat.tree,!notat.tree]
pruned.tree.all<- drop.tip(treefile,treefile$tip.label[-match( names.at.tree, treefile$tip.label)])
phylo.dist.all.at.tree<- cophenetic.phylo(pruned.tree.all)
phylo.dist.all.at.tree<- phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)]

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


MatrixCor(as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)], mx.compare$BS.RS$correlations[-c(41,43:44), -c(41,43:44)])

Phylo.Dist<- phylo.dist.all.at.tree[rownames( mx.all.at.tree), rownames( mx.all.at.tree)]
Phylo.Dist <- Phylo.Dist[upper.tri(Phylo.Dist, diag=F)]
Sim.Mx <- mx.compare$BS.RS$correlations[-c(41, 43:44), -c(41, 43:44)]
Sim.Mx <- Sim.Mx[upper.tri(t(Sim.Mx), diag = F)]
Harm.M <- as.matrix(harm_matrix)[-c(41, 43:44), -c(41, 43:44)]
Harm.M  <- Harm.M[upper.tri(Harm.M, diag = F)]
str(Sim.Mx)
str(Phylo.Dist)
str(Harm.M)

COr.PHyHM <- cbind(Sim.Mx, Harm.M, Phylo.Dist)
COr.PHyHM <- as.data.frame(COr.PHyHM)
names(COr.PHyHM) <- c("Matrix.Similarity.RS", "Harmonic.Mean", "Phylogenetic.Distance")

 Plot.Phy.SIM<- 
 COr.PHyHM %>% 
  ggplot( ., aes(x = Phylogenetic.Distance, y = Matrix.Similarity.RS), varwidth = T) +
  geom_point(color = "grey") +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_text(x = 26,  y = 0.4, label = "r.squared = 0.35 \n 818 DF,  p-value: < 2.2e-16 \n MatrixCor = -0.5770781 ", color = "red", size = 2) + 
  #ggtitle("Geometric mean x PC1.percent") +
  theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

summary(lm(Matrix.Similarity.RS ~ Phylogenetic.Distance , data = COr.PHyHM))


Plot.HM.SIM<- 
  COr.PHyHM %>% 
  ggplot( ., aes(x = log(Harmonic.Mean), y = Matrix.Similarity.RS), varwidth = T) +
  geom_point(color = "grey") +
  stat_smooth(method="lm", aes(group=1), color = "red") +
  #scale_x_continuous(limits = c(0.1, 0.7)) +
  #scale_y_continuous(limits = c(1, 7)) +
  theme_bw() +
  geom_text(x = 3.8,  y = 0.4, label = "r.squared = 0.328 \n 818 DF,  p-value: < 2.2e-16 \n MatrixCor = 0.602", color = "red", size = 2) + 
  #ggtitle("Geometric mean x PC1.percent") +
  #theme(legend.position  ="none") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

summary(lm(Matrix.Similarity.RS ~ log(Harmonic.Mean) , data = COr.PHyHM))


plot_grid(Plot.GM.PC1, Plot.R2.PC1, Plot.HM.SIM, Plot.Phy.SIM, labels = LETTERS[1:4], ncol = 2)
