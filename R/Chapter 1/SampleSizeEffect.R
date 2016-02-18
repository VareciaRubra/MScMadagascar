### Remover efeito do tamanho amostral da matriz de similaridade ###

#Perdoe-me pelo script portugues... Tenho certeza que voce pode fazer melhor!

#Mas da para voce ter uma ideia: mandei um for para calcular as medias harmonicas (n/(soma dos reciprocos)) para cada
#especie

harmonic.n.size <- c(n.size[mask,2], 130, 230) ## vetor com o tamanho amostral de cada especie que entrou nas comparaçao

names(harmonic.n.size) <- c(n.size[mask,1], 'Saguinus_P', 'Saguinus_G')

#matriz de medias harmonicas
harm_matrix <- matrix(NA,nrow = length(harmonic.n.size), ncol = length(harmonic.n.size) )
for (j in 1:length(harmonic.n.size)) {
  

    for(i in 1:length(harmonic.n.size)) { harm_matrix[i,j] <- 2/((1/harmonic.n.size[j]) + (1/harmonic.n.size[i]))} 
}

### Testando o efeito do n na similaridade e removendo
#correlacao entre os valores de similaridade e os valores de media harmonica
cor.test(mx.compare$BS.RS$correlations[upper.tri(mx.compare$BS.RS$correlations, diag = F)], harm_matrix[upper.tri(harm_matrix, diag = F)])

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

model.harmonic.mean <- lm(mx.compare$BS.RS$correlations [upper.tri(mx.compare$BS.RS$correlations, diag = F)]~harm_matrix[upper.tri(harm_matrix, diag = F)])
#modelo linear: similaridade por RS explicado por medias harmonicas

temp5 <- array(0,c(44,44)) #array de dez linhas por dez colunas (tenho 10 especies)
temp5[upper.tri(temp5)] <- model.harmonic.mean$residuals #pegando os residuos da regressao: variacao na similaridade nao explicada
#pela media harmonica dos n
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

################## Ponderar os valores de comparação pelo tamanho amostral : média armonica
