### Remover efeito do tamanho amostral da matriz de similaridade ###

#Perdoe-me pelo script portugues... Tenho certeza que voce pode fazer melhor!

#Mas da para voce ter uma ideia: mandei um for para calcular as medias harmonicas (n/(soma dos reciprocos)) para cada
#especie

#sample size = vetor com o tamanho amostral de cada especie

temp <- rep(NA,10)
for(i in 1:11) { temp[i] <- 2/((1/sample.size[1]) + (1/sample.size[i]))} #no primeiro elemento da soma, coloquei o valor 1
#referindo-se a primeira especie, centralis

cent <- temp #vetor com as medias harmonicas somente para centralis

#Depois troquei para o numero 2, humboldti e criei o vetor humb = medias harmonicas somente para humboldti
#E assim por diante, ate ter uma matriz par a par de medias harmonicas

harm_matrix <- cbind(cent,humb,meri,gra,mira,maj,ber,pyg,dorb,fer,marg) #juntei todos os vetores de especie 
harm_matrix #matriz de medias harmonicas

### Testando o efeito do n na similaridade e removendo

cor.test(sim.cov.RS.nonadj_nosize[-11,-11] [upper.tri(sim.cov.RS.nonadj_nosize[-11,-11])],harm_matrix[upper.tri(harm_matrix)])
#correlacao entre os valores de similaridade e os valores de media harmonica
###NAO ESQUECA DE CONFERIR QUE A ORDEM DAS MATRIZES EH A MESMA!

model5 <- lm(sim.cov.RS.nonadj_nosize[-11,-11] [upper.tri(sim.cov.RS.nonadj_nosize[-11,-11])]~harm_matrix[upper.tri(harm_matrix)])
#modelo linear: similaridade por RS explicado por medias harmonicas

temp5 <- array(0,c(10,10)) #array de dez linhas por dez colunas (tenho 10 especies)
temp5[upper.tri(temp5)] <- model5$residuals #pegando os residuos da regressao: variacao na similaridade nao explicada
#pela media harmonica dos n
sim.cov.RS.nonadj_nosize.residual5 <- temp5 + t(temp5) #tornando os residuos uma matriz quadrada

#Agora voce pode correlacionar uma matriz de similaridade sem efeito do n com qualquer fator