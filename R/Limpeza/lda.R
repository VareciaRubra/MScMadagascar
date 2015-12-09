library(MASS)

ordered<- info[info$Tombo %in% dimnames(MEAN.ind)[[1]], ]

restricted.info<- ordered[ordered$Especie %in% c("Avahi laniger","Avahi peyrierasi "), ]
restricted.ed<- MEAN.ind[rownames(MEAN.ind) %in% restricted.info$Tombo,]
table(ordered$Especie)
table(restricted.info$Especie)


lfunction<-lda(restricted.ed, restricted.info$Especie)
# com o resultado disso pegar o objeto $rotation (ou algo do tipo) que tem o peso de cada caractere 
#Isso vai tar no espaço das canonicas. pra voltar pro espaço das variaveis originais tem que pegar essa mx de rotação pela inversa da mx v/cv
predicted.scores  <- predict(lfunction, restricted.ed)

write.csv2(coef(lfunction), file="coef_lda.csv")

source('confusion.R')

cross.v<-lda(restricted.ed, restricted.info$Especie, CV=TRUE)
confu<-confusion(restricted.info$Especie, predicted.scores$class, printit=FALSE)
confu.CV<-confusion(restricted.info$Especie, cross.v$class, printit=FALSE)

ordered$Especie <- as.factor(ordered$Especie)
restricted.info$Especie<- as.factor(restricted.info$Especie)

confusion.ld<-matrix(0,length(levels(restricted.info$Especie)),length(levels(restricted.info$Especie)))
rownames(confusion.ld)<-levels(restricted.info$Especie)
colnames(confusion.ld)<-levels(restricted.info$Especie)

nsp<-length(levels(restricted.info$Especie))
for(i in 1:nsp) for(j in 1:nsp){
  confusion.ld[i,j]<-paste(round(confu[i,j],2),"(", round(confu.CV[i,j],2),")", sep="")
}

write.csv2(confusion.ld, file="confusion.csv")

confusion
