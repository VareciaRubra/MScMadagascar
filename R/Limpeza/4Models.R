###################### testando modelos para possiveis correções de efeitos fixos ##################
options(contrasts=c("contr.sum", "contr.poly"))

current.data <- main.data$Propithecus_verreauxi
  y = vector("list", 4)
  y[[1]] <- manova(as.matrix(current.data$ed) ~ Museu, data = as.data.frame(current.data$info))
  y[[2]] <- Manova(y[[1]] , type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly"))
  y[[3]] <- apply(as.matrix(current.data$ed), 2, function (x) return (Anova(lm(x ~ Especie, data = as.data.frame(current.data$info) ), type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))))
  y[[4]] <- CalculateMatrix(y[[1]]) 
  names(y)[1:3] <- c("fit", "multi", "uni", "cov.mx")

require (plotrix)

RandomSkewers (CalculateMatrix(y[[1]]), 
               var (current.data$ed [current.data$info$Museu != 'MCZ', ]))

RandomSkewers (CalculateMatrix(y[[1]]), 
               var (current.data$ed))

par(mfrow =c(1,2) )
color2D.matplot (cov2cor(CalculateMatrix(y[[1]])), axes = F)
axis(1, 1:39, rownames(cor (current.data$ed)), las = 3, cex.axis = 0.51)
axis(2, 39:1, rownames(cor (current.data$ed)), las = 1, cex.axis = 0.51)
color2D.matplot( cor (current.data$ed), axes = F) 
axis(1, 1:39, rownames(cor (current.data$ed)), las = 3, cex.axis = 0.51)