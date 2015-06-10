###################### testando modelos para possiveis correções de efeitos fixos ##################
options(contrasts=c("contr.sum", "contr.poly"))

current.data <- all.main.data$Hapalemur
  y = vector("list", 4)
  y[[1]] <- manova(as.matrix(current.data$ed) ~ Especie, data = as.data.frame(current.data$info))
  y[[2]] <- Manova(y[[1]] , type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly"))
  y[[3]] <- apply(as.matrix(current.data$ed), 2, function (x) return (Anova(lm(x ~ Especie, data = as.data.frame(current.data$info) ), type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))))
  y[[4]] <- CalculateMatrix(y[[1]]) 
  names(y)[1:3] <- c("fit", "multi", "uni", "cov.mx")