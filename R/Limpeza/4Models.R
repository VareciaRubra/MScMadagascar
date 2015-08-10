
#############################################
############# PRIMEIRA PARTE ################
##### testando modelos de efeitos fixos #####
#############################################
options(contrasts=c("contr.sum", "contr.poly"))

current.data <- sp.main.data$Varecia_variegata
  y = vector("list", 7)
  y[[1]] = vector("list", 10)
  y[[1]] [[1]] <- manova(as.matrix(current.data$ed)  ~ Sexo, data = as.data.frame(current.data$info) ) 
  y[[1]] [[2]] <- Manova(y[[1]] [[1]], type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly"))
  y[[1]] [[3]]  <- apply(as.matrix(current.data$ed ), 2, function (x) return (Anova(lm(x ~ Sexo, data = as.data.frame(current.data$info) ), type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))))
  y[[1]] [[4]] <- CalculateMatrix(y[[1]] [[1]]) 
  y[[1]] [[5]] <- current.data$matrix$cov
  names(y)[1] <- "factors"
  names(y[[1]])[1:5] <- c("fit", "multi", "uni",  "cov.fit.mx", "cov.mx")
  y$RS <- RandomSkewers (cov.x = y$factors$cov.fit.mx, cov.y = y$factors$cov.mx , num.vectors = 10000, parallel = TRUE)
  y$KRZ <- KrzProjection(cov.x = y$factors$cov.fit.mx, cov.y = y$factors$cov.mx )$total.variation
  y$PCA.s <- PCAsimilarity(cov.x = y$factors$cov.fit.mx, cov.y = y$factors$cov.mx )
  y$cor.KRZ <- KrzProjection(cov.x = cov2cor(y$factors$cov.fit.mx), cov.y = cov2cor(y$factors$cov.mx) ) $total.variation
  y$cor.PCA.s <- PCAsimilarity(cov.x = cov2cor(y$factors$cov.fit.mx), cov.y = cov2cor(y$factors$cov.mx) )
  
  MeanMatrixStatistics(y$factors$cov.mx)
  plot.matrix.cor(y[[1]] [[5]])
  

  sig.test <-  Manova(manova(as.matrix(x$ed) ~ Sexo, data = x$info), type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly"))
  summary.Anova <- summary(Anova(lm(as.matrix(sp.main.data$Microcebus_griseorufus$ed)~ Sexo, data = sp.main.data$Microcebus_griseorufus$info) , type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly") ) )
  Linhas <- cat(capture.output(summary.Anova)[], sep = "\n")
#############################################
############## SEGUNDA PARTE ################
## MATRIZ de genero: residual de especies  ##
#############################################


require (plotrix)

RandomSkewers (CalculateMatrix(y[[1]]), 
               var (current.data$ed [current.data$info$Museu != 'MCZ', ]))

RandomSkewers (CalculateMatrix(y[[1]]), 
               var (current.data$ed))

plot.matrix.cor<- function(cor.mx = NULL, main = "", title = "", brewer = "BrBG", show.values = TRUE )
  {
  library("RColorBrewer")
  paleta.of.choice  <- rev(brewer.pal(11, brewer))
  colfunc <- colorRampPalette(paleta.of.choice)
  paleta  <- colfunc(90)
  
  nCores  <- length(paleta)
  intervalMarks  <- seq(from = 0.1, to = 1, length.out = nCores)
  dados<- as.vector(cor.mx)
  intervals  <- findInterval(dados , intervalMarks, rightmost.closed = T, all.inside = T)
  cores  <- paleta[intervals]
  
  mx.dimentions <- dim(cor.mx)[1]
  
  color2D.matplot(x = cor.mx, axes = F, cellcol = cores, show.values= show.values, vcex= 0.8, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  title(main = main , sub =  title)
  axis(1, 1:mx.dimentions, paste("(n=", sample.40,")",1:length(rownames(cor.mx)), sep = " "), las = 2, cex.axis = 0.7, tick = FALSE, line = 0)
  #axis(2, mx.dimentions:1, length(rownames(cor.mx)):1, las = 1, cex.axis = 1, tick = FALSE, line = 0)
  
  axis(2, mx.dimentions:1, paste(rownames(cor.mx), 1:18), las = 1, cex.axis = 1, tick = FALSE, line = 0)
  }
plot.matrix.cor(RS.compare.40$correlations, main = "V/CV Matrix compared by Random Skewers", brewer = "RdBu")


par(mfrow= c(1,2))
plot.matrix.cor(cor.mx = cov2cor(CalculateMatrix(y[[1]][[1]])) )

plot.matrix.cor(cor.mx = all.main.data$Indri$matrix$cor)




RS.compare <- RandomSkewers(cov.x = cov.no.na, repeat.vector = rep.no.na, num.vectors = 1000, parallel = TRUE)
plot.matrix.cor(RS.compare$correlations, main = "V/CV Matrix compared by Random Skewers", brewer = "BrBG")

RS.compare.40 <- RandomSkewers(cov.x = cov.40, repeat.vector = rep.40, num.vectors = 1000, parallel = TRUE)
plot.matrix.cor(RS.compare.40$correlations, main = "V/CV Matrix compared by Random Skewers", brewer = "RdBu")


