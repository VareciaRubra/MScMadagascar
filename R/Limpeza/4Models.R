
#############################################
############# PRIMEIRA PARTE ################
##### testando modelos de efeitos fixos #####
#############################################
options(contrasts=c("contr.sum", "contr.poly"))

current.data <-Gen.raw.main.data$Archaeolemur


  y = vector("list", 4)
  y[[1]] = vector("list", 4)
  y[[1]] [[1]] <- manova(as.matrix(current.data$ed)  ~ Especie, data = as.data.frame(current.data$info) ) 
  y[[1]] [[2]] <- Manova(y[[1]] [[1]], type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly"))
  y[[1]] [[3]]  <- apply(as.matrix(current.data$ed ), 2, function (x) return (Anova(lm(x ~ Especie, data = as.data.frame(current.data$info) ), type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))))
  y[[1]] [[4]] <- CalculateMatrix(y[[1]] [[1]]) 
  names(y)[[1]][1:4] <- c("fit", "multi", "uni", "cov.mx")




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
  paleta  <- colfunc(200)
  
  nCores  <- length(paleta)
  intervalMarks  <- seq(from = -1, to = 1, length.out = nCores)
  dados<- as.vector(cor.mx)
  intervals  <- findInterval(dados , intervalMarks, rightmost.closed = T, all.inside = T)
  cores  <- paleta[intervals]
  
  mx.dimentions <- dim(cor.mx)[1]
  
  color2D.matplot(x = cor.mx, axes = F, cellcol = cores, show.values= show.values, vcex= 0.5, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  title(main = main , sub =  title)
  axis(1, 1:mx.dimentions, rownames(cor.mx), las = 3, cex.axis = 0.4, tick = FALSE, line = 0)
  axis(2, mx.dimentions:1, rownames(cor.mx), las = 1, cex.axis = 0.4, tick = FALSE, line = 0)
  }
  

par(mfrow= c(1,2))
plot.matrix.cor(cor.mx = cov2cor(CalculateMatrix(y[[1]][[1]])) )

plot.matrix.cor(cor.mx = all.main.data$Indri$matrix$cor)




RS.compare <- RandomSkewers(cov.x = cov.no.na, repeat.vector = rep.no.na, num.vectors = 1000, parallel = TRUE)
plot.matrix.cor(RS.compare$correlations, main = "V/CV Matrix compared by Random Skewers", brewer = "BrBG")

RS.compare.40 <- RandomSkewers(cov.x = cov.40, repeat.vector = rep.40, num.vectors = 1000, parallel = TRUE)
plot.matrix.cor(RS.compare.40$correlations, main = "V/CV Matrix compared by Random Skewers", brewer = "BrBG")
