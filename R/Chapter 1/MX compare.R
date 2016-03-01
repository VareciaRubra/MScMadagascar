  current.data <- sp.main.data
  current.data <- gen.main.data
  
  cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
  #cov.sizeless <- current.data %>% llply(function(x) x$matrix$cov.sizeless ) #### só pra sizeless
  cor.mx <- current.data %>% llply(function(x) x$matrix$cor)
  mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
  #mx.rep.sizeless <- current.data %>% ldply(function(x) x$size.mx.rep.rs$cov.sizeless) #### só pra sizeless
  mx.rep.mc <- sp.main.data.t%>% ldply(function(x) x$Mx.Rep$MCRep) 
  mx.rep.all <- cbind(mx.rep, mx.rep.mc[,-1])
  n.size <- current.data %>% ldply(function(x) x$sample.size) 
  
  #criando mascaras para selecionar só parte do dataset.
  #Todas que tem matriz, mesmo as mal estimadas
  mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
  #Só mtrizes com mais de 40 indivíduos
  mask.n.size <- n.size[,2]>40
  
  # Add as matrizes P e G de Saguinus e suas respectivas informações 
  #Saguinus_G.cov <- read.csv(file = "Data/Saguinus_G.csv", header = F)
  #Saguinus_G.cov <- as.matrix(Saguinus_G.cov)
  #dimnames(Saguinus_G.cov) <- dimnames(cov.list[[1]])
  #Saguinus_G.cor <- cov2cor(Saguinus_G.cov)
  
  #Saguinus_P.cov <- read.csv(file = "Data/Saguinus_P.csv", header = F)
  #Saguinus_P.cov <- as.matrix(Saguinus_P.cov)
  #Saguinus_P.cov.t <- t(Saguinus_P.cov)
  #Saguinus_P.cov[upper.tri(diag(39))] <- Saguinus_P.cov.t[upper.tri(diag(39))]
  
  #dimnames(Saguinus_P.cov) <- dimnames(cov.list[[1]])
  #Saguinus_P.cor <- cov2cor(Saguinus_P.cov)
  
  cov.list <- cov.mx[mask]
  cov.list <- cov.sizeless[mask]
  cov.list$Saguinus_P.cov <- Saguinus_P.cov
  cov.list$Saguinus_G.cov <- Saguinus_G.cov
  
  cor.list <- cor.mx[mask]
  cor.list$Saguinus_P.cor <- Saguinus_P.cor
  cor.list$Saguinus_G.cor <- Saguinus_G.cor
  
  #safe.copy.rep.list <- rep.list 
  mx.rep.all <- cbind(mx.rep, mx.rep.mc[,-1])
  rep.list <- mx.rep.all[mask, ]
  #rep.list <-  mx.rep.sizeless[mask,] #### só pra sizeless
  row.names(rep.list) <- rep.list$Especie
  rep.list$Especie <- as.character(rep.list$Especie)
  rep.list[dim(rep.list)[1] +1, ] <- c("Saguinus_P", rep(0.97, 7) )
  rep.list[dim(rep.list)[1] +1, ] <- c("Saguinus_G", rep(0.75, 7) )
  row.names(rep.list) <- rep.list$Especie
  #colnames(rep.list) <- c("Especie", "BS.rs")  #### só pra sizeless
  rep.list$Especie <- factor(rep.list$Especie, levels = unique(rep.list$Especie) )
  rep.list$BS.rs <- as.numeric(rep.list$BS.rs)
  rep.list$BS.krz <- as.numeric(rep.list$BS.krz)
  rep.list$pcas <- as.numeric(rep.list$pcas)
  rep.list$cor.mantel <- as.numeric(rep.list$cor.mantel)
  rep.list$cor.krz <- as.numeric(rep.list$cor.krz)
  rep.list$MC.rs <- as.numeric(rep.list$MC.rs)
  rep.list$MC.krz <- as.numeric(rep.list$MC.krz)
  
  str(rep.list)

  sample.size.list <- c(n.size[mask,2], 130, 230)

  safe.copy.mx.compare <- mx.compare
  mx.compare = vector("list", 7)
  mx.compare[1:7] <- NA
  mx.compare[[1]] <- RandomSkewers(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 2] )
  mx.compare[[2]]$correlations <- as.matrix(KrzCor(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 3] ))
  mx.compare[[3]]$correlations <- as.matrix(PCAsimilarity(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 4] ) )
  mx.compare[[4]]$correlations <- as.matrix(MatrixCor(cor.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 5] ) )
  mx.compare[[5]]$correlations <- as.matrix(KrzCor(cov.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 6] ))
  mx.compare[[6]] <- RandomSkewers(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 7] )
  mx.compare[[7]]$correlations <- as.matrix(KrzCor(cov.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 8] ))
  names(mx.compare)[1:7] <-  c('BS.RS', 'BS.KRZ','PCA.s', 'Mantel.Cor', 'KRZ.Cor', 'MC.RS', 'MC.KRZ')
  mx.class<- c('V/CV', 'V/CV','V/CV', 'COR', 'COR','V/CV','V/CV')
  for (i in 1:7)  {mx.compare[[i]]$method <- names(mx.compare)[i]}
  for (i in 1:7)  {mx.compare[[i]]$mx.class <- mx.class[i]}
  
  mat_data <- mx.compare$RS$correlations
  mat_data_raw <-t(mat_data)
  mat_data_raw[lower.tri(mat_data)] <- mx.compare$KRZ$correlations[lower.tri(mx.compare$KRZ$correlations)]
  
  mat_data[lower.tri(mat_data)] <- t(mx.compare$KRZ$correlations)[lower.tri(mx.compare$KRZ$correlations)]
  range.values<- range(mat_data, na.rm = T)
  diag(mat_data) <- sample.size.list
  
  ##################################################################################################
  Combine.Mx.Plot <- function(Mx1, Mx2, diag.info = NA, prob = NULL, titulo = "blé"){
    ## Mx1 = Matriz que vai pegar os valores da diagonal de cima e manter na diagonal de cima (geralmente os results de RS)
    ## Mx2 = Matriz que vai pegar os valores da diagonal de cima e passar pra diagonal de baixo (geralmente os results de KRZ)
    ## prob = Matriz com as probabilidades do result comparação de RS
    
    dimnames(Mx1)[[1]] %<>% gsub("_", ' ',.)
    dimnames(Mx1)[[2]] %<>% gsub("_", ' ',.)
    
    mat_data <- Mx1
    mat_data[lower.tri(mat_data)] <- t(Mx2)[lower.tri(Mx2)]
    diag(mat_data)<- NA
    range.values<- range(mat_data, na.rm = T)
    range.values<- range.values- c(0.01, -0.01)
    diag(mat_data) <- diag.info
    
    sig.mx <- prob >= 0.04
    sig.mx[upper.tri(sig.mx, diag = T)] <- FALSE
    sig.mx [sig.mx == TRUE] <- "white"
    sig.mx [sig.mx == FALSE] <- "black"
    sig.mx <- t(sig.mx)
    melted.sig <- melt(sig.mx)
    names(melted.sig) <- c("Sp1", "Sp2", "prob")
    melted.sig$black
    Prob <- melted.sig$prob
    
    
    if method !=RS 

    mixed.mx = melt(mat_data) 
    mixed.mx.position =  mixed.mx
    mixed.mx.position$value = round( mixed.mx.position$value, 2)
    
    
    myPalette <- colorRampPalette(rev(brewer.pal(5, 'Spectral')), space = 'Lab')(n = 3)
    mixed.mx.cute.plot <- 
      ggplot (mixed.mx.position) +
      geom_tile(aes(x = Var2, y = Var1, fill = value)) +
      scale_fill_gradientn(name = '', colours = myPalette, limits = range.values, na.value = "white") +
      ylab ('') + xlab ('') + labs(title = titulo) + theme(plot.title = element_text(face = "bold", size = 30)) +
      geom_text(aes(x = Var2, y = Var1, label = value), size = 3, color = as.character(Prob) ) +
      scale_y_discrete(limits = rev(levels(mixed.mx.position$Var1))) +
      # scale_x_discrete() +
      theme_minimal() +  
      theme(axis.text.x = element_text(angle = 270, hjust = 0, face = 'italic', size =13),
            axis.text.y = element_text(face = "italic", size =13),
            axis.ticks = element_line(size = 0),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            rect = element_blank(), line   Steppan$Plot$KRZ = element_blank())
    
    return(mixed.mx.cute.plot)
    
  }
  
  ########################## Montando os plots de combinação valores de comparação corrigidos, nao corrigidos e por método

  Combine.Mx.Plot(Mx1 =   mx.compare.sizeless$RS$correlations, Mx2 = t(mx.compare.log$RS$correlations), prob =   mx.compare.sizeless$RS$probabilities, diag.info = sample.size.list, titulo = "Species matrices comparisons")
  Steppan <- list("RS" = mx.compare$BS.RS$correlations - (-6.437 /harm_matrix),
                  "KRZ" = mx.compare$BS.KRZ$correlations - (-6.437 /harm_matrix))
  
    Steppan$Plot$RS <- Combine.Mx.Plot(Mx1 = t(mx.compare$BS.RS$correlations), Mx2 = t(Steppan$RS) , prob = mx.compare$BS.RS$probabilities, diag.info = c(n.size[mask,2], 130, 12), titulo = "Species matrices comparisons via Random Skewers")
  Steppan$Plot$RS
  Steppan$Plot$KRZ <- Combine.Mx.Plot(Mx1 = (mx.compare$BS.KRZ$correlations), Mx2 = t(Steppan$KRZ) , prob = mx.compare$BS.RS$probabilities, diag.info = c(n.size[mask,2], 130, 12), titulo = "Species matrices comparisons via Krzanowski")
  Steppan$Plot$KRZ   
  
  
  plot.mix.raw<- Combine.Mx.Plot(Mx1 = t(mx.compare$BS.RS$correlations), Mx2 = t(mx.compare$BS.KRZ$correlations), prob = mx.compare$BS.RS$probabilities, diag.info = sample.size.list, titulo = "Species matrices comparisons")
  plot.mix.correct.rep.mx<- Combine.Mx.Plot(Mx1 = mx.compare$BS.RS$correlations, Mx2 = mx.compare$BS.KRZ$correlations, diag.info = sample.size.list, titulo = "Matrices comparisons")
  
  plot.b<- Combine.Mx.Plot(Mx1 = mx.compare$PCA.s$correlations, Mx2 = mx.compare$KRZ$correlations, diag.info = sample.size.list, titulo = "Covariance matrices comparison values via KRZ and PCA Similarity ")
  plot.c<- Combine.Mx.Plot(Mx1 = mx.compare$Mantel.Cor$correlations, Mx2 = mx.compare$KRZ.Cor$correlations, diag.info = sample.size.list, titulo = "B. Correlation matrices comparison values via KRZ and Mantel ")
  
  plot_grid(plot.a, plot.c, nrow = 2, vjust = 0.1)
  
  plot.a.log <- Combine.Mx.Plot(Mx1 = mx.compare.log$RS$correlations, Mx2 = mx.compare.log$KRZ$correlations, diag.info = n.size[mask,2], titulo = "A. Covariance matrices comparison values via KRZ and RS ")
  plot.b.log <- Combine.Mx.Plot(Mx1 = mx.compare$PCA.s$correlations, Mx2 = mx.compare$KRZ$correlations, diag.info = sample.size.list, titulo = "Covariance matrices comparison values via KRZ and PCA Similarity ")
  plot.c.log <- Combine.Mx.Plot(Mx1 = mx.compare$Mantel.Cor$correlations, Mx2 = mx.compare$KRZ.Cor$correlations, diag.info = sample.size.list, titulo = "B. Correlation matrices comparison values via KRZ and Mantel ")
  
  plot_grid(plot.a, plot.c, nrow = 2, vjust = 0.1)
  
  plot_grid(plot.a, plot.a.log, nrow = 2, vjust = 0.1)
  
   
 mat_data2 <- array (0, dim (mat_data))
 mat_data2 [lower.tri (mat_data2)] <- mat_data [lower.tri (mat_data)]
 mat_data2 <- mat_data2 + t(mat_data2)
 diag(mat_data2) <- 1
 dimnames (mat_data2) <- dimnames (mat_data)
  
 eigen1matdata <- Re(eigen (mat_data2) $ vectors [,1])
 
 mat_data2 <- mat_data2 [order(eigen1matdata, decreasing = T), order(eigen1matdata, decreasing = T)]
 diag(mat_data2) <- sample.size.list

  log.cov.mx <- current.data %>% llply(function(x) x$matrix$cov.log)
  log.cor.mx <- llply(log.cov.mx[mask], cov2cor)
  mx.compare.log = vector("list", 5)
  mx.compare.log[1:5] <- NA
  mx.compare.log[[1]] <- RandomSkewers(cov.x= log.cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 2] )
  mx.compare.log[[2]]$correlations <- as.matrix(KrzCor(cov.x= log.cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 3] ))
  mx.compare.log[[3]]$correlations <- as.matrix(PCAsimilarity(cov.x= log.cov.mx[mask], num.vectors = 1000, repeat.vector = mx.rep[mask, 4] ) )
  mx.compare.log[[4]]$correlations <- as.matrix(MatrixCor(cor.x= log.cor.mx, correlation = TRUE, num.vectors = 1000, repeat.vector = mx.rep[mask, 5] ) )
  mx.compare.log[[5]]$correlations <- as.matrix(KrzCor(cov.x= log.cor.mx, correlation = TRUE, num.vectors = 1000, repeat.vector = mx.rep[mask, 6] ))
  names(mx.compare.log)[1:5] <-  c('RS', 'KRZ','PCA.s', 'Mantel', 'KRZ')
  mx.class<- c('V/CV', 'V/CV','V/CV', 'COR', 'COR')
  for (i in 1:5)  {mx.compare.log[[i]]$method <- names(mx.compare.log)[i]}
  for (i in 1:5)  {mx.compare.log[[i]]$mx.class <- mx.class[i]}

  
  
##################################
#Função para plotar matrizes 
#com respectivos valores seguindo 
# palleta de padrao de cor
###################################
myPalette <- colorRampPalette(c("yellow", "white", "purple"))(n = 100)
myPalette <- "RdBu"  
plot.matrix<- function(mx = NULL, sample.size = NULL, brewer = myPalette, show.values = TRUE)
{ require (plotrix)
  library("RColorBrewer")
  paleta.of.choice  <- rev(brewer.pal(11, brewer))
  colfunc <- colorRampPalette(paleta.of.choice)
  paleta  <- colfunc(90)
  nCores  <- length(paleta)
  intervalMarks  <- seq(from = 0.1, to = 1, length.out = nCores)
  
  mx.cor <- mx
  mx.dimentions <- dim(mx.cor)[1]
  
  dados<- as.vector(mx.cor)
  intervals  <- findInterval(dados , intervalMarks, rightmost.closed = T, all.inside = T)
  cores  <- paleta[intervals]
  
  color2D.matplot(x = mx.cor, axes = F, cellcol = cores, show.values= 2, vcex= 0.6, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  #title(main = paste(mx$mx.class, "Matrices compared by", mx$method, "method"))
  title(main = "Matrix comparison via KRZ (lower) and RS (upper)")
  axis(1, 1:mx.dimentions, paste("(n=", sample.size,")",1:length(rownames(mx.cor)), sep = " "), las = 2, cex.axis = 0.9, tick = FALSE, line = 0)
  #axis(2, mx.dimentions:1, length(rownames(mx.cor)):1, las = 1, cex.axis = 1, tick = FALSE, line = 0)
  axis(2, mx.dimentions:1, paste(rownames(mx.cor), 1:length(rownames(mx.cor))), las = 1, cex.axis = 0.9, tick = FALSE, line = 0)

}
par(mar = c(7,17,3,2))
par(mfrow = c(3,1))
  lapply(mx.compare, plot.matrix)
  lapply(mx.compare.log, plot.matrix)
  
  
  plot.matrix(mx = mat_data, sample.size = sample.size.list, show.values = TRUE)
  #Calculando a similaridade média das matrizes.
mean.sim<-function(x) {
    x[upper.tri(x)]<-t(x)[upper.tri(x)]
    diag(x)<-NA
    x<-rowMeans(x,na.rm=TRUE)
    return(x)
}

##########Plotando os valores médios de comparaçao por metodo por matriz
current.mx.list <- mx.compare
current.mx.list <- mx.compare.log

mean.comp.values.cov<- current.mx.list[c(1:3)] %>% llply (function (x) x$correlations ) %>% llply ( function (x) x[-c(2,39), -c(2,39)]) %>% ldply (function (x) mean.sim(x) )
mean.comp.values.cor<- current.mx.list[c(4:5)] %>% llply (function (x) x$correlations ) %>% llply ( function (x) x[-c(2,39), -c(2,39)]) %>% ldply (function (x) mean.sim(x) )
mean.comp.values <- mean.comp.values.cov
names(mean.comp.values)[1] <- "method"

method.mean <- data.frame("method"  = mean.comp.values[,1], "mean"=apply(mean.comp.values[,-1], 1, mean))
Mean.Mx.Plots <- mean.comp.values %>% 
  gather(key=.sp, value = value, 2:45) %>%
  ggplot( .,aes(x= method, y = value, shape = method), varwidth = T) +
  geom_violin(aes(label = .sp, color = method, shape = method), alpha = 0.5) +
  geom_text( aes(label = .sp), size =2, vjust = 1, alpha = 0.4)  +
  geom_jitter(aes(shape = method, color = method)) +
  scale_shape( guide = "none", name = "Mean value by method") +
  ggtitle("Mean value of comparisson") +
  theme_bw() +
  theme(plot.title = element_text(lineheight=.8, face="bold"))


Mean.Mx.Plots <- mean.comp.values %>% 
  gather(key=.sp, value = value, 2:45) %>% 
  mutate(sample = rep(sample.size.list, each = 3)) %>%
  filter(sample < 120) %>%
  filter(sample > 20) %>%
  filter(method != "PCA.s") %>%
  ggplot( .,aes(x = 1, y = value, shape = method), varwidth = T) +
  geom_violin(aes(label = .sp, size = sample, shape = method)) +
  geom_text( aes(label = .sp), size =7, vjust = 1, alpha = 0.4)  +
  geom_jitter(aes(shape = method, size = sample, color = interaction(sample, method)) ) +
  ggtitle("Comparação de matrizes de espécie por método") +
  facet_grid(~method) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.ticks = element_line(size = 0),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 30),
        plot.title = element_text(lineheight=.8, size = 20, face="bold"),
        legend.position = c(.06, .15)) +
  labs (y = "Valor médio de comparação") +
  guides(colour = FALSE,
         shape = FALSE,
         size = guide_legend(title = "Tamanho\nAmostral", override.aes = list(alpha = 0.5) ) )



Mean.Mx.Plots + geom_point(aes (x=method, y=mean, color = method, shape = method), data = method.mean, size = 6 )


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
current.data <- sp.main.data
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
PC1 <- current.data[mask] %>% ldply(function(x) eigen(x$matrix$cov)$vectors[,1])
names.sp <-  PC1[,1]
PC1<- as.matrix(PC1[,-1])
dimnames(PC1)[[1]] <- names.sp
# se estiver correlacionado negativamente com o vetor isometrico, multiplicar por -1
ed.names <- names(current.data[mask][[1]]$ed)
Iso.Compare <- function(x) {
  isometrico<- rep( (1/sqrt(39)) , 39)
  corr = rep(NA, length(x) )
  re.oriented = x
    for (i in 1:length(x)) {
    corr[i]<- x[[i]] %*% isometrico
      
      if(corr[i] <0) {re.oriented[[i]] <- as.numeric(as.matrix(x[[i]]) %*% -1)}  else x[[i]]
    }
    
  names (corr) <- paste ("PC", 1:length(x), sep = "")
  names (re.oriented) <- paste ("PC", 1:length(x), sep = "") 
  
  return(results= list("corr" = corr, "re.oriented" = re.oriented)) }
PCs1to4<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )
PCs1to4.log<- current.data[mask] %>% llply(function(x) as.data.frame(eigen(x$matrix$cov.log)$vectors[,1:4]) ) %>% llply(function(x) as.list(x) )

Iso.Compare.reoriented <- llply(PCs1to4.log, Iso.Compare) %>% ldply(function(x) as.data.frame(x$re.oriented)) 

names(Iso.Compare.reoriented) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.reoriented$.ed <- ed.names
Iso.Compare.reoriented %>% gather(key=".pcScore", value=value, 2:5 )
Iso.Compare.reoriented$.sp <- factor(Iso.Compare.reoriented$.sp, levels = unique(Iso.Compare.reoriented$.sp)) 
Iso.Compare.reoriented$.ed <- factor(Iso.Compare.reoriented$.ed, levels = rev(ed.names) ) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
Iso.Compare.reoriented %>% gather(key=".pcScore", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = .pcScore, y = .ed, fill = value)) +
  facet_wrap(~.sp) +
  theme_bw() +
  scale_fill_gradientn(name = 'PC.Score', colours = myPalette(1000)) +
  ylab ('PC.Scores') + xlab ('') + labs(title = "First 4 PC's scores") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(size= 7),
        strip.text= element_text(size=7)) +
  coord_flip()

Iso.Compare.cor.iso <- llply(PCs1to4, Iso.Compare) %>% ldply(function(x) x$corr) 
names(Iso.Compare.cor.iso) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
Iso.Compare.cor.iso$.sp <- factor(Iso.Compare.cor.iso$.sp, levels = rev(unique(Iso.Compare.cor.iso$.sp))) 
Iso.Compare.cor.iso %>% gather(key="Isometric.Correlation", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = Isometric.Correlation, y = .sp, fill = abs(value) ) ) +
  theme_bw() +
  scale_fill_gradientn(name = 'Isometric.Correlation', colours = myPalette(1000)) +
  ylab ('') + xlab ('') + labs(title = "Matrix of log values correlation with isometric vector") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10),
        axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(size= 10)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))

PCs1to4<- current.data[mask] %>% ldply(function(x) as.data.frame(eigen(x$matrix$cov)$vectors[,1:4]) )  
names(PCs1to4) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
PCs1to4$.ed <- ed.names
PCs1to4 %>% gather(key=".pcScore", value=value, 2:5 )
PCs1to4$.sp <- factor(PCs1to4$.sp, levels = unique(PCs1to4$.sp)) 
PCs1to4$.ed <- factor(PCs1to4$.ed, levels = rev(ed.names) ) 
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
PCs1to4 %>% gather(key=".pcScore", value=value, 2:5 ) %>%
  ggplot (.) +
  geom_tile(aes(x = .sp, y = .ed, fill = value)) +
  facet_wrap(~.pcScore) +
  theme_bw() +
  scale_fill_gradientn(name = 'PC.Score', colours = myPalette(1000)) +
  ylab ('PC.Scores') + xlab ('') + labs(title = "First 4 PC's scores") +
  #scale_y_discrete(limits = levels(PCs1to4$.pcScore))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0))


################## Grafico que nao serve pra nada, nao sai nenum padrao
####################### nao faz sentido olhar só pros valores das gm, tenho que pensar no quanto de variação que tem em tamanho.

gm.mean <- sp.main.data[mask] %>% ldply(function(x) x$gm.mean)
sd.gm <- sp.main.data[mask] %>% ldply(function(x) sd(x$gm.ind, na.rm = T) )
PcPercent.1to4<- cov.mx[mask] %>% ldply(function(x) eigen(x)$values[1:4]/sum(eigen(x)$values) )
names(PcPercent.1to4) <- paste ("percent", names(PcPercent.1to4), sep = "")

names(sd.gm)[2] <- "sd.gm"

data.frame ('sd.gm' = sd.gm[2],  'gm' = gm.mean[,2], 'Sample' = n.size [mask, -1], Iso.Compare.cor.iso, PcPercent.1to4[2:5]) %>%
  gather(key = '.pc', value=value, c(5:8)) %>%
  ggplot (.) +
  geom_point(aes (x = sd.gm^2, y = percentPC1, color = Sample ), size = 3) +
  scale_colour_gradientn(colours = terrain.colors(7) ) +
  geom_text(aes (x =  sd.gm^2, y = percentPC1, label = .sp), alpha = 0.4, size =4) +
  #scale_color_brewer(name = 'PC.Score' , type = "div"  , palette =  myPalette(1000)) +
  ylab ('Absolute value of correlation') + xlab ('sd of traits geometric mean (Skull size)') + labs(title = "First 4 PC's correlation with Isometric vector") +
  facet_wrap(~ .pc, scales = "fixed") +
  theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold"),
        axis.text.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=17),
        axis.title.x = element_text(size=17),
        strip.text= element_text(size=15) )

#######################################################################################################################################
########################################################################################################################################
########################################################################################################################################


########### quantidade de variação espalhada nos primeiros 4 PC
  PcPercent.1to4<- cov.mx[mask] %>% ldply(function(x) eigen(x)$values[1:4]/sum(eigen(x)$values) )
  names(PcPercent.1to4) <- c(".sp", "PC1", "PC2", "PC3", "PC4")
  PcPercent.1to4$.sp <- factor(PcPercent.1to4$.sp, levels = unique(PcPercent.1to4$.sp)) 
  library(scales)
  PcPercent.1to4 %>%
  gather(key="Percent_var_PC", value=value, 2:5 ) %>%
    ggplot( ., aes(x= Percent_var_PC, y = value, color = Percent_var_PC, label = .sp), varwidth = T) +
    scale_y_continuous(labels=percent) +
    geom_text(size = 2, vjust = 1, alpha = 0.4 )  +
    theme_bw() +
    geom_jitter() +
    geom_violin(alpha = 0) +
    ggtitle("% of variance on first 4 PC by specie") +
    theme(plot.title = element_text(lineheight=.8, face="bold"), axis.title.x = element_blank()) +
    theme(legend.position="none")


  PhyloW(tree = tree.indridae, 
         tip.data = cov.indridae, 
         tip.sample.size = n.size.indridae )
