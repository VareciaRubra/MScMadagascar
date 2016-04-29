
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
  #cov.list <- cov.sizeless[mask]
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
  colnames(rep.list) <- c("Especie", "BS.rs", "BS.krz", "pcas", "cor.mantel", "cor.krz", "MC.rs", "MC.krz") 
  colnames(rep.list)
  rep.list$Especie <- factor(rep.list$Especie, levels = unique(rep.list$Especie) )
  rep.list$BS.rs <- as.numeric(rep.list$BS.rs)
  rep.list$BS.krz <- as.numeric(rep.list$BS.krz)
  rep.list$pcas <- as.numeric(rep.list$pcas)
  rep.list$cor.mantel <- as.numeric(rep.list$cor.mantel)
  rep.list$cor.krz <- as.numeric(rep.list$cor.krz)
  rep.list$MC.rs <- as.numeric(rep.list$MC.rs)
  rep.list$MC.krz <- as.numeric(rep.list$MC.krz)
  
  str(rep.list)

  sample.size.list <- c(n.size[mask,2], 130, 36)

  safe.copy.mx.compare <- mx.compare
  mx.compare = vector("list", 7)
  mx.compare[1:7] <- NA
  mx.compare[[1]] <- RandomSkewers(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 2] )
  mx.compare[[2]]$correlations <- as.matrix(KrzCor(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 3] ))
  mx.compare[[3]]$correlations <- as.matrix(PCAsimilarity(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 4] ) )
  mx.compare[[4]]$correlations <- as.matrix(MatrixCor(cor.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 5] ) )
  mx.compare[[5]]$correlations <- as.matrix(KrzCor(cov.x= cor.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 6] ))
  mx.compare[[6]] <- RandomSkewers(cov.x= cov.list, num.vectors = 1000, repeat.vector = rep.list[, 7] )
  mx.compare[[7]]$correlations <- as.matrix(KrzCor(cov.x= cov.list, correlation = TRUE, num.vectors = 1000, repeat.vector = rep.list[, 8] ))
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
  Combine.Mx.Plot <- function(Mx1, Mx2, diag.info = NA, prob = NULL, titulo = "blé", method = NULL){
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
    ############# se for plotar duas mx de RS coment essa linha de cima e gera essa mx mista por fora
    # combined.prob <- MX1$probabilities
    # combined.prob [upper.tri(combined.prob, diag = F)] <- t(MX2$probabilities)[upper.tri(MX2$probabilities, diag = F) ]
    
    sig.mx [sig.mx == TRUE] <- "white"
    sig.mx [sig.mx == FALSE] <- "black"
    sig.mx <- t(sig.mx)
    melted.sig <- melt(sig.mx)
    names(melted.sig) <- c("Sp1", "Sp2", "prob")
    melted.sig$black <- melted.sig$prob
    melted.sig$black [melted.sig$black == "white"] <- "black"

    Prob <-  if (method == "Random Skewers") melted.sig$prob else melted.sig$black

    mixed.mx = melt(mat_data) 
    mixed.mx.position =  mixed.mx
    mixed.mx.position$value = round( mixed.mx.position$value, 2)
    mixed.mx.position$texto <-  mixed.mx.position$value %>% gsub("0.", ".", .)
   
    myPalette <- colorRampPalette(rev(brewer.pal(5, 'Spectral')), space = 'Lab')(n = 40)
    mixed.mx.cute.plot <- 
      ggplot (mixed.mx.position) +
      geom_tile(aes(x = Var2, y = Var1, fill = value), alpha = 0.9) +
      scale_fill_gradientn(name = "Value of \ncomparison", colours = myPalette, limits = range.values, na.value = "white") +
      geom_text(aes(x = Var2, y = Var1, label = texto), size = 4, color = as.character(Prob)) +
      scale_y_discrete(limits = rev(levels(mixed.mx.position$Var1))) +
      ylab ('') + xlab ('') + labs(title = paste(titulo, method) ) + 
      theme_minimal() +  
      theme(plot.title = element_text(face = "bold", size = 20),
            axis.text.x = element_text(angle = 270, hjust = 0, face = 'italic', size =13),
            axis.text.y = element_text(face = "italic", size =13),
            axis.ticks = element_line(size = 0),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10),
            rect = element_blank(), 
            line = element_blank())
    
    return(mixed.mx.cute.plot)
    
  }
  
  ########################## Montando os plots de combinação valores de comparação corrigidos, nao corrigidos e por método

  Combine.Mx.Plot(Mx1 =   mx.compare.sizeless$RS$correlations,
                  Mx2 = t(mx.compare.log$RS$correlations),
                  prob =   mx.compare.sizeless$RS$probabilities,
                  diag.info = sample.size.list[-c(43:44)],
                  titulo = "Species matrices comparisons",
                  method = "RS")

  Steppan$Plot$BS.RS <- 
    Combine.Mx.Plot(Mx1 = t(mx.compare$BS.RS$correlations), 
                    Mx2 = Steppan$BS.RS, 
                    prob = mx.compare$BS.RS$probabilities, 
                    diag.info = c(n.size[mask,2], 130, 12), 
                    method = "Random Skewers",
                    titulo = "Species matrices comparisons via " )

  Steppan$Plot$BS.KRZ <- 
   Combine.Mx.Plot(Mx1 = t(mx.compare$BS.KRZ$correlations), 
                  Mx2 = Steppan$BS.KRZ, 
                  prob =  mx.compare$BS.RS$probabilities, 
                  diag.info = c(n.size[mask,2], 130, 12), 
                  method = "Krzanowski",
                  titulo = "Species matrices comparisons via " )
  
  Steppan$Plot$rawRaw.mc <- 
#     Combine.Mx.Plot(Mx1 = t(mx.compare$MC.RS$correlations), 
#                     Mx2 = t(mx.compare$MC.KRZ$correlations), 
#                     prob =  mx.compare$MC.RS$probabilities, 
#                     diag.info = c(n.size[mask,2], 130, 12), 
#                     method = "Random Skewers",
#                     titulo = "Raw values of comparison via Krzanowski and") # tamanho salvo 15 x 12
  
  Steppan$Plot$correct <- 
#    Combine.Mx.Plot(Mx1 = Steppan$BS.RS, 
#                     Mx2 = Steppan$BS.KRZ, 
#                     prob =  mx.compare$BS.RS$probabilities, 
#                     diag.info = c(n.size[mask,2], 130, 12), 
#                     method = "Random Skewers",
#                     titulo = "Steppan's correction values of comparison via Krzanowski and") # tamanho salvo 15 x 12
Steppan$Plot$Corrected.mc <-     
    Combine.Mx.Plot(Mx1 = mx.compare$MC.RS$correlations, 
                    Mx2 = mx.compare$MC.KRZ$correlations, 
                    prob = mx.compare$MC.RS$probabilities,
                    diag.info = c(n.size[mask,2], 130, 36),
                    method = "Random Skewers",
                    titulo = "Corrected repetabilities Species matrices comparisons via Krzanowski and" )

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
    
Steppan$Plot$SizelesAndLog <- 
    Combine.Mx.Plot(Mx1 = t(mx.compare.log$RS$correlations), 
                    Mx2 = t(mx.compare.log$KRZ$correlations), 
                    prob = mx.compare.log$RS$probabilities, 
                    diag.info = c(n.size[mask,2]), 
                    method = "Random Skewers",
                    titulo = "Species matrices comparisons via " )

combined.prob <- mx.compare.log$RS$probabilities
combined.prob [upper.tri(combined.prob, diag = F)] <- t(mx.compare.sizeless$RS$probabilities)[upper.tri(mx.compare.sizeless$RS$probabilities, diag = F) ]

Log.cov.mx <- sp.main.data[mask] %>% llply(function(x) x$matrix$cov.log )
Sizeless.cov.mx <- sp.main.data[mask] %>% llply(function(x) x$matrix$cov.sizeless )

ChechSizeMatters <- sp.main.data[mask] %>% 
  llply(function(x) x$rs.size.comparisson$correlations ) %>% 
  ldply (function(x) x[lower.tri(x, diag= F)] ) 
names (ChechSizeMatters) <- c(".id", "ED x GM", "ED x log(ED)", "GM x log(ED)")

ChechSizeMatters<- list("comparison" = ChechSizeMatters, 
                        "Plot" = "Plot missing")

ChechSizeMatters$Plot <- Combine.Mx.Plot(Mx1 = t(mx.compare.log$RS$correlations), 
                Mx2 = t(mx.compare.sizeless$RS$correlations), 
                prob = t(combined.prob) , 
                diag.info = ChechSizeMatters$comparison[[2]], 
                method = "Random Skewers",
                titulo = "ed/GM and Log(ed) Matrices comparisons via " )

ChechSizeMatters$comparison$.id %<>% gsub("_", " ", .)
ChechSizeMatters$comparison$.id <- factor(ChechSizeMatters$comparison$.id, unique (ChechSizeMatters$comparison$.id)[42:1])

colore <- c("#FF9326",  "#9999FF", "#D92121", "#21D921", "#FFFF4D", "#2121D9")

ChechSizeMatters$comparison %>% melt %>%
  ggplot(aes(y = value, x = .id, color = variable)) + 
  geom_linerange(aes(ymin =0.35, ymax = value), size = 4, alpha = 0.8) +
  #scale_colour_brewer(palette = "Set1", direction = 1) +
  scale_colour_manual (values=colore) +
  theme_bw() +
  coord_flip() +
  facet_wrap(~variable) + xlab("") + ylab("") +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(face =  "italic", size =10),
        axis.text.x = element_text(size =13),
        legend.position = "none")


 mat_data2 <- array (0, dim (mat_data))
 mat_data2 [lower.tri (mat_data2)] <- mat_data [lower.tri (mat_data)]
 mat_data2 <- mat_data2 + t(mat_data2)
 diag(mat_data2) <- 1
 dimnames (mat_data2) <- dimnames (mat_data)
  
 eigen1matdata <- Re(eigen (mat_data2) $ vectors [,1])
 
 mat_data2 <- mat_data2 [order(eigen1matdata, decreasing = T), order(eigen1matdata, decreasing = T)]
 diag(mat_data2) <- sample.size.list


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
# Eigen Tensors decomposition
mask.n.size <- n.size[,2]>40
cov.list.39 <- cov.mx[mask.n.size]
cov.list.39$AncestrAll <-  W.matrix
cov.list.39$Saguinus_P.cov <- Saguinus_P.cov
# cov.list.39$Saguinus_G.cov <- Saguinus_G.cov

eigenLemur <- EigenTensorDecomposition(cov.list.39 [1:(length(cov.list.39) -1)])
dimnames(eigenLemur$projection)$X1 <-  names(cov.list.39)[-length(cov.list.39)]
eigenLemur.SaguinusP <- ProjectMatrix(Saguinus_P.cov,eigenLemur)
eigenLemur.SaguinusG <- ProjectMatrix(Saguinus_G.cov,eigenLemur)

eigenLS <- rbind (eigenLemur$projection, eigenLemur.SaguinusP, eigenLemur.SaguinusG)
plot(eigenLS[, 1], eigenLS[, 2] )
text(eigenLS[, 1], eigenLS[, 2] , labels = rownames(eigenLS))

plot(eigenLemur$projection[,1], eigenLemur$projection[,2])
text(eigenLemur$projection[,1], eigenLemur$projection[,2], labels = names(cov.list.39)[-20])
eigenLemur.Saguinus

