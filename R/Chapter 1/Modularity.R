Modular.hyp <- Aux$def.hyp[c(2:20, 1, 21:39),] # pegando as hipoteses modulares do arquivo auxiliat do Pato
Modular.hyp[21,] <- c(0,0,1,0,0,1,1,1) # mas ele tem logCS e nao tem PT.TSP, enfiemolhos!
rownames(Modular.hyp)[20] <- "PT.TSP"

Modulemurs <- vector("list")
Modulemurs$test.modularity <- sp.main.data[mask.no.na.cov] %>% 
  llply (function (x) TestModularity(cor.matrix = x$matrix$cor, modularity.hipot = Modular.hyp , permutations = 10000, MHI = F), .progress = "text")

Modulemurs$Probability<- Modulemurs$test.modularity %>% ldply(function (x) x$Probability)
names(Modulemurs$Probability) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
Modulemurs$Probability$sp <- factor (Modulemurs$Probability$sp, levels = rev(unique(Modulemurs$Probability$sp)))
Modulemurs$Probability$sp %<>% gsub("_", " ",.)
Modulemurs$Probability$sp <- factor (Modulemurs$Probability$sp, levels = rev(unique(Modulemurs$Probability$sp)))

Modulemurs$Probs<- Modulemurs$test.modularity %>% ldply(function (x) x$Probability)
names(Modulemurs$Probs) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
rownames(Modulemurs$Probs) <- Modulemurs$Probs[,1]
Modulemurs$Probs[Modulemurs$Probs <= 0.05] <- "red"
Modulemurs$Probs[Modulemurs$Probs != "red"] <- "darkgrey"
Modulemurs$Probs[,1] <- rownames(Modulemurs$Probs)

Modulemurs$Probability %>% gather(key = "Probability", value = probability, 2:10) %>%
ggplot (.) +
  geom_tile(aes(x = Probability, y = sp, fill = as.numeric(probability<= 0.05) ), alpha = 0.9,  color = "grey") +
  #geom_point(aes(x =Probability, y = sp, color = as.numeric(probability<= 0.05) ) ) +
  theme_bw() +
  scale_fill_gradient(name = 'Probability', low = "transparent", high = "red") +
  #scale_color_gradient(name = 'Probability', low = "transparent", high = "red" ) +
  ylab ('') + xlab ('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 8),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        legend.position = "none",
        rect = element_blank(), 
        line = element_blank())

Modulemurs$AVG.ratio<- Modulemurs$test.modularity %>% ldply(function (x) x$`AVG Ratio`)
names(Modulemurs$AVG.ratio) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
Modulemurs$AVG.ratio$sp <- factor (Modulemurs$AVG.ratio$sp, levels = rev(unique(Modulemurs$AVG.ratio$sp)))
Modulemurs$AVG.ratio$sp %<>% gsub("_", " ",.)
Modulemurs$AVG.ratio$sp <- factor (Modulemurs$AVG.ratio$sp, levels = rev(unique(Modulemurs$AVG.ratio$sp)))
rownames(Modulemurs$AVG.ratio) <- Modulemurs$AVG.ratio[,1]
Modulemurs$AVG.ratio[,1] <- rownames(Modulemurs$AVG.ratio)

Modulemurs$AVG.ratio %>% gather(key = "Region", value = AVG.ratio, 2:10) %>% 
  ggplot (.) +
  geom_tile(aes(x = Region, y = sp, fill = abs(AVG.ratio)) , alpha = 0.9,  color = "grey") +
  #geom_point(aes(x =AVG.ratio, y = sp, color = as.numeric(AVG.ratio<= 0.05) ) ) +
  theme_bw() +
  scale_fill_gradient2(name = 'AVG ratio', low = "transparent", mid = "#f46d43", high ="#d73027", midpoint = 3) +
  ylab ('') + xlab ('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 8),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        #legend.position = "none",
        rect = element_blank(), 
        line = element_blank())

Modulemurs$AVG.ratio<- Modulemurs$test.modularity %>% ldply(function (x) x$`AVG Ratio`)
names(Modulemurs$AVG.ratio) <- c("sp", Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
Modulemurs$AVG.ratio$sp <- factor (Modulemurs$AVG.ratio$sp, levels = rev(unique(Modulemurs$AVG.ratio$sp)))
Modulemurs$AVG.ratio$sp %<>% gsub("_", " ",.)
Modulemurs$AVG.ratio$sp <- factor (Modulemurs$AVG.ratio$sp, levels = rev(unique(Modulemurs$AVG.ratio$sp)))
rownames(Modulemurs$AVG.ratio) <- Modulemurs$AVG.ratio[,1]
Modulemurs$AVG.ratio[,1] <- rownames(Modulemurs$AVG.ratio)

Modulemurs$AVG.ratio %>% gather(key = "Region", value = AVG.ratio, 2:10) %>% 
  ggplot (.) +
  geom_tile(aes(x = Region, y = sp, fill = abs(AVG.ratio)) , alpha = 0.9,  color = "grey") +
  #geom_point(aes(x =AVG.ratio, y = sp, color = as.numeric(AVG.ratio<= 0.05) ) ) +
  theme_bw() +
  scale_fill_gradient2(name = 'AVG ratio', low = "transparent", mid = "#f46d43", high ="#d73027", midpoint = 3) +
  ylab ('') + xlab ('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1, size = 8),
        axis.ticks = element_line(size = 0),
        axis.text.y = element_text(face = "italic", size= 10),
        #legend.position = "none",
        rect = element_blank(), 
        line = element_blank())

Modulemurs$cor.pc1.flex <- cbind.data.frame("cor.Flex" = cor(abs(Modulemurs$AVG.ratio[, 2:10]), Flex$V1), 
                                            "cor.PC1" = cor(abs(Modulemurs$AVG.ratio[, 2:10]), PC1.percent$V1))

Flex <- Variae %>% ldply( function(x) mean(x$flex.dist))
PC1.percent <- Variae %>% ldply( function(x) x$intervalo.mc.pc$observed[1] )
AVG.r <- Modulemurs$AVG.ratio %>% gather(key = "Region", value = AVG.ratio, 2:10) 
AVG.r$Flex <- rep(Flex$V1, 9) 
AVG.r$PC1.percent <- rep(PC1.percent$V1, 9) 
AVG.r$cor.flex <- rep(Modulemurs$cor.pc1.flex$cor.Flex, each = length(unique(AVG.r$sp)) )
AVG.r$cor.PC1<- rep(Modulemurs$cor.pc1.flex$cor.PC1, each = length(unique(AVG.r$sp)) )
Modulemurs$Plot$AVG.Flex <- AVG.r %>% tbl_df %>%
  ggplot (aes(group = Region)) +
  geom_smooth(aes(y = AVG.ratio, x = Flex), method = "lm") +
  geom_text(aes(y = 3, x = 0.27 , label = round(cor.flex,2), group = Region)) +
  #geom_text(aes(y = AVG.ratio, x = Flex, label = sp), size = 2) + 
  geom_point(aes(y = AVG.ratio, x = Flex), size = 1, color = "darkgrey") + 
  facet_wrap(~Region, scale = "free")
Modulemurs$Plot$AVG.PC1 <- AVG.r %>% tbl_df %>%
  ggplot (aes(group = Region )) +
  geom_smooth(aes(y = AVG.ratio, x = PC1.percent), method = "lm") +
  geom_text(aes(y = 3, x = 0.65 , label = round(cor.PC1, 2), group = Region)) +
  #geom_text(aes(y = AVG.ratio, x = PC1.percent, label = sp), size = 2) + 
  geom_point(aes(y = AVG.ratio, x = PC1.percent), size = 1, color = "darkgrey") + 
  facet_wrap(~Region, scale = "free")

plot_grid(Modulemurs$Plot$AVG.Flex, Modulemurs$Plot$AVG.PC1, ncol = 2) #########################################################################################


  
class(Modulemurs$test.modularity$Tarsius_bancanus$hypothesis)
Variae$Tarsius_bancanus


DistModular <- function (x, simulations = 1000, modularity.hipot = Modular.hyp)
  {
  hypothesis <- c(dimnames(Modular.hyp)[[2]], "Full.integration")
  Rsquared <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]] +1)
  colnames(Rsquared) <- hypothesis
  AVG.plus <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]] + 1)
  colnames(AVG.plus) <- hypothesis
  AVG.minus <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]]+ 1)
  colnames(AVG.minus) <- hypothesis
  Probability <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]]+ 1)
  colnames(Probability) <- hypothesis
  AVG.rat <- matrix(data = NA, nrow = simulations, ncol = dim(Modular.hyp)[[2]]+ 1)
  colnames(AVG.rat) <- hypothesis
 
  
  for (i in 1:simulations)  {
    Mx.Cor <- var (rmvnorm (x$sample.size, sigma = x$matrix$cor, method = 'svd') )
    modular <- TestModularity(cor.matrix = Mx.Cor, modularity.hipot = Modular.hyp , permutations = 200, MHI = F) 
    
    Rsquared [i, ] <- modular$Rsquared
    AVG.plus [i, ] <- modular$`AVG+`
    AVG.minus [i, ] <- modular$`AVG-`
    Probability [i, ] <-  modular$Probability
    AVG.rat [i, ] <-  modular$`AVG Ratio`
  }
  plotaisso <- cbind(rep (unique(x$info.raw$Especie), simulations),
  melt(Probability, value.name = "Probability")[2:3],
  melt(AVG.plus, value.name = "AVG+")[3],
  melt(AVG.minus, value.name = "AVG-")[3] )
  names(plotaisso)[1:2] <- c("sp", "hyp")
  
  who <- unique(x$info.raw$Especie)
  Probs <- cbind(rep (unique(x$info.raw$Especie), simulations), melt(Probability, value.name = "Probability")[2:3])
  names(Probs)[1:2] <- c("sp", "hyp")
  AVG.r <- cbind(rep (unique(x$info.raw$Especie), simulations), melt(AVG.rat, value.name = "AVG.rat")[2:3])
  names(AVG.r)[1:2] <- c("sp", "hyp")
  AVG.r.plot <- AVG.r %>% 
    ggplot(., aes (x = hyp, y = AVG.rat, group = hyp)) +
    geom_violin() +
    ggtitle(who) 
  
  plot <- plotaisso %>% melt %>%
    ggplot(., aes(x = variable, y = value, group = interaction(hyp, variable))) +
    geom_boxplot(aes(group = variable)) +
    ggtitle(who) +
    facet_wrap(~hyp, ncol = 3)
  
  return(list ("plotaisso" = plotaisso,
               "plot" = plot,
               "probability" =  Probs, 
               "AVG.r" =  AVG.r,
               "AVG.r.plot" = AVG.r.plot
               ) )
}

temp <- DistModular(x = sp.main.data$Varecia_variegata, simulations = 1000, modularity.hipot =  Modular.hyp )
boxplot(temp$probability$Probability)
temp$plot
temp$AVG.r %>% ggplot(., aes (x = hyp, y = AVG.rat, group = hyp)) +
  geom_violin() +
  geom_boxplot() +
  geom_point(aes (x = hyp, y = mean(AVG.rat), group = hyp) )
  geom_errorbar(aes(x = hyp, ymin = min(AVG.r), ymax = max(AVG.r), group = hyp ))+
  ggtitle(who) 

Modulemurs$test.modularity.dist <- sp.main.data[mask.no.na.cov] %>% llply(function (x) DistModular(x = x, simulations = 100, modularity.hipot =  Modular.hyp) , .progress = "text")
save(Modulemurs, file = "Data/Modularity.RData")


Modulemurs$test.modularity.dist %>% ldply(function(x) x$probability) %>% melt %>% 
  ggplot(., aes(x = variable, y = value, group = interaction(variable, hyp, sp))) +
  geom_violin() +
  facet_wrap(~hyp, ncol = 3)


Modulemurs$test.modularity.dist %>% ldply(function(x) x$plotaisso[, 3]) %>% melt %>% 
  ggplot(., aes(x = variable, y = value, group = interaction(variable, hyp, sp))) +
  geom_boxplot(aes(group = interaction (variable, hyp, sp))) +
  facet_wrap(~hyp, ncol = 3)

Modulemurs$test.modularity.dist %>% ldply(function(x) x$plotaisso[, -3]) %>% melt %>% head
  ggplot(., aes(x = variable, y = value, group = interaction(variable, hyp, sp))) +
  geom_boxplot(aes(group = interaction (variable, hyp, sp))) +
  facet_wrap(~hyp, ncol = 3)

  temp <- 
    Modulemurs$test.modularity.dist %>% llply(function(x) apply (x, 1 , FUN = x$plotaisso$`AVG+`/x$plotaisso$`AVG-`) ) %>% head



Modulemurs$test.modularity.dist$Tarsius_bancanus$plotaisso$`AVG+`/ Modulemurs$test.modularity.dist$Tarsius_bancanus$plotaisso$`AVG-`


