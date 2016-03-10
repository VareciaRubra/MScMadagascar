#Montando as tabelas de Repetibilidades


# faltou calcular as repetibilidades por monte carlo
AddMCReps <- function(x){
  x[[16]][[2]] <- matrix(NA, nrow = 1, ncol= 2)
  x[[16]][[2]][1] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = RandomSkewers, iterations = 1000, parallel = T) else NA
  x[[16]][[2]][2] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = KrzCor, iterations = 1000, parallel = T) else NA
  #x[[16]][[2]][3] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = PCAsimilarity, iterations = 1000, parallel = T) 
  #x[[16]][[2]][4] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = MantelCor, iterations = 1000, correlation = TRUE, parallel = T) 
  #x[[16]][[2]][5] <- if(x[[9]]>14) MonteCarloRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = KrzCor, iterations = 1000, correlation = TRUE, parallel = T)
  colnames(x[[16]][[2]])[1:2] <- c('rs', 'krz')
  names(x[[16]])[2] <- c('MCRep')
  return(x)
  
}
registerDoParallel(cores = 2)
current.data <- sp.main.data
sp.main.data.t <- current.data %>% llply(AddMCReps, .progress = "text")

mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
names(mx.rep) <- c("Especie", "BS.rs", "BS.krz", "pcas", "cor.mantel", "cor.krz")
mx.rep.mc <- sp.main.data.t%>% ldply(function(x) x$Mx.Rep$MCRep) 
names(mx.rep.mc) <- c("Especie", "MC.rs", "MC.krz")
tab.rep.mx <- cbind(mx.rep[,1:3], mx.rep.mc[,2:3])
tab.rep.mx$Especie %<>% gsub("_", ' ', .)
row.names(tab.rep.mx) <- tab.rep.mx[,1]
tab.rep.mx$n <- n.size$V1

cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]

rep.list <- mx.rep[mask, ]
row.names(rep.list) <- rep.list$Especie
rep.list$Especie <- as.character(rep.list$Especie)
rep.list[dim(rep.list)[1] +1, ] <- c("Saguinus_P", rep(0.97, 5) )
rep.list[dim(rep.list)[1] +1, ] <- c("Saguinus_G", rep(0.75, 5) )
row.names(rep.list) <- rep.list$Especie
rep.list$Especie <- factor(rep.list$Especie, levels = unique(rep.list$Especie) )
rep.list$rs <- as.numeric(rep.list$rs)
rep.list$krz <- as.numeric(rep.list$krz)
rep.list$pcas <- as.numeric(rep.list$pcas)
rep.list$cor.mantel <- as.numeric(rep.list$cor.mantel)
rep.list$cor.krz <- as.numeric(rep.list$cor.krz)
str(rep.list)

####################### calculando o valor médio de comparação pra cada especie ##########

mean.sim<-function(x) {
  x[upper.tri(x)]<-t(x)[upper.tri(x)]
  diag(x)<-NA
  x<-rowMeans(x,na.rm=TRUE)
  return(x)
}
mean.compare.cov<- mx.compare[c(1:2)] %>% llply (function (x) x$correlations ) %>% ldply (function (x) mean.sim(x) )

dimnames(mean.compare.cov)[[1]] <- c(paste("Mean", mean.compare.cov[,1] ))
dimnames(mean.compare.cov)[[2]] <- dimnames(mean.compare.cov)[[2]]

row.names(mean.compare.cov)

mean.compare.cov<- t(mean.compare.cov[,-1])
tab.rep.mx.mix<- cbind(tab.rep.mx[mask,-1] , mean.compare.cov[-c(43:44),] )

xtable(tab.rep.mx.mix, digits = 3, caption = "Matrix repetabilities calculated via Bootstrap and Monte Carlo using RS and KRZ methods of comparisson.")


MonteCarloCompleteRep <- function (cov.matrix, sample.size, ComparisonFunc, ..., iterations = 1000, 
          correlation = FALSE, parallel = TRUE) 
{
  if (correlation) {
    StatFunc <- cov
    c2v <- cov2cor
  }
  else {
    StatFunc <- cov
    c2v <- function(x) x
  }
  repeatability <- MonteCarloStat(cov.matrix, sample.size, 
                                  iterations, ComparisonFunc = function(x, y) ComparisonFunc(c2v(x), 
                                                                                             c2v(y), ...), StatFunc = StatFunc, parallel = parallel)
  return(repeatability)
}

AddMCRepDists <- function(x){
  x[[16]][[3]] <- vector("list", 2)
  x[[16]][[3]][[1]]<- if(x[[9]]>14) MonteCarloCompleteRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = RandomSkewers, iterations = 100, parallel = T)$correlation else NA
  x[[16]][[3]][[2]] <- if(x[[9]]>14) MonteCarloCompleteRep( x[[11]][[1]], sample.size = x[[9]], ComparisonFunc = KrzCor, iterations = 100, parallel = T)$V1 else NA
  names(x[[16]][[3]])[1:2] <- c('rs', 'krz')
  names(x[[16]])[3] <- c('MCRepDist')
  return(x)
  
}

registerDoParallel(cores = 2)
sp.main.data.t <- sp.main.data.t  %>% llply(AddMCRepDists, .progress = "text")

mx.rep.mc.dist <- sp.main.data.t[mask] %>% ldply(function(x) data.frame("rs"= x$Mx.Rep$MCRepDist$rs, "krz" = x$Mx.Rep$MCRepDist$krz ) )

# mx.rep.mc.dist <- sp.main.data.t[mask] %>% ldply(function(x) data.frame("min.rs"= min(x$Mx.Rep$MCRepDist$rs), 
#                                                                         "mean.rs" = mean(x$Mx.Rep$MCRepDist$rs),
#                                                                         "max.rs"= max(x$Mx.Rep$MCRepDist$rs), 
#                                                                         "min.krz" = min(x$Mx.Rep$MCRepDist$krz),
#                                                                         "mean.krz" = mean(x$Mx.Rep$MCRepDist$krz),
#                                                                         "max.krz" = max(x$Mx.Rep$MCRepDist$krz) ) )

#1 pegar os ranges das repetibilidades 

#mx.rep.mc.dist$Especie <- as.factor(mx.rep.mc.dist$.id)
mx.rep.mc.dist$.id %<>% gsub("_", " ",.) 
mx.rep.mc.dist$.id <- factor(mx.rep.mc.dist$.id, levels = unique(mx.rep.mc.dist$.id)[42:1])
str(mx.rep.mc.dist)



mx.rep.mc.dist.sumarry <- mx.rep.mc.dist %>% melt(variable.name = "rep") %>% group_by(rep, .id) %>%
  summarise(mean = mean(value), min = min(value), max= max(value) )  
RepPlotMC <- 
  mx.rep.mc.dist.sumarry  %>% 
  ggplot(.) + 
  geom_linerange(aes(x = .id, ymin = min, ymax = max, color = rep), alpha = 0.1, size = 2) +
  geom_point(aes(x = .id, y = mean, color = rep )) +
  #facet_wrap(~ PC, scales = 'free') +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  xlab("") + ylab("") + labs(title = "Matrices repetabilities via MC") +
  coord_flip() +
  theme_bw() +
  facet_grid(~rep) + 
  theme(legend.position="none",
        axis.text.y = element_text(face =  "italic", size =10),
        axis.text.x = element_text(size =19),
        strip.text= element_text(size=19),
        plot.title = element_text(face = "bold", size = 20)) 

 mx.rep.mc.dist %>% melt(variable.name = "rep") %>% group_by(rep, .id)  %>%
  ggplot(.) +
   geom_violin(aes(x = .id, y =  value, color = rep), alpha = 0.4) +
   geom_jitter(aes(x = .id, y =  value, color = rep), alpha = 0.3, size = 1, shape = 4) +
  coord_flip() +
  theme_bw() +
  facet_grid(~rep) + 
  scale_colour_brewer(palette = "Set1", direction = -1) +
   xlab("")+
  ylab( "Matrix repetabilities (MC)")+
  theme(legend.position="none",
        axis.text.y = element_text(face =  "italic", size =10),
        axis.text.x = element_text(size =19),
        strip.text= element_text(size=19),
        plot.title = element_text(face = "bold", size = 20)) 


SummaySim <- function (x, tri = lower.tri) 
  { ## x = similarity matrix
    ## tri = qual triangulo da matriz usar (valores corrigidos estao na upper.tr)
      x[tri(x)]<-t(x)[tri(x)]
      diag(x)<-NA
      diag(x)<-rowMeans(x, na.rm=TRUE)
      x<-summary(x, na.rm=TRUE)
return(x) } 

SimSim <- function (x, tri = upper.tri) 
{ ## x = similarity matrix
  ## tri = qual triangulo da matriz excluir (valores corrigidos estao na upper.tr)
  x[tri(x)]<-t(x)[tri(x)]
  diag(x)<-NA
  diag(x)<-rowMeans(x, na.rm=TRUE)
  return(x) } 
 
mx.violins <- mx.compare[6:7] %>% llply (function (x) t(x$correlations) )
mx.violins$Steppan.RS <- Steppan$MC.RS
mx.violins$Steppan.KRZ <- Steppan$MC.KRZ

mx.comp.values <- mx.compare[6:7] %>% llply (function (x) x$correlations[-c(43:44), -c(43:44)] ) %>%  ldply (function (x) SimSim(x) ) %>%  melt(variable.name = ".sp") %>% group_by(.id, .sp) 

mx.comp.values$.sp %<>% gsub("_", " ", .)
mx.comp.values$.sp <- factor(mx.comp.values$.sp, levels = unique(mx.comp.values$.sp)[42:1])
mx.comp.values$.id %<>% gsub("MC.", "", .) %>% gsub("KRZ", "raw.krz", .)
mx.comp.values$.id %<>% gsub("MC.", "", .) %<>% gsub("RS", "raw.rs", .)
mx.comp.values$.id<- factor(mx.comp.values$.id, levels = unique(mx.comp.values$.id)[1:2])
names( mx.comp.values) <- c("rep", ".sp", "sim")
head(mx.comp.values) 

mx.violins<- mx.violins %>% llply (function (x) x[-c(43:44), -c(43:44)] ) %>% ldply (function (x) SimSim(x) ) %>%  melt(variable.name = ".sp") %>% group_by(.id, .sp) 

mx.violins$.sp %<>% gsub("_", " ", .)
mx.violins$.id %<>% gsub("Steppan.KRZ", "sample.size.krz", .) 
mx.violins$.id %<>% gsub("Steppan.RS", "sample.size.rs", .) 
mx.violins$.id %<>% gsub("MC.KRZ", "mc.rep.krz", .) 
mx.violins$.id %<>% gsub("MC.RS", "mc.rep.rs", .) 
mx.violins$.id<- factor(mx.violins$.id, levels = unique(mx.violins$.id)[c(2,1,4,3)])
mx.violins$.sp <- factor(mx.violins$.sp, levels = unique(mx.violins$.sp)[42:1])
names( mx.violins) <- c("rep", ".sp", "sim")
mx.violins %>% head
mx.comp.values <- rbind(mx.comp.values, mx.violins)
mx.comp.values$rep <- factor(mx.comp.values$rep, levels = unique(mx.comp.values$rep)[1:6])

mx.comp.values.summary <- mx.comp.values %>% group_by(rep, .sp) %>% summarise(mean = mean(sim), min = min(sim), max= max(sim) ) 
names(mx.comp.values.summary) <- c( "rep", ".sp", "sim.mean", "sim.min", "sim.max")
mx.comp.values.summary$rep <- factor(mx.comp.values.summary$rep, levels = unique(mx.comp.values.summary$rep)[1:6])


mx.rep.mc.dist <- sp.main.data.t[mask] %>% ldply(function(x) data.frame("rs"= x$Mx.Rep$MCRepDist$rs, "krz" = x$Mx.Rep$MCRepDist$krz ) )
mx.rep.mc.dist <-  mx.rep.mc.dist %>% melt(variable.name = "rep") %>% group_by(rep, .id) %>%  summarise(mean = mean(value), min = min(value), max= max(value) )
mx.rep.mc.dist$rep %<>% gsub("rs", "mc.rep.rs", . ) %>% gsub ("krz", "mc.rep.krz", .)
mx.rep.mc.dist$rep <- factor(mx.rep.mc.dist$rep, levels = unique(mx.comp.values.summary$rep)[1:6])
mx.rep.mc.dist$.id %<>% gsub("_", " ", .)
mx.rep.mc.dist$.id <- factor(mx.rep.mc.dist$.id, levels =  levels(mx.violins$.sp))
str(mx.rep.mc.dist)

medianas<- mx.comp.values %>% group_by(rep, .sp) %>% summarise("median" = median(sim), "quantile.4" = quantile(sim)[4])
names(medianas)
medianas$sample.size <- n.size[mask,2][42:1]
medianas$sample.size.1 <- c(rep(NA, 210), n.size[mask,2][42:1])
medianas$sample.size.2 <- c(n.size[mask,2][42:1],n.size[mask,2][42:1], rep(NA, length(medianas$sample.size) - 84)) 
m.lines<- medianas %>% group_by(rep) %>% summarise(median(median)) %>% as.data.frame
q.lines <- medianas %>% group_by(rep) %>% summarise(median(quantile.4)) %>% as.data.frame
medianas$line  <- rep(m.lines[1:6,2], each=42) 
medianas$q.line  <- rep(q.lines[1:6,2], each=42) 
medianas$rep <- factor(medianas$rep, levels = unique(medianas$rep)[1:6])

PlotSumarySym<- 
mx.comp.values %>% #head
  ggplot(.) + 
  #geom_point(data = medianas, aes(x = .sp, y =  median, size = sample.size.2, color = rep), shape = 16, alpha = 0.6) +
  geom_text(data = medianas, aes(x = .sp, y =  1, label = sample.size.1),size = 4, color = "darkgrey" ) +
  geom_hline(data = medianas, aes(x = .sp, yintercept = line, color = rep), size = 0.8, linetype = "dashed", alpha = 0.02) +
  geom_hline(data = medianas, aes(x = .sp, yintercept = q.line, color = rep), size = 1, linetype = "solid", alpha = 0.02) +
  geom_violin(aes(x = .sp, y =  sim, color = rep), alpha = 0.4) +
  geom_jitter(aes(x = .sp, y =  sim, color = rep), alpha = 0.3, size = 1, shape = 4) +
  geom_linerange(data = mx.rep.mc.dist, aes(x = .id, ymin = min, ymax = max), color = "darkgrey", alpha = 0.7, size = 2.5) +
  geom_point(data = mx.rep.mc.dist, aes(x = .id, y = mean), color = "white" , shape = 20, size = 1) +
  scale_colour_brewer(palette = "Set1", direction = 1) +
  xlab("") + ylab("") + labs(title = "Summary of matrix comparisons")  +
  #scale_x_continuous(limits = c(0, 0.9), breaks = c(0.2, 0.4, 0.6, 0.8)) 
  coord_flip() +
  theme_gray()  +  facet_wrap(~rep, nrow = 1) + 
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face =  "italic", size =9),
        axis.text.x = element_text(size =9), 
        strip.text = element_text(face = "bold", size = 10),
        legend.position="none") 


 
  
