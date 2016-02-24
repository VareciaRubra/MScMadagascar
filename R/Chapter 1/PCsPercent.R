########################################################################################################################################
############ OLHANDO PARA OS PCs DAS MATRIZES: INTERPRETANDO OS PCs E VENDO A DISTRIBUIÇÃO DE VARIAÇÃO #################################
########################################################################################################################################
current.data <- sp.main.data
current.data <- gen.main.data

cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
cor.mx <- current.data %>% llply(function(x) x$matrix$cor)
mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
mx.rep.mc <- sp.main.data.t%>% ldply(function(x) x$Mx.Rep$MCRep) 
mx.rep.all <- cbind(mx.rep, mx.rep.mc[,-1])
n.size <- current.data %>% ldply(function(x) x$sample.size)

mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
names (cov.mx[mask])

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

######################################################
###### CHUPINHADORA DE MONTE CARLO ###################
######################################################
## OLHANDO A CARA DA MonteCarloR2

require (mvtnorm)

registerDoParallel(cores = 2)

MonteCarloPCPercent <- function (otu, iterations = 1000, parallel = FALSE) 
{## B_A = t(E) B E
  ## B é a sua matriz a ser rotacionada
  ## A é a sua matriz original referencia
  ## E sao os autovalores de A
  with (otu, 
        {
          Aevec <- eigen (matrix$cov) $ vectors
          Aeval <- eigen (matrix$cov) $ values
          PCpercent <- as.matrix(Aeval/sum (Aeval) )
          
          projec.B.A <- aaply (1:iterations, 1, .fun = 
                 function (i)
                   {  B <- var (rmvnorm (sample.size, sigma = matrix$cov, method = 'svd') ) ## gerando valores a partir de uma distribuiçao que tem a cara da sua matriz
                      B.rot <- t(Aevec) %*% B %*% Aevec ## calculando a projeçao de B na A
                      diag(B.rot) / sum(Aeval) ## normalizando pelos autovalores da matriz original
                      
                   }, .parallel = parallel) 
          
          gms <- aaply (1:iterations, 1, 
                 function (i)
                 {
                   eds <- rmvnorm (sample.size, mean= ed.means, sigma = matrix$cov, method = 'svd') ## gerando valores a partir duma distribuiçao com o jeitao da sua matriz e com média igual a da sp
                   mean(aaply (eds, 1, function (x) exp (mean (log (x) ))) )
                   # calculando a média geométrica média naquela populaçao amostrada nessa iteraçao
                   }, .parallel = parallel)

          intervalo.mc.pc <- projec.B.A %>% adply(., 2, function (x) sort(x) [ c(3,97) ] )
          names(intervalo.mc.pc) <- c("PC", "min", "max")
          intervalo.mc.pc$observed <- PCpercent   
        
        
         intervalo.mc.gm <- data.frame("min" = sort(gms)[3], "max" = sort(gms)[97])
         intervalo.mc.gm$observed <- gm.mean  
         temp$intervalo.mc.gm <- as.data.frame(temp$intervalo.mc.gm)
         
Plot.PC <- temp$intervalo.mc.pc %>%
           ggplot(.) +
           geom_errorbar( aes(x= PC, ymin = min, ymax = max )) +
           geom_point(aes(x= PC, y = observed), color = "red") +
           theme_bw()
         
Plot.gm <- temp$intervalo.mc.gm %>%
           ggplot(.) +
           geom_errorbar( aes(x= 1, ymin = min, ymax = max )) +
           geom_point(aes(x= 1, y = observed), color = "red") +
           theme_bw()

Plotao <- plot_grid( Plot.PC, Plot.gm)
          return(list("intervalo.mc.pc" = intervalo.mc.pc,
                      "intervalo.mc.gm" = intervalo.mc.gm,
                      "gms" = gms,
                      "Plot" = Plotao))
         
                 })
}
 
temp <- MonteCarloPCPercent (sp.main.data$Microcebus_griseorufus, iterations = 100, parallel = TRUE)


Vairiae <- llply(sp.main.data[mask], .fun =  MonteCarloPCPercent, .progress = "text")


MonteCarloPCPercent.nowith <- function (main.data, iterations = 1000, parallel = FALSE) 
{## B_A = t(E) B E
  ## B é a sua matriz a ser rotacionada
  ## A é a sua matriz original referencia
  ## E sao os autovalores de A
  main.data <- sp.main.data$Microcebus_griseorufus
  
  Aevec <- eigen (main.data$matrix$cov) $ vectors
  Aeval <- eigen (main.data$matrix$cov) $ values
  PCpercent <- as.matrix(Aeval/sum (Aeval) )
  
  projec.B.A <- aaply (1:iterations, 1, .fun = 
                         function (i)
                         {  B <- var (rmvnorm (main.data$sample.size, sigma = main.data$matrix$cov, method = 'svd') ) ## gerando valores a partir de uma distribuiçao que tem a cara da sua matriz
                         B.rot <- t(Aevec) %*% B %*% Aevec ## calculando a projeçao de B na A
                         diag(B.rot) / sum(Aeval) ## normalizando pelos autovalores da matriz original
                         
                         }, .parallel = parallel) 
  
  gms <- aaply (1:iterations, 1, 
                function (i)
                {
                  eds <- rmvnorm (main.data$sample.size, mean= main.data$ed.means, sigma = main.data$matrix$cov, method = 'svd') ## gerando valores a partir duma distribuiçao com o jeitao da sua matriz e com média igual a da sp
                  mean(aaply (eds, 1, function (x) exp (mean (log (x) ))) )
                  # calculando a média geométrica média naquela populaçao amostrada nessa iteraçao
                }, .parallel = parallel)
  
  intervalo.mc.pc <- projec.B.A %>% adply(., 2, function (x) sort(x) [ c(3,97) ] )
  names(intervalo.mc.pc) <- c("PC", "min", "max")
  intervalo.mc.pc$observed <- PCpercent   
  
  
  intervalo.mc.gm <- data.frame("min" = sort(gms)[3], "max" = sort(gms)[97])
  intervalo.mc.gm$observed <- gm.mean  
  temp$intervalo.mc.gm <- as.data.frame(temp$intervalo.mc.gm)
  
  Plot.PC <- temp$intervalo.mc.pc %>%
    ggplot(.) +
    geom_errorbar( aes(x= PC, ymin = min, ymax = max )) +
    geom_point(aes(x= PC, y = observed), color = "red") +
    theme_bw()
  
  Plot.gm <- temp$intervalo.mc.gm %>%
    ggplot(.) +
    geom_errorbar( aes(x= 1, ymin = min, ymax = max )) +
    geom_point(aes(x= 1, y = observed), color = "red") +
    theme_bw()
  
  Plotao <- plot_grid( Plot.PC, Plot.gm)
  return(list("intervalo.mc.pc" = intervalo.mc.pc,
              "intervalo.mc.gm" = intervalo.mc.gm,
              "gms" = gms,
              "Plot" = Plotao))
  
  
}

temp <- MonteCarloPCPercent (sp.main.data$Indri_indri, iterations = 100, parallel = TRUE)
  

















MonteCarloVarSize <- function (otu, iterations = 1000, parallel = FALSE) 
{
  with (otu, 
        {
          return (aaply (1:iterations, 1, 
                         function (i)
                         {
                           eds <- (rmvnorm (sample.size, mean= ed.means, sigma = matrix$cov, method = 'chol'))
                           gms <- aaply (eds, 1, function (x) exp (mean (log (x) )))
                           var (gms)
                         }, .parallel = parallel))
          
          
        })
}
tmp <- MonteCarloVarSize (sp.main.data$Indri_indri, iterations = 100, parallel = TRUE)

gm.var.dist <- llply(sp.main.data[mask], .fun = MonteCarloVarSize, .progress = "text")

tab.rep.mx.mix

exp (mean (log (sp.main.data$Indri_indri$ed.means))) / sqrt



tmp * 39





