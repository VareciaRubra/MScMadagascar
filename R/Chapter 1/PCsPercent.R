########################################################################################################################################
#############OLHANDO PARA OS PCs DAS MATRIZES: INTERPRETANDO OS PCs E VENDO A DISTRIBUIÇÃO DE VARIAÇÃO##################################
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

registerDoParallel(cores = 4)

MonteCarloPCPercent <- function (otu, iterations = 1000, parallel = FALSE) 
{
  with (otu, 
        {
          return (aaply (1:iterations, 1, 
                 function (i)
                   {
                      covmat <- var (rmvnorm (sample.size, sigma = matrix$cov, method = 'chol'))
                      eval <- eigen (covmat) $ values
                      eval/sum (eval)
                   }, .parallel = parallel))
                 })
}
 
tmp <- MonteCarloPCPercent(sp.main.data$Tarsius_bancanus, iterations = 100, parallel = TRUE)

l

tab.rep.mx.mix



evolqg:::MonteCarloStat







