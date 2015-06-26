current.data <- sp.main.data

cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) 
n.size <- current.data %>% ldply(function(x) x$sample.size) 
ed.means <- current.data %>% llply(function(x) x$ed.means) 
gm.mean <- current.data %>% ldply(function(x) x$gm.mean) 

mask.n.size <- n.size[,2]>40
mask.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
mask.rep <- !is.na(mx.rep)[,2]

distance <- llply(cov.mx[mask.na.cov], MatrixDistance, distance = "OverlapDist")
RS<- llply(cov.x= cov.mx[mask.n.size], RandomSkewers, num.vectors = 1000, repeat.vector = mx.rep[mask.n.size] )
plot.matrix.cor(RS$correlations, main = "V/CV Matrix compared by Random Skewers", brewer = "RdBu")
mask2 <- ldply(stats, function(x) !is.na(x[1]))[,2]
stats <- stats[mask1]

R2<- stats %>% ldply(CalcR2, .progress = progress_text(char = "."), .inform = T)
R2 %>% ggplot() +geom_point( aes(x= .id, y = V1))
plot(R2)

RS.compare <- RandomSkewers(cov.x = cov.40, repeat.vector = rep.40, num.vectors = 1000, parallel = TRUE)
coord