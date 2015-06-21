cov.mx <- sp.master.main.data %>% llply(function(x) x$matrix$cov)
mx.rep <- sp.master.main.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) 
n.size <- sp.master.main.data %>% ldply(function(x) x$sample.size) 

mask.n.size <- n.size[,2]>40
mask.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]
mask.rep <- !is.na(mx.rep)[,2]

distance <- llply(cov.mx[mask.na.cov], MatrixDistance, distance = "OverlapDist")
RS<- llply(cov.mx[mask.n.size], RandomSkewers, num.vectors = 1000, repeat.vector = )
mask2 <- ldply(stats, function(x) !is.na(x[1]))[,2]
stats <- stats[mask1]

R2<- stats %>% ldply(CalcR2, .progress = progress_text(char = "."), .inform = T)
R2 %>% ggplot() +geom_point( aes(x= .id, y = V1))
plot(R2)