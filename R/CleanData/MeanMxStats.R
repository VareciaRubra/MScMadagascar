sp.master.main.data <- Strepsirrhini.image$specie.lists

cov.mx <- sp.master.main.data %>% llply(function(x) x$matrix$cov)
mx.rep <- sp.master.main.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) 
n.size <- sp.master.main.data %>% ldply(function(x) x$sample.size) 
ed.means <- sp.master.main.data %>% llply(function(x) x$ed.means) 
gm.mean <- sp.master.main.data %>% ldply(function(x) x$gm.mean) 

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


###Sample sizes #####
sample.no.na <- n.size[mask.na.cov,2]
names(sample.no.na) <- n.size[,1][mask.na.cov]
##### Cov mx ########
cov.no.na <- cov.mx[mask.na.cov]
##### Mx repetabilities via RS ######
rep.no.na <- mx.rep[mask.na.cov,2]
names(rep.no.na) <- mx.rep[,1][mask.na.cov]
##### ED means ######
ed.means <- ed.means[mask.na.cov]

names(ed.means) == names(rep.no.na) == names(cov.mx[mask.na.cov]) == mx.rep[,1][mask.na.cov]





