cov.mx<- sp.master.main.data %>% llply(function(x) x$matrix$cov) %>% {. %*% t(.)}
mx.rep<- sp.master.main.data %>% ldply(function(x) x$Mx.Rep$BootsRep[1]) %>% {. %*% t(.)}
n.size<- sp.master.main.data %>% ldply(function(x) x$sample.size) %>% {. %*% t(.)}


mask1 <- ldply(sp.master.main.data, function(x) x$sample.size >40) [,2]
stats <- llply(sp.master.main.data, function(x) if (!is.na(x$matrix$cov[1])  ) x$matrix$cov else NA)
mask2 <- ldply(stats, function(x) !is.na(x[1]))[,2]
stats <- stats[mask1]

R2<- stats %>% ldply(CalcR2, .progress = progress_text(char = "."), .inform = T)
R2 %>% ggplot() +geom_point( aes(x= .id, y = V1))
plot(R2)