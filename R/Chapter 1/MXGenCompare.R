current.data <- gen.main.data

gen.cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
str(gen.cov.mx)
gen.mx.rep <- current.data %>% ldply(function(x) x$Mx.Rep$BootsRep) 
gen.n.size <- current.data %>% ldply(function(x) x$sample.size) 

gen.mask <- ldply(gen.cov.mx, function(x) !is.na(x[1]))[,2]

gen.cov.list <- gen.cov.mx[gen.mask]
gen.cov.list$Saguinus_P.cov <- Saguinus_P.cov
gen.cov.list$Saguinus_G.cov <- Saguinus_G.cov
str(gen.cov.list)
#names(gen.cov.list) %<>% gsub("_", ' ', .)


gen.rep.list <- gen.mx.rep[gen.mask,]
row.names(gen.rep.list) <- gen.rep.list$Genero
as.data.frame(gen.rep.list)
gen.rep.list[dim(gen.rep.list)[1] +1, ] <- c("Saguinus_P", rep(0.97, 5) )
gen.rep.list[dim(gen.rep.list)[1] +1, ] <- c("Saguinus_G", rep(0.75, 5) )
row.names(gen.rep.list) <-  c(row.names(gen.rep.list)[-c(22:23)], "Saguinus_P", "Saguinus_G")
gen.rep.list$Genero <- row.names(gen.rep.list)
names(gen.rep.list) %<>% gsub("_", ' ', .)

gen.sample.size.list <- c(gen.n.size[gen.mask,2], 130, 230)


gen.mx.compare = vector("list", 2)
gen.mx.compare[1:2] <- NA
gen.mx.compare[[1]] <- RandomSkewers(cov.x = gen.cov.list, num.vectors = 1000, repeat.vector = as.numeric(gen.rep.list[, 2]) )
gen.mx.compare[[2]]$correlations <- as.matrix(KrzCor(cov.x= gen.cov.list, num.vectors = 1000, repeat.vector = as.numeric(gen.rep.list[, 3]) ))
names(gen.mx.compare)[1:2] <-  c('BS.RS', 'BS.KRZ')
mx.class<- c('V/CV', 'V/CV')
for (i in 1:2)  {gen.mx.compare[[i]]$method <- names(gen.mx.compare)[i]}
for (i in 1:2)  {gen.mx.compare[[i]]$mx.class <- mx.class[i]}

Combine.Mx.Plot(Mx1 = t(gen.mx.compare$BS.RS$correlations), Mx2 = t(gen.mx.compare$BS.KRZ$correlations), prob = gen.mx.compare$BS.RS$probabilities, diag.info = gen.sample.size.list, titulo = "Genus matrices comparisons")


