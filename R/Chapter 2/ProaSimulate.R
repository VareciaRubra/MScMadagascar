#Drift test dos dois jeitos, com 3 conjudos de dados  #vai porra!

cov.mx<- sp.main.data %>% llply(function(x) x$matrix$cov) # todas as matrizes de cov
mask.extant <- ldply(sp.main.data, function(x) unique(x$info$Status) == "Extant")[,2] #mascara para os vventes
mask.is.na.cov <- ldply(cov.mx, function(x) is.na(x[1]))[,2] # mascara dos sem matriz
mask.no.na.cov <- ldply(cov.mx, function(x) !is.na(x[1]))[,2] #mascara dos com matriz, mesmo que seja tosca
mask.madagascar <- ldply(sp.main.data, function(x) unique(x$info$Regiao) == "Madagascar")[,2] 
n.size <- sp.main.data %>% ldply(function(x) x$sample.size) %>% .[,2] #tamanho amostral de todos
mask.n.size <- n.size> 30 # mascara de tamanho amostral >30
ed.means <- sp.main.data %>% llply(function(x) x$ed.means) #médias de todos
mask.madagascar & mask.no.na.cov

B.var <- sp.main.data[mask.madagascar & mask.no.na.cov] %>% ldply(function(x) x$ed.means) %>% .[, -1] %>% var

B.all <- sp.main.data[mask.madagascar & mask.no.na.cov] %>% ldply(function(x) x$ed) 

B.all.lm <- lm(as.matrix(B.all[,-1]) ~ B.all[, 1])
B.all.lm.1 <- lm(as.matrix(B.all[,-1]) ~ 1)

sumsqr.W <- t(B.all.lm$residuals) %*% B.all.lm$residuals 
W.lm <- CalculateMatrix(B.all.lm)
W.lm - Ancestral.Matrices$`45`
sumsqr.W.t <- t(B.all.lm.1$residuals) %*% B.all.lm.1$residuals 

sp.main.data[mask.madagascar & mask.no.na.cov] %>% ldply(function(x) x$ed) %>% dim

B.sumsqr <- (sumsqr.W.t - sumsqr.W) / (1318 -1)
MatrixCompare(B.sumsqr, B.var)
all.equal(B.sumsqr, B.var)
MatrixCompare(B.sumsqr, Ancestral.Matrices$`42`)
registerDoParallel(cores = 2)
Proa <- vector("list", 5)

# Quando rodei o teste de deriva mandei ele salvar todos os pares de B e W por nó.
# isso está armazenado em:
Drift.results$extant.sp$BW.compare$

Proa$Sim1.lower <- simulateGP(method = "sim1", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe =1.25, 
                        pop = length(n.size[mask.extant, 2]),
                        n = 40,
                        sim.n = 1000,
                        G = Drift.results$extant.sp$BW.compare$`71`$W, 
                        P = Drift.results$extant.sp$BW.compare$`71`$B.ed)
Proa$Sim1.lower$TypeIerror
Proa$Sim1.lower$Beta.ic

Proa$Sim1.upper <- simulateGP(method = "sim1", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe =700, 
                        pop = length(n.size[mask.extant, 2]),
                        n = 40,
                        sim.n = 1000,
                        G = Drift.results$extant.sp$BW.compare$`71`$W, 
                        P = Drift.results$extant.sp$BW.compare$`71`$B.ed)
Proa$Sim1.upper$TypeIerror
Proa$Sim1.upper$Beta.ic


Proa$Sim2.lower <- simulateGP(method = "sim2", 
                              m = dim(cov.mx[[1]])[1], 
                              tNe =1.25, 
                              pop = length(n.size[mask.extant, 2]),
                              n = 40,
                              sim.n = 1000,
                              G = Drift.results$extant.sp$BW.compare$`71`$W, 
                              P = Drift.results$extant.sp$BW.compare$`71`$B.ed)
Proa$Sim2.lower$TypeIerror
Proa$Sim2.lower$Beta.ic

Proa$Sim2.upper <- simulateGP(method = "sim2", 
                              m = dim(cov.mx[[1]])[1], 
                              tNe =700, 
                              pop = length(n.size[mask.extant, 2]),
                              n = 40,
                              sim.n = 1000,
                              G = Drift.results$extant.sp$BW.compare$`71`$W, 
                              P = Drift.results$extant.sp$BW.compare$`71`$B.ed)
Proa$Sim2.upper$TypeIerror
Proa$Sim2.upper$Beta.ic

save.image("~/MScMadagascar/.RData")

