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

Proa$Sim1 <- 
  temp <- simulateGP(method = "sim1", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe = 0.001, 
                        pop = length(n.size[mask.madagascar & mask.no.na.cov]),
                        n = 40,
                        sim.n = 1000,
                        G = B.sumsqr, 
                        P = Ancestral.Matrices$`45`)

Proa$Sim2 <- simulateGP(method = "sim2", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe = 0.001, 
                        pop = length(n.size[mask.madagascar & mask.no.na.cov]),
                        n = 40,
                        sim.n = 1000,
                        G = B.sumsqr, 
                        P = Ancestral.Matrices$`45`)

Proa$Sim3 <- simulateGP(method = "sim3", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe = 0.001, 
                        pop = length(n.size[mask.madagascar & mask.no.na.cov]),
                        n = 40,
                        sim.n = 1000,
                        G = B.sumsqr, 
                        P = Ancestral.Matrices$`45`)
Proa$Sim4 <- simulateGP(method = "sim4", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe = 0.001, 
                        pop = length(n.size[mask.madagascar & mask.no.na.cov]),
                        n = 40,
                        sim.n = 1000,
                        G = B.sumsqr, 
                        P = Ancestral.Matrices$`45`)
Proa$Sim5 <- simulateGP(method = "sim5", 
                        m = dim(cov.mx[[1]])[1], 
                        tNe = 0.001, 
                        pop = length(n.size[mask.madagascar & mask.no.na.cov]),
                        n = 40,
                        sim.n = 1000,
                        G = B.sumsqr, 
                        P = Ancestral.Matrices$`45`)

save.image("~/MScMadagascar/.RData")

