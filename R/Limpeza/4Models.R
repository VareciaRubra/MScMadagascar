
#############################################
############# PRIMEIRA PARTE ################
##### testando modelos de efeitos fixos #####
################ SEXO ######################
#############################################
current.data <- sp.main.data$
SexSig <- function (current.data) 
{
  y = vector("list", 12)
  y [[1]] <-  if (sum(table(current.data$info$Sexo))>= 43) manova(as.matrix(current.data$ed) ~ Sexo, data = as.data.frame(current.data$info) ) else NA
  y [[2]] <- if (sum(table(current.data$info$Sexo))>= 43) Manova(y[[1]], type=3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly")) else NA
  y [[3]]  <- if (table(current.data$info$Sexo)[1] >= 3 & table(current.data$info$Sexo)[2] >= 3) 
              apply(as.matrix(current.data$ed ), 2, function (x) return (Anova(lm(x ~ Sexo, data = as.data.frame(current.data$info) ), type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))) ) else NA
  y [[4]] <- if (sum(table(current.data$info$Sexo))>=43)  CalculateMatrix(y[[1]]) else current.data$matrix$cov
  y [[5]] <- current.data$matrix$cov
  names(y)[1:5] <- c("fit", "multi", "uni",  "cov.fit", "cov")
  y[[6]] <- if (!is.na(y$cov.fit[1]) ) RandomSkewers (cov.x = y$cov.fit, cov.y = y$cov , num.vectors = 10000, parallel = TRUE)$correlation else NA
  y[[7]] <- if (!is.na(y$cov.fit[1]) ) KrzProjection(cov.x = y$cov.fit,  cov.y = y$cov)$total.variation else NA
  y[[8]] <- if (!is.na(y$cov.fit[1]) ) PCAsimilarity(cov.x = y$cov.fit, cov.y = y$cov ) else NA
  y[[9]] <- if (!is.na(y$cov.fit[1]) ) KrzProjection(cov.x = cov2cor(y$cov.fit), cov.y = cov2cor(y$cov) ) $total.variation else NA
  y[[10]] <- if (!is.na(y$cov.fit[1]) ) PCAsimilarity(cov.x = cov2cor(y$cov.fit), cov.y = cov2cor(y$cov) ) else NA
  names(y)[6:10] <- c("RS", "KRZ", "PCA.s",  "cor.KRZ", "cor.PCA.s")
  sig <- as.data.frame(cbind(names(current.data$ed), names(current.data$ed)))
  sig[,2] <- NA
  names(sig) <- c(".ed", "V1")
  y[[11]] <- if(!is.na( y [[4]][1]) & !is.na(y$uni)) ldply(y$uni, .fun = function(x) x[4][2,] ) else sig
  names(y)[11] <- "sig.sex.uni"
  y[[12]] <- 
  names(y)[12] <- "sig.sex.multi"
  print (unique(current.data$info$Especie))
  return(y)
  
}

sex.sig <- llply(sp.main.data, SexSig, .progress = 'text')

sex.sig %>% ldply(function(x) x$RS[1])
sex.sig %<>% ldply(function(x) x$sig.sex[,2]) 
xtable(na.omit(sex.sig), digits = 2)
multi.sig <- sex.sig %>% llply(function(x) summary(x$multi) ) 
multi.sig <- sex.sig %>% llply(function(x) x$multi )
multi.sig %>% !is.na() %>% ldply(!is.na(multi.sig), getTable)
#############################Tabela pra meter no relatório ################################
rownames(sex.sig) <- sex.sig[,1]
colnames(sex.sig) <- c("Especie", dimnames(sp.main.data$Tarsius_bancanus$matrix$cov)[[1]]) 
table.sex.sig <- xtable(na.omit(sex.sig)[,-1], digits = 3)


getTable <- function (x, ...) 
{
  if ((!is.null(x$singular)) && x$singular) 
    stop("singular error SSP matrix; multivariate tests unavailable\ntry summary(object, multivariate=FALSE)")
  test <- x$test
  repeated <- x$repeated
  ntests <- length(x$terms)
  tests <- matrix(NA, ntests, 4)
  if (!repeated) 
    SSPE.qr <- qr(x$SSPE)
  for (term in 1:ntests) {
    eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else SSPE.qr, 
                             x$SSP[[term]]), symmetric = FALSE)$values)
    tests[term, 1:4] <- switch(test, Wilks = car:::Wilks(eigs, x$df[term],  x$error.df), x$error.df)
  }
  ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
  ok <- !is.na(ok) & ok
  tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], 
                                 tests[ok, 4], lower.tail = FALSE))
  rownames(tests) <- x$terms
  colnames(tests) <- c("Df", "test stat", "approx F", "num Df", 
                       "den Df", "Pr(>F)")
  tests <- structure(as.data.frame(tests), heading = paste("\nType ", 
                                                           x$type, if (repeated) 
                                                             " Repeated Measures", " MANOVA Tests: ", test, " test statistic", 
                                                           sep = ""), class = c("anova", "data.frame"))
  tests
}

################# Student t test pra sexo e GM das variaveis cranianas

t.test.sex.gm <- sp.main.data %>% ldply(function (x) table(x$info$Sexo, useNA = "always") ) %>% llply()
sp.main.data %>% llply(function(x) t.test(x$gm.ind[!is.na(x$info$Sexo)] ~ x$info$Sexo[!is.na(x$info$Sexo)]) ) 

t.test(sp.main.data$Tarsius_bancanus$ ~ )

