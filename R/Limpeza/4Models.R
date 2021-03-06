
#############################################
############# PRIMEIRA PARTE ################
##### testando modelos de efeitos fixos #####
################ SEXO ######################
#############################################
current.data <- sp.main.data$Microcebus_griseorufus
SexSig <- function (current.data) 
{
  y = vector("list", 12)
  y[1:12] <- NA
  y [[1]] <-  if (sum(table(current.data$info$Sexo)) >= 43) manova(as.matrix(current.data$ed) ~ Sexo, data = as.data.frame(current.data$info) ) else NA
  y [[2]] <- if (sum(table(current.data$info$Sexo)) >= 43) Manova(y[[1]], type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly")) else NA
  y [[3]]  <- if (table(current.data$info$Sexo)[1] >= 3 & table(current.data$info$Sexo)[2] >= 3) 
              apply(as.matrix(current.data$ed ), 2, function (x) return (Anova(lm(x ~ Sexo, data = as.data.frame(current.data$info) ), type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))) ) else NA
  y [[4]] <- if (sum(table(current.data$info$Sexo)) >= 43)  CalculateMatrix(y[[1]]) else current.data$matrix$cov
  y [[5]] <- current.data$matrix$cov
  names(y)[1:5] <- c("fit", "multi", "uni",  "cov.fit", "cov")
  y[[6]] <- if (!is.na(y$cov.fit[1]) ) RandomSkewers (cov.x = y$cov.fit, cov.y = y$cov , num.vectors = 10000, parallel = TRUE)[1] else NA
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
  y[[12]] <- getTable(y[[2]])[[6]][2] 
  names(y)[12] <- "sig.sex.multi"
  print (unique(current.data$info$Especie))
  return(y)
  
}

sex.sig <- llply(sp.main.data, SexSig, .progress = 'text')

sex.sig %>% ldply(function(x) x$RS[1])
sex.sig %>% llply(function(x) x$sig.sex.uni[,2][!is.na(x$sig.sex.uni[,2])]) %>% 

xtable(na.omit(sex.sig), digits = 2)
multi.sig <- sex.sig %>% llply(function(x) summary(x$multi) ) 
multi.sig <- sex.sig %>% llply(function(x) x$multi )
multi.sig[!is.na(multi.sig)] %>% ldply(function(x) getTable(x)[[6]][2]) %>% xtable(digits = 4)

fixp <- function(x, dig=3){
  x <- getTable(x)
  
  if(substr(names(x)[ncol(x)],1,2) != "Pr")
    warning("The name of the last column didn't start with Pr. This may indicate that p-values weren't in the last row, and thus, that this function is inappropriate.")
  x[,ncol(x)] <- round(x[,ncol(x)], dig)
  for(i in 1:nrow(x)){
    if(x[i,ncol(x)] == 0)
      x[i,ncol(x)] <- paste0("< .", paste0(rep(0,dig-1), collapse=""), "1")
  }
  
  x
}
multi.sig[!is.na(multi.sig)] %>% ldply(fixp) 

#############################Tabela pra meter no relatório ################################
rownames(sex.sig) <- sex.sig[,1]
colnames(sex.sig) <- c("Especie", dimnames(sp.main.data$Tarsius_bancanus$matrix$cov)[[1]]) 
table.sex.sig <- xtable(na.omit(sex.sig)[,-1], digits = 3)

##################3 Comparação entre matrizes de sexo e raw e residual ####################

SexCompare <- function (current.data, sex.sig.list){
  
  quantos<- table(current.data$info$Sexo, useNA = "always") 
  
  i.femeas <- if (quantos[1] >3)  current.data$info$Sexo == "F" | is.na(current.data$info$Sexo) else current.data$info$Sexo == "nenhum"
  femeas <- as.matrix(current.data$ed.raw)[i.femeas,]
  mx.fm <- cov(femeas)
  
  i.machos <- if (quantos[2] >3)  current.data$info$Sexo == "M" | is.na(current.data$info$Sexo) else current.data$info$Sexo == "nenhum"
  machos <- as.matrix(current.data$ed.raw)[i.machos,]
  mx.mac <- cov(machos)
  
  lista.sex <- list("Males" = mx.mac,
                    "Females" = mx.fm,
                    "Sex Fit" = sex.sig.list$cov.fit, 
                    "Raw Data" =  current.data$matrix$cov)
  
  RS <- RandomSkewers(lista.sex)$correlations
  KRZ<- KrzCor(lista.sex, ret.dim = 16)
  MMstat<- ldply(lista.sex, MeanMatrixStatistics)
  
  results <- list("Quantos" = quantos,
                  "RS" = RS,
                  "KRZ" =  KRZ,
                  "MMstat" = MMstat)
  return(results)
  
}

tabele.quantos.sexo <- sp.main.data %>% ldply(function (x) table(x$info$Sexo, useNA = "always") ) 
tabele.quantos.sexo[,1] %<>% gsub("_", " ", .)
rownames(tabele.quantos.sexo) <-tabele.quantos.sexo[,1]
tabele.quantos.sexo[,-1] %>% xtable

mask.sex.sample.size <- tabele.quantos.sexo[,2] != 0 & tabele.quantos.sexo[,3] != 0

sex.MX.compare <- SexCompare(sp.main.data[[34]], sex.sig[[34]])
Microcebus.sex.comparison <- SexCompare(sp.main.data$Microcebus_griseorufus, sex.sig$Microcebus_griseorufus)
Combine.Mx.Plot(Mx1 = t(Microcebus.sex.comparison$RS), Mx2 = t(Microcebus.sex.comparison$KRZ), diag.info = c(32,34, 68, 68), titulo = "Covariance matrices comparisons via KRZ and RS")

mat_data <- t(Microcebus.sex.comparison$RS)
mat_data[lower.tri(mat_data)] <- t(t(Microcebus.sex.comparison$KRZ))[lower.tri(t(Microcebus.sex.comparison$KRZ))]
diag(mat_data)<- NA
range.values<- range(mat_data, na.rm = T) - c(0.01, -0.01)
diag.info = c(32,34, 68, 68)
diag(mat_data) <- diag.info

mixed.mx = melt(mat_data) 
mixed.mx.position =  mixed.mx
mixed.mx.position$value= round( mixed.mx.position$value, 2)

myPalette <- colorRampPalette(rev(brewer.pal(7, 'Spectral')), space = 'Lab')(n = 10)
mixed.mx.cute.plot <- 
  ggplot (mixed.mx.position) +
  geom_tile(aes(x = Var2, y = Var1, fill = value)) +
  scale_fill_gradientn(name = '', colours = myPalette, limits = range.values, na.value = "white") +
  ylab ('') + xlab ('') + labs(title = "Covariance matrices comparisons via KRZ and RS") + theme(plot.title = element_text(face = "bold", size = 50)) +
  geom_text(aes(x = Var2, y = Var1, label = value), size = 10) +
  scale_y_discrete(limits = rev(levels(mixed.mx.position$Var1)) ) +
  # scale_x_discrete() +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle =0, hjust = 0.5, face = 'italic', size =15),
        axis.text.y = element_text(face = "italic", size =20),
        axis.ticks = element_line(size = 0),
        #legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        rect = element_blank(), line = element_blank())



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

quantos.sexo <- sp.main.data %>% ldply(function (x) table(x$info$Sexo, useNA = "always") ) 
mask.sex <- quantos.sexo[,2] >=2 & quantos.sexo[,3] >= 2 
sp.main.data[mask.sex] %>% llply(function(x) t.test(x$gm.ind[!is.na(x$info$Sexo)] ~ x$info$Sexo[!is.na(x$info$Sexo)]) ) %>% ldply( function (x) x$estimate  ) 
sp.main.data[mask.sex] %>% llply(function(x) t.test(x$gm.ind[!is.na(x$info$Sexo)] ~ x$info$Sexo[!is.na(x$info$Sexo)]) ) %>% llply( function (x) x$p.value ) %>% ldply (function (x) x >=0.05 ) %>% .[,2] %>% table


gm.plot <-  sp.main.data[mask] %>% ldply(function(x) cbind(as.character(x$info$Sexo), t(x$gm.ind) ) ) 
names(gm.plot) <- c(".sp", "sex", "gm")
gm.plot$sex %<>% {gsub("M", "Male", .)} %>% {gsub("F", "Female", .)} 
gm.plot$sex[is.na(gm.plot$sex)] <- "Undefined"
gm.plot$gm <- as.character(gm.plot$gm)
gm.plot$gm <- as.numeric(gm.plot$gm)
gm.plot$gm.log <- log(gm.plot$gm)
gm.plot$gm <- round(gm.plot$gm, digits = 2)
gm.plot$.sp <- factor(gm.plot$.sp, levels = unique(gm.plot$.sp)[1:42])



ggplot(gm.plot, aes(x = sex, y = gm, group = interaction(.sp, sex), fill = sex)) + 
  geom_violin() + 
  facet_wrap(~.sp) + 
 # background_grid(major = 'y', minor = "none") + 
  labs(y = "Geometric mean of cranial traits", x = "geometric mean")



gm.plot %>% 
  ggplot(.,aes(y = log(gm), x = sex  )) +
  geom_violin(aes(group = sex, fill = sex, color = sex) ,alpha= 0.1) +
  geom_jitter(aes(group = sex, fill = sex, color = sex)) +
  #scale_fill_manual(values = c("red", "grey")) +
  #scale_color_manual(values = c("red", "grey"))  +
  #geom_text(aes(group = "Trait", label = Specie)) +
  facet_wrap(~.sp, scale="free") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=17),
        legend.text= element_text(size=10), 
        plot.title = element_text(lineheight=.8, face="bold", size = 8)) + labs(title = "Cranial traits loged geometric mean")
  




cvs = data.frame(t(laply(sp.main.data[mask], function(x) sqrt(diag(x$matrix$cov))/ x$ed.means)))
names(cvs) <- names(sp.main.data[mask])
names(cvs) %<>% gsub("_", " ", . )
cvs$traits = factor(rownames(cvs), levels = rownames(cvs)[1:39])
cvs %<>% melt 
cvs$variable <- factor(cvs$variable, levels = unique(cvs$variable))
cvs$traits <- factor(cvs$traits, levels= rev(unique(cvs$traits)) )
cvs$cagated <- cvs$value < 0.21

colore <- c("#FF9326",  "#9999FF", "#D92121", "#21D921", "#FFFF4D", "#2121D9", "white")

cv.plot <- cvs %>%  
  ggplot(aes(x = traits, y = value)) +
  #geom_violin() +
  geom_jitter(alpha = 0.8,color = "black", size = 2, shape = 4) +
  geom_text(aes(label = variable, alpha = value - 0.1, colour = factor(cagated) ), fontface = "italic" ) +
  coord_flip() +
  scale_colour_manual (values= c("red","darkgrey") ) +
  theme_bw() + 
  #facet_wrap(~variable, nrow = 10, ncol = 4) +
  theme(legend.position = "none") +
  labs(x = "", y = "Coeficient of variation") 
  

REP$Trait <- factor(REP$Trait, rev(unique(REP$Trait)))

rep.plot <- REP %>% na.omit %>% group_by(Trait)  %>%
  ggplot(.) +
  geom_violin(aes(x = Trait, y =  Repetability), alpha = 0.4) +
  geom_jitter(aes(x = Trait, y =  Repetability), alpha = 0.8,color = "black", size = 2, shape = 4) +
  geom_text(aes(x = Trait, y =  Repetability, label = Specie, alpha = abs(1.1-Repetability)), color = "darkgrey", size = 4) +
  coord_flip() +
  theme_bw() +
  scale_colour_brewer(palette = "Set1", direction = -1) +
  xlab("")+
  #ylab( "Traits repetabilities")+
  #labs(title = "Traits repetabilities") +
  theme(legend.position="none",
        axis.text.y = element_text( size =10),
        axis.text.x = element_text(size =19),
        strip.text= element_text(size=19),
        plot.title = element_text(face = "bold", size = 20)) 

plot_grid(rep.plot, cv.plot, labels = LETTERS[1:2] )

