if(!require(mvoutlier)) {install.packages('mvoutlier'); library(mvoutlier)}
if(!require(MVN)) {install.packages('MVN'); library(MVN)} # achar outliers de uma distribuição multivariada normal usando a mvOutlier {MVN}
if(!require(nortest)) {install.packages('nortest'); library(nortest)} # testar normalidade univariados por Kolmogorov-Smirnov usando a lillie.test {nortest}
result.outlier<- mvoutlier.CoDa(x = sp.main.data$Microcebus_griseorufus$ed,quan = 0.75)
plot.mvoutlierCoDa(x = result.outlier, onlyout=TRUE,bw=FALSE,which="parallel",symb=TRUE,symbtxt=TRUE,transp=0.3)

result.mv.out<- mvOutlier(data = sp.main.data$Microcebus_griseorufus$ed, qqplot = T, method = 'adj.quan', label = sp.main.data$Microcebus_griseorufus$info$Tombo)
mvOutlier(data = sp.main.data$Microcebus_griseorufus$ed[!(result.outlier$outliers),], qqplot = T, method = 'adj.quan')

dif <- colMeans(sp.main.data$Microcebus_griseorufus$ed[(as.logical(result.mv.out$outlier$Outlier)),]) - colMeans(sp.main.data$Microcebus_griseorufus$ed[!(as.logical(result.mv.out$outlier$Outlier)),])
dif / apply(sp.main.data$Microcebus_griseorufus$ed, 2, sd)
result.outlier$outliers

RandomSkewers(cov.list, cov(sp.main.data$Microcebus_griseorufus$ed[!(as.logical(result.mv.out$outlier$Outlier)),]))

roystonTest(data = sp.main.data$Microcebus_griseorufus$ed[!(result.outlier$outliers),], qqplot = T)
mar.test <- mardiaTest(data = sp.main.data$Microcebus_griseorufus$ed[!(result.outlier$outliers),], qqplot = T, cov = T)

lillie.test(x = sp.main.data$Microcebus_griseorufus$ed)

cov.mx <- current.data %>% llply(function(x) x$matrix$cov)
eds.sp <- sp.main.data %>% llply(function(x) x$ed)

mask <- ldply(cov.mx, function(x) !is.na(x[1]))[,2]eds.sp [mask]

p.uni.norm <- eds.sp[mask] %>% ldply(function (x) apply(x, 2, FUN = function (x) round(lillie.test(x)$p.value, 4 )) )
names (p.uni.norm)[1] <- "Species"
sig.p.uni.norm <- p.uni.norm[,-1] >= 0.05


Uni.Normal.Density.Test <- function (RawSpList){
  ## RawSpList =  Lista com o rawdata frame por especie
  eds <- dplyr::select(RawSpList, c(IS_PM:BA_OPI)) # pegando as ed de cada replica
  gm <- apply(eds, 1, FUN = function(x) exp( mean( as.matrix(log (x)) ) )  ) 
  eds$GM <- gm
  eds.log <- log(eds)
  
  
  uni.normal.c <- eds %>% apply(. , 2, FUN = function (x) round(lillie.test(x)$p.value, 4 ))
  uni.normal.l <- eds.log %>% apply(. , 2, FUN = function (x) round(lillie.test(x)$p.value, 4 ))
  uni.normal<- uni.normal.c<=0.05
  uni.normal.log <- uni.normal.l<=0.05
  uni.normal[uni.normal == TRUE] <- "p < 0.05" 
  uni.normal[uni.normal == FALSE] <- "p > 0.05" 
  uni.normal <- as.character(as.data.frame(uni.normal)[,1])
  
  uni.normal.log[uni.normal.log == TRUE] <- "p < 0.05" 
  uni.normal.log[uni.normal.log == FALSE] <- "p > 0.05" 
  uni.normal.log <- as.character(as.data.frame(uni.normal.log)[,1])
  
  RawSpList$GM <- apply(RawSpList[16:54], 1, FUN = function(x) exp( mean( as.matrix(log (x)) ) )  ) 
  RawSpList$Sexo <-  as.character(RawSpList$Sexo)
  RawSpList$Sexo[is.na(RawSpList$Sexo)] <- "U"
  RawSpList$Sexo <-  as.factor(RawSpList$Sexo)
  
  RawSpList$Planilha <- as.factor(RawSpList$Planilha)
  Density <- RawSpList[,c(1:54, 63)]  %>% melt() 
  names(Density)[16] <- "Trait" 
  Density.log <- Density
  Density.log$Trait <- log(Density.log$value)

  Univariate_Normality <- rep(uni.normal, each =  dim(RawSpList)[1]) 
  Univariate_Normality.log<- rep(uni.normal.log, each =  dim(RawSpList)[1]) 
  
  specie <- as.character(unique(RawSpList$Especie))

  Plot <- Density%>% 
    ggplot(.,aes(x = Sexo, y = value)) +
    geom_violin(aes(group = interaction(Trait, Sexo), fill =Univariate_Normality, color = Univariate_Normality) ,alpha= 0.1) +
    scale_fill_manual(values = c("red", "grey")) +
    scale_color_manual(values = c("red", "grey"))  +
    geom_jitter(aes(group = interaction(Trait, Sexo), shape= Sexo), alpha = 0.3 ) +
    #geom_text(aes(group = "Trait", label = Specie)) +
    facet_wrap(~Trait, scale="free", ncol =5, nrow = 8) +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 7),
          #axis.title.x = element_blank(),
          #axis.title.y = element_text(size=17),
          legend.text= element_text(size=10), 
          plot.title = element_text(lineheight=.8, face="bold", size = 8),
          title = element_text(specie),
          legend.position = "bottom") + labs(title = paste(specie, "traits distribuition") )
  
  Plot.log <- Density.log%>% 
    ggplot(.,aes(value)) +
    geom_density(aes(group = Trait, fill = Univariate_Normality.log, color = Univariate_Normality.log) ,alpha= 0.1) +
    scale_fill_manual(values = c("red", "grey")) +
    scale_color_manual(values = c("red", "grey"))  +
    #geom_text(aes(group = "Trait", label = Specie)) +
    facet_wrap(~Trait, scale="free", ncol =5, nrow = 8) +
    theme_bw() +
    ggtitle(paste(c(specie, "loged values")) ) +
    theme(axis.text.x = element_text(size =8), 
          axis.text.y = element_text(size = 7),
          #axis.title.x = element_blank(),
          #axis.title.y = element_text(size=17),
          legend.text= element_text(size=10), 
          plot.title = element_text(lineheight=.8, face="bold", size = 8)) + 
    labs(title = paste(specie, "log ed values") )
  
  return(list("Plot" = Plot,
              "UniTestLillie" = uni.normal.c,
              "Plot.log" = Plot.log ))
}

test <- Uni.Normal.Density.Test(RawSpList = Sp.raw.main.data$Varecia_rubra)
test$Plot
test$Plot.log
Uni.Tests <- llply(Sp.raw.main.data[mask], .fun = Uni.Normal.Density.Test)
Uni.Tests$Microcebus_griseorufus$Plot



euoticus.uni.normal<- eds.sp$Euoticus_elegantulus %>% apply(. , 2, FUN = function (x) round(lillie.test(x)$p.value, 4 ))
euoticus.uni.normal<- euoticus.uni.normal<=0.05
euoticus.uni.normal[euoticus.uni.normal == TRUE] <- "p.value < 0.05" 
euoticus.uni.normal[euoticus.uni.normal == FALSE] <- "p.value > 0.05"
euoticus.uni.normal <- as.character(as.data.frame(euoticus.uni.normal)[,1])

RawData$Planilha <- as.factor(RawData$Planilha)
DensityEuoticus <- RawData[, 1:54] %>% filter(., Especie == "Euoticus_elegantulus") %>% melt() 
names(DensityEuoticus)[16] <- "Trait" 

dim(RawData[, 1:54] %>% filter(., Especie == "Euoticus_elegantulus"))

col.uni.normal <- rep(euoticus.uni.normal, each =  30) 

DensityEuoticus %>% 
  ggplot(.,aes(value)) +
  geom_density(aes(group = Trait, fill = col.uni.normal, color = col.uni.normal) ,alpha= 0.1) +
  scale_fill_manual(values = c("red", "grey")) +
  scale_color_manual(values = c("red", "grey")) +
  #geom_text(aes(group = "Trait", label = Specie)) +
  facet_wrap(~Trait, scale="free", ncol =5, nrow = 8) +
  theme_bw() +
  theme(axis.text.x = element_text(size =8), 
        axis.text.y = element_text(size = 7),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(size=17),
        legend.text= element_text(size=10)) +
  ggtitle("Euoticus elegantulus Univariate distribuition") + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = 8)) 


DensityGM <- sp.main.data %>% llply(function(x) x$gm.ind ) 

DensityGmPlot <- function (x)
  { x <- x %>% melt()
    p.normal.test<- round(lillie.test(x$value)$p.value, 4 )
 Plot <- x %>%
    ggplot(.,aes(value)) +
    geom_density() +
    theme_bw() +
    theme(axis.text.x = element_text(size =15), 
          axis.text.y = element_text(size = 15),
          #axis.title.x = element_blank(),
          #axis.title.y = element_text(size=17),
          legend.text= element_text(size=10)) +
    theme(plot.title = element_text(lineheight=.8, face="bold", size = 8)) 
 
 return (list("Plot" = Plot,
              "GM.sd" = sd(x$value),
              "GM.mean" = mean(x$value),
              "p.lillie.test" = p.normal.test) )
}

GMDensityPlots <- llply(DensityGM[mask], .fun = DensityGmPlot)
