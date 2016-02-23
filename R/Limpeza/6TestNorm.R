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

apply(sig.p.uni.norm,2, FUN = function(x) table(x == TRUE))
quantos.normais <- apply(sig.p.uni.norm,1, FUN = function(x) table(x == TRUE))
names(quantos.normais) <- p.uni.norm$Species



RawData<- raw.data %>% tbl_df()
RawData[, 1:54] <- filter(., Genero == "Microcebus") as.data.frame() %>% ggscatmat(., columns = 16:54, color = "Especie" )
RawData[, 1:54] %>% filter(., Familia == "Cheirogaleidae") %>% as.data.frame() %>% ggscatmat(., columns = 16:54, color = "Especie" )

RawData[, 1:54] %>% filter(., Especie == "Euoticus_elegantulus") %>% as.data.frame() %>% ggscatmat(., columns = 16:54 )
########  graficos de densidade 
RawData[, 1:54] %>% filter(., Especie == "Euoticus_elegantulus") %>% gather(key=ed, value=value, (16:54) ) %>% 
  ggplot(., aes(x=value) ) +
  geom_density(aes(group = ed, color = ed)) +
  facet_wrap(~ed,  scales="free")
  