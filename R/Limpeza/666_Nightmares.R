 apply(as.matrix(sp.main.data$Euoticus_elegantulus$ed ), 2, function (x) return (Anova(lm(x ~ Museu, data = as.data.frame(sp.main.data$Euoticus_elegantulus$info) ), 
                                                                 type= 3, test.statistic="Wilks",icontrasts=c("contr.sum", "contr.poly" ))) )
correct.euoticus <- manova(as.matrix(sp.main.data$Euoticus_elegantulus$ed) ~ Sexo, data = as.data.frame(sp.main.data$Euoticus_elegantulus$info) )

correct.euoticus.mx <- CalculateMatrix(correct.euoticus)

max(RandomSkewers(cov.x = cov.list, cov.y = correct.euoticus.mx, num.vectors = 100)$correlation)


########## Microcebus griseorrufus
current.data <- Sp.raw.main.data$Microcebus_griseorufus
current.data %>% 
  filter (Tombo != "174522") %>%
  filter (Tombo != "174496") %>%
  filter (Tombo != "174515") %>%
  filter (Tombo != "174434")
new.microcebus.griseorrufus <- llply(current.data, specie = TRUE, compare.size = TRUE, final = TRUE, makeMainData, .progress = progress_text(char = "."), .inform = T)

  removed <- c("174522", "174496", "174515", "174434")

removed  <- (sp.main.data$Microcebus_griseorufus$info$Tombo %in% removed) == FALSE
temp.microcebus <- cov(sp.main.data$Microcebus_griseorufus$ed[removed,])


as.matrix(eigen (temp.microcebus) $ values/sum (eigen (temp.microcebus) $ values) )
as.matrix(eigen (sp.main.data$Microcebus_griseorufus$matrix$cov) $ values/sum (eigen (sp.main.data$Microcebus_griseorufus$matrix$cov) $ values) )

mean(RandomSkewers(cov.x = cov.list[names(cov.list) != "Microcebus_griseorufus"], cov.y =  temp.microcebus, num.vectors = 100)$correlation)
temp.rar.rs <- Rarefaction(ind.data = sp.main.data$Microcebus_griseorufus$ed[removed,], num.reps = 64, ComparisonFunc = RandomSkewers, parallel = F)
temp.rar.krz <- Rarefaction(ind.data = sp.main.data$Microcebus_griseorufus$ed[removed,], num.reps = 64, ComparisonFunc = KrzCor, parallel = F)

dim(sp.main.data$Microcebus_griseorufus$ed)
