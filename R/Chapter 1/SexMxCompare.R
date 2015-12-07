SexCompare <- function (current.data){
  
  i.machos <- current.data$info$Sexo == "M" | is.na(current.data$info$Sexo)
  i.femeas <- current.data$info$Sexo == "F" | is.na(current.data$info$Sexo)
  machos <- as.matrix(current.data$ed.raw)[i.machos,]
  femeas <- as.matrix(current.data$ed.raw)[i.femeas,]
  regular<- current.data$matrix$cov
  residuals <- sex.sig$Microcebus_griseorufus
  
  mx.mac <- cov(machos)
  mx.fm <- cov(femeas)
  
  lista.sex <- list("machos" = mx.mac,
                               "femeas" = mx.fm,
                               "sex.fit" = sex.sig.$Microcebus_griseorufus$cov.fit,
                               "raw.cov" =  current.data$matrix$cov)
  RS <- RandomSkewers(lista.sex)$correlations
  KRZ<- KrzProjection(lista.sex)
  MMstat<- ldply(lista.sex, MeanMatrixStatistics)
  
  results <- list("RS" = RS,
                  "KRZ" =  KRZ,
                  "MMstat" = MMstat)
  return(results)
  
}


i.machos <- sp.main.data$Microcebus_griseorufus$info$Sexo == "M" | is.na(sp.main.data$Microcebus_griseorufus$info$Sexo)
i.femeas <- sp.main.data$Microcebus_griseorufus$info$Sexo == "F" | is.na(sp.main.data$Microcebus_griseorufus$info$Sexo)
machos <- as.matrix(sp.main.data$Microcebus_griseorufus$ed.raw)[i.machos,]
femeas <- as.matrix(sp.main.data$Microcebus_griseorufus$ed.raw)[i.femeas,]
regular<- sp.main.data$Microcebus_griseorufus$matrix$cov
residuals <- sex.sig$Microcebus_griseorufus

mx.mac <- cov(machos)
mx.fm <- cov(femeas)

lista.sex.microcebus <- list("machos" = mx.mac,
                             "femeas" = mx.fm,
                             "sex.fit" = sex.sig.$Microcebus_griseorufus$cov.fit,
                             "raw.cov" =  sex.sig.$Microcebus_griseorufus$cov)
RandomSkewers(lista.sex.microcebus)$correlations
KrzProjection(lista.sex.microcebus)
ldply(lista.sex.microcebus, MeanMatrixStatistics)



