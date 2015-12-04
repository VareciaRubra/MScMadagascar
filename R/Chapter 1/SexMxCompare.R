SexCompare <- function (current.data, )

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



