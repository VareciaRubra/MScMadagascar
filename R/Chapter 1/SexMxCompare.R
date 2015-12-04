i.machos <- sp.main.data$Microcebus_griseorufus$info$Sexo == "M" | is.na(sp.main.data$Microcebus_griseorufus$info$Sexo)
i.femeas <- sp.main.data$Microcebus_griseorufus$info$Sexo == "F" | is.na(sp.main.data$Microcebus_griseorufus$info$Sexo)
machos <- as.matrix(sp.main.data$Microcebus_griseorufus$ed.raw)[i.machos,]
femeas <- as.matrix(sp.main.data$Microcebus_griseorufus$ed.raw)[i.femeas,]
regular<- sp.main.data$Microcebus_griseorufus$matrix$cov
residuals <- sex.sig$Microcebus_griseorufus

mx.mac <- cov(machos)
mx.fm <- cov(femeas)

lista.sex.microcebus <- list("machos" = mx.mac,
                             "femeas" = mx.fm)
RandomSkewers(lista.sex.microcebus)
KrzProjection(lista.sex.microcebus)
llply(lista.sex.microcebus, CalcR2)
llply(lista.sex.microcebus, CalcICV)
llply(lista.sex.microcebus, Pc1Percent)
llply(lista.sex.microcebus, function (x) mean(Flexibility(x)) )
llply(lista.sex.microcebus, function (x) mean(Evolvability(x)) )
            