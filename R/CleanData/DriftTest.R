

Wmat.All.Fuckers = CalculateMatrix(manova(as.matrix(all.main.data$All$ed)  ~ Especie, data = as.data.frame(all.main.data$All$info) ) )
Wmat.extant.Fuckers = CalculateMatrix(manova(as.matrix(extant.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.main.data$All$info) ) )
Wmat.madagascar.Fuckers = CalculateMatrix(manova(as.matrix(madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(madagascar.main.data$All$info) ) )
Wmat.extant.madagascar.Fuckers = CalculateMatrix(manova(as.matrix(extant.madagascar.main.data$All$ed)  ~ Especie, data = as.data.frame(extant.madagascar.main.data$All$info) ) )

mask.extant<- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Status != "Extinct"])
mask.madagascar<- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Regiao == "Madagascar"])
mask.madagascar.extant <- unique(All.raw.main.data$All$Especie[All.raw.main.data$All$Status != "Extinct" &All.raw.main.data$All$Regiao == "Madagascar"])

means.All.Fuckers<- sp.master.main.data %>% llply(function(x) x$ed.means )
means.extant.Fuckers<- means.All.Fuckers[mask.extant] 
means.madagascar.Fuckers<- means.All.Fuckers[mask.madagascar]
means.extant.madagascar.Fuckers<- means.All.Fuckers[mask.madagascar.extant]
  
DriftTest(means = means.All.Fuckers, cov.matrix = Wmat.All.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.Fuckers, cov.matrix = Wmat.extant.Fuckers, show.plot = TRUE)
DriftTest(means = means.madagascar.Fuckers, cov.matrix = Wmat.madagascar.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.madagascar.Fuckers, cov.matrix = Wmat.extant.madagascar.Fuckers, show.plot = TRUE)
DriftTest(means = means.extant.madagascar.Fuckers, cov.matrix = Wmat.madagascar.Fuckers, show.plot = TRUE)
