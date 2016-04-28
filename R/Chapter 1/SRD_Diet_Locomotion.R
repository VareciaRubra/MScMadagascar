SRD.selected.Diet.Lemuridae.plot <- plot_grid(PC.Plots$Varecia[[7]] +ggtitle("PC1 Varecia"), ggplot() ,ggplot() ,ggplot() ,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Hapalemur), SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "Varecia x Hapalemur")$plot.muDev,
                                              PC.Plots$Hapalemur[[7]] +ggtitle("PC1 Hapalemur"), ggplot(), ggplot(),
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Lemur), SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "Varecia x Lemur")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lemur, gen.cov.list$Hapalemur), SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "Lemur x Hapalemur")$plot.muDev,
                                              PC.Plots$Lemur[[7]]+ggtitle("PC1 Lemur"), ggplot(), 
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Eulemur), SHAPE = Shapes.sym$eulemur, ROTACIONI =  c(1,-1,1), TTL = "Varecia x Eulemur")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Hapalemur, gen.cov.list$Eulemur), SHAPE = Shapes.sym$lemur, ROTACIONI =  c(1,-1,1), TTL = "Hapalemur x Eulemur")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lemur, gen.cov.list$Eulemur), SHAPE = Shapes.sym$lemur, ROTACIONI =  c(1,-1,1), TTL = "Lemur x Eulemur")$plot.muDev,
                                              PC.Plots$Eulemur[[7]] + ggtitle("PC1 Eulemur"),
                                              ncol = 4)
SRD.selected.Diet.Lemuridae.plot

SRD.selected.Diet.Indridae.plot  <- plot_grid(PC.Plots$Indri[[7]] + ggtitle("PC1 Indri"), ggplot() ,ggplot(),
          SRD.plot.wire(SRD.result = SRD(gen.cov.list$Indri, gen.cov.list$Propithecus), SHAPE = Shapes.sym$indri, ROTACIONI =  c(-1,-1,1), TTL = "Indri x Propithecus")$plot.muDev,
          PC.Plots$Propithecus[[7]] +ggtitle("PC1 Propithecus"), ggplot(), 
          SRD.plot.wire(SRD.result = SRD(gen.cov.list$Indri, gen.cov.list$Avahi), SHAPE = Shapes.sym$avahi, ROTACIONI =  c(-1,-1,1), TTL = "Indri x Avahi")$plot.muDev,
          SRD.plot.wire(SRD.result = SRD(gen.cov.list$Avahi, gen.cov.list$Propithecus), SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "Avahi x Propithecus")$plot.muDev,
          PC.Plots$Avahi[[7]]+ggtitle("PC1 Avahi"),
          ncol = 3)
SRD.selected.Diet.Indridae.plot


SRD.selected.Lepilemuridae.plot  <- plot_grid(SRD.plot.wire(SRD.result = SRD(Drift.results$extant.sp$Regression.test$`96`$W.mx, Drift.results$extant.sp$Regression.test$`89`$W.mx), SHAPE = Shapes.sym$lepilemur, ROTACIONI =  c(1,-1,1), TTL = "Lepilemuridae")$plot.muDev, 
                                              PC.Plots$Lepilemur[[7]], ncol = 2)



plot_grid(PC.Plots$Daubentonia[[7]], ggplot(), ggplot(), ggplot(), ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Daubentonia), SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Varecia[[7]], ggplot(), ggplot(), ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Hapalemur, gen.cov.list$Daubentonia), SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Hapalemur, gen.cov.list$Varecia), SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Hapalemur[[7]], ggplot(), ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Propithecus, gen.cov.list$Daubentonia), SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Propithecus, gen.cov.list$Varecia), SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Propithecus, gen.cov.list$Hapalemur), SHAPE = Shapes.sym$propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Propithecus[[7]] , ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lepilemur, gen.cov.list$Daubentonia), SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Lepilemur), SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lepilemur, gen.cov.list$Hapalemur), SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lepilemur, gen.cov.list$Propithecus), SHAPE = Shapes.sym$lepilemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Lepilemur[[7]] , ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Daubentonia), SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Microcebus), SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Hapalemur), SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Propithecus), SHAPE = Shapes.sym$microcebus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Lepilemur), SHAPE = Shapes.sym$microcebus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Microcebus[[7]] , ggplot(), 
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Daubentonia), SHAPE = Shapes.sym$daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Phaner), SHAPE = Shapes.sym$varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Hapalemur), SHAPE = Shapes.sym$hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Propithecus), SHAPE = Shapes.sym$cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Lepilemur), SHAPE = Shapes.sym$phaner, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Microcebus), SHAPE = Shapes.sym$cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Phaner[[7]],
ncol = 7, scale = c(1.2,1.3,1.3,1.3,1.3,1.3,1.3,
                   1.3,1.2,1.3,1.3,1.3,1.3,1.3,
                   1.3,1.3,1.2,1.3,1.3,1.3,1.3,
                   1.3,1.3,1.3,1.2,1.3,1.3,1.3,
                   1.3,1.3,1.3,1.3,1.2,1.3,1.3,
                   1.3,1.3,1.3,1.3,1.3,1.2,1.3,
                   1.3,1.3,1.3,1.3,1.3,1.3,1.2))









  