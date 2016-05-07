SRD.selected.Diet.Lemuridae.plot <- plot_grid(PC.Plots$Varecia[[7]], ggplot() ,ggplot() ,ggplot() ,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Hapalemur),hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              PC.Plots$Hapalemur[[7]], ggplot(), ggplot(),
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Lemur),varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lemur, gen.cov.list$Hapalemur),hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              PC.Plots$Lemur[[7]], ggplot(), 
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Eulemur),eulemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Hapalemur, gen.cov.list$Eulemur),lemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lemur, gen.cov.list$Eulemur),eulemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              PC.Plots$Eulemur[[7]],
                                              ncol = 4, scale = c(1.2,1.3,1.3,1.3,
                                                                  1.3,1.2,1.3,1.3,
                                                                  1.3,1.3,1.2,1.3,
                                                                  1.3,1.3,1.3,1.2))
SRD.selected.Diet.Lemuridae.plot

SRD.selected.Diet.Indridae.plot  <- plot_grid(PC.Plots$Indri[[7]], ggplot() ,ggplot(),
          SRD.plot.wire(SRD.result = SRD(gen.cov.list$Indri, gen.cov.list$Propithecus),indri, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
          PC.Plots$Propithecus[[7]] , ggplot(), 
          SRD.plot.wire(SRD.result = SRD(gen.cov.list$Indri, gen.cov.list$Avahi),avahi, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
          SRD.plot.wire(SRD.result = SRD(gen.cov.list$Avahi, gen.cov.list$Propithecus),propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
          PC.Plots$Avahi[[7]],
          ncol = 3)
SRD.selected.Diet.Indridae.plot


SRD.selected.Lepilemuridae.plot  <- plot_grid(SRD.plot.wire(SRD.result = SRD(Drift.results$extant.sp$Regression.test$`96`$W.mx, Drift.results$extant.sp$Regression.test$`89`$W.mx),lepilemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev, 
                                              PC.Plots$Lepilemur[[7]], ncol = 2, scale = c(1.1,1.2))



SRD.selected.Diet.Cheirogaleidae.plot <- plot_grid(PC.Plots$Phaner[[7]], ggplot() ,ggplot() ,ggplot() ,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Cheirogaleus),cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              PC.Plots$Cheirogaleus[[7]], ggplot(), ggplot(),
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Microcebus),microcebus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Cheirogaleus),cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              PC.Plots$Microcebus[[7]], ggplot(), 
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Mirza),phaner, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Cheirogaleus, gen.cov.list$Mirza),cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
                                              SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Mirza),mirza, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
                                              PC.Plots$Mirza[[7]],
                                              ncol = 4, scale = c(1,1.3,1.3,1.3,
                                                                  1.3,1,1.3,1.3,
                                                                  1.3,1.3,1,1.3,
                                                                  1.3,1.3,1.3,1))


plot_grid(PC.Plots$Daubentonia[[7]], ggplot(), ggplot(), ggplot(), ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Daubentonia),daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Varecia[[7]], ggplot(), ggplot(), ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Hapalemur, gen.cov.list$Daubentonia),hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Hapalemur, gen.cov.list$Varecia),varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Hapalemur[[7]], ggplot(), ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Propithecus, gen.cov.list$Daubentonia),propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Propithecus, gen.cov.list$Varecia),propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Propithecus, gen.cov.list$Hapalemur),propithecus, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Propithecus[[7]] , ggplot(), ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lepilemur, gen.cov.list$Daubentonia),daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Lepilemur),varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lepilemur, gen.cov.list$Hapalemur),hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Lepilemur, gen.cov.list$Propithecus),lepilemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Lepilemur[[7]] , ggplot(),ggplot(),
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Daubentonia),daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Microcebus),varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Hapalemur),hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Propithecus),microcebus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Microcebus, gen.cov.list$Lepilemur),microcebus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Microcebus[[7]] , ggplot(), 
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Daubentonia),daubentonia, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Varecia, gen.cov.list$Phaner),varecia, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Hapalemur),hapalemur, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Propithecus),cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Lepilemur),phaner, ROTACIONI =  c(-1,-1,1), TTL = "")$plot.muDev,
SRD.plot.wire(SRD.result = SRD(gen.cov.list$Phaner, gen.cov.list$Microcebus),cheirogaleus, ROTACIONI =  c(1,-1,1), TTL = "")$plot.muDev,
PC.Plots$Phaner[[7]],
ncol = 7, scale = c(1.2,1.3,1.3,1.3,1.3,1.3,1.3,
                   1.3,1.2,1.3,1.3,1.3,1.3,1.3,
                   1.3,1.3,1.2,1.3,1.3,1.3,1.3,
                   1.3,1.3,1.3,1.2,1.3,1.3,1.3,
                   1.3,1.3,1.3,1.3,1.2,1.3,1.3,
                   1.3,1.3,1.3,1.3,1.3,1.2,1.3,
                   1.3,1.3,1.3,1.3,1.3,1.3,1.2))


################ refazer o PCs do Microcebus, que tá com -
################ faz algum sentido fazer K multi nos valores de SRD?

          SRD_Diet_list <- list (SRD(gen.cov.list$Varecia, gen.cov.list$Daubentonia), 
          SRD(gen.cov.list$Hapalemur, gen.cov.list$Daubentonia), 
          SRD(gen.cov.list$Hapalemur, gen.cov.list$Varecia),
          SRD(gen.cov.list$Propithecus, gen.cov.list$Daubentonia), 
          SRD(gen.cov.list$Propithecus, gen.cov.list$Varecia), 
          SRD(gen.cov.list$Propithecus, gen.cov.list$Hapalemur), 
          SRD(gen.cov.list$Lepilemur, gen.cov.list$Daubentonia), 
          SRD(gen.cov.list$Varecia, gen.cov.list$Lepilemur), 
          SRD(gen.cov.list$Lepilemur, gen.cov.list$Hapalemur), 
          SRD(gen.cov.list$Lepilemur, gen.cov.list$Propithecus), 
          SRD(gen.cov.list$Microcebus, gen.cov.list$Daubentonia), 
          SRD(gen.cov.list$Varecia, gen.cov.list$Microcebus), 
          SRD(gen.cov.list$Microcebus, gen.cov.list$Hapalemur), 
          SRD(gen.cov.list$Microcebus, gen.cov.list$Propithecus), 
          SRD(gen.cov.list$Microcebus, gen.cov.list$Lepilemur), 
          SRD(gen.cov.list$Phaner, gen.cov.list$Daubentonia), 
          SRD(gen.cov.list$Varecia, gen.cov.list$Phaner), 
          SRD(gen.cov.list$Phaner, gen.cov.list$Hapalemur), 
          SRD(gen.cov.list$Phaner, gen.cov.list$Propithecus), 
          SRD(gen.cov.list$Phaner, gen.cov.list$Lepilemur), 
          SRD(gen.cov.list$Phaner, gen.cov.list$Microcebus))





  