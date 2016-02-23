############# Rarefação no calculo das matrizes 
# pegando as maiores amostras de duas espécies nos extremos da variaçao de tamanho : Microcebus griseorrufus e Indri indri
# por dois métodos de comparação RS e KRZ 
# tudo ta calculado já no sp.main.data

p1 <- PlotRarefaction(sp.main.data$Microcebus_griseorufus$rarefaction$krz) + 
  labs(title = "Microcebus griseorufus") + 
  theme(plot.title = element_text(face = "italic", size = 20), axis.title.x = element_blank(), axis.text.x =element_blank(), axis.ticks.x = element_line(size =0), axis.text= element_text(size = 15)) + 
  coord_cartesian(ylim=c(0.4, 1)) + scale_y_continuous("KRZ", breaks=seq(0.4, 1, 0.2))
p2 <- PlotRarefaction(sp.main.data$Indri_indri$rarefaction$krz) + 
  labs(title = "Indri indri") + 
  theme(plot.title = element_text(face = "italic", size = 20), axis.title = element_blank(), axis.text =element_blank(), axis.ticks = element_line(size =0)) +
  coord_cartesian(ylim=c(0.4, 1)) + scale_y_continuous( breaks=seq(0.40, 1, 0.2))
p3 <- PlotRarefaction(sp.main.data$Microcebus_griseorufus$rarefaction$rs) + 
  theme(axis.text= element_text(size = 15)) +
  coord_cartesian(ylim=c(0, 1)) + scale_y_continuous("RS", breaks=seq(0.25, 1, 0.25))
p4 <- PlotRarefaction(sp.main.data$Indri_indri$rarefaction$rs) + 
  theme(axis.text= element_text(size = 15), axis.title.y = element_blank(), axis.text.y =element_blank(), axis.ticks.y = element_line(size =0)) +
  coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=seq(0.25, 1, 0.25))
plot_grid(p1, p2, p3, p4)  


