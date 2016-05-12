pc_plot14 <- ggplot(plot.W, aes(PC1, PC4)) +
  #geom_polygon(aes(hPC1, hPC4, fill = Family, color = Family, group= Especie, log = "x"), data = hulls, alpha=0.2) + 
  geom_point(data = points.info,
             aes (PC1, PC4, group= Genero, shape = Regiao, color = Family, log = "x"), 
             alpha = 0.7, size = 2) + 
  #theme(legend.text = element_text(angle = 90, size = 10, hjust = 3, vjust = 3, face = 'italics')) +
  #scale_shape(name = "Localização e status", labels = c("Madagascar - viventes", "Madagascar - Extintos", "Fora de Madagascar - Viventes")) +
  geom_text(data = ddply(points.info, .(Genero, Family), numcolwise(mean)),
            aes(PC1, PC4, label= Genero, color = Family),  alpha = 0.4 , size = 6) +
  #        scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
  #                labels = c("Small", "Medium", "Large") ) +
  scale_colour_brewer( palette="Set1") +
  theme_bw() + 
  #labs(x = "PC1 = Tamanho") +
  #labs(y = "PC4 = Contraste entre comprimento\nda face e volume do neurocrânio ") +
  #coord_cartesian(ylim=c(-2, 58), xlim = c(-10, 580)) +
  coord_fixed() +
  guides(size = FALSE, 
         text = FALSE, 
         colour = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ), 
         fill = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ),
         shape = guide_legend(title.position = "top", override.aes = list(alpha = 1, size = 5), nrow = 1) ) +
  theme(legend.position = "bottom", 
        axis.title= element_text(size=15),
        plot.title = element_text(lineheight=.8, size = 13))

pc_plot12 <- ggplot(plot.W, aes(PC1, PC2)) +
  #geom_polygon(aes(hPC1, hPC4, fill = Family, color = Family, group= Especie, log = "x"), data = hulls, alpha=0.2) + 
  geom_point(data = points.info,
             aes (PC1, PC2, group= Genero, shape = Regiao, color = Family, log = "x"), 
             alpha = 0.7, size = 2) + 
  #theme(legend.text = element_text(angle = 90, size = 10, hjust = 3, vjust = 3, face = 'italics')) +
  #scale_shape(name = "Localização e status", labels = c("Madagascar - viventes", "Madagascar - Extintos", "Fora de Madagascar - Viventes")) +
  geom_text(data = ddply(points.info, .(Genero, Family), numcolwise(mean)),
            aes(PC1, PC2, label= Genero, color = Family),  alpha = 0.4 , size = 6) +
  #        scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
  #                labels = c("Small", "Medium", "Large") ) +
  scale_colour_brewer( palette="Set1") +
  theme_bw() + 
  #labs(x = "PC1 = Tamanho") +
  #labs(y = "PC4 = Contraste entre comprimento\nda face e volume do neurocrânio ") +
  ggtitle("Species projection in ancestral Strepsirrhini W") +  #coord_cartesian(ylim=c(-2, 58), xlim = c(-10, 580)) +
  coord_fixed() +
  guides(size = FALSE, 
         label = FALSE, 
         colour = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ), 
         fill = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ),
         shape = guide_legend(title.position = "top", override.aes = list(alpha = 1, size = 5), nrow = 1) ) +
  theme(legend.position = "bottom", 
        axis.title= element_text(size=15),
        plot.title = element_text(lineheight=.8, size = 13))

pc_plot13 <- ggplot(plot.W, aes(PC1, -PC3)) +
  #geom_polygon(aes(hPC1, hPC4, fill = Family, color = Family, group= Especie, log = "x"), data = hulls, alpha=0.2) + 
  geom_point(data = points.info,
             aes (PC1, -PC3, group= Genero, shape = Regiao, color = Family, log = "x"), 
             alpha = 0.7, size = 2) + 
  #theme(legend.text = element_text(angle = 90, size = 10, hjust = 3, vjust = 3, face = 'italics')) +
  #scale_shape(name = "Localização e status", labels = c("Madagascar - viventes", "Madagascar - Extintos", "Fora de Madagascar - Viventes")) +
  geom_text(data = ddply(points.info, .(Genero, Family), numcolwise(mean)),
            aes(PC1, -PC3, label= Genero, color = Family),  alpha = 0.4 , size = 6) +
  #        scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
  #                labels = c("Small", "Medium", "Large") ) +
  scale_colour_brewer( palette="Set1") +
  theme_bw() + 
  #labs(x = "PC1 = Tamanho") +
  #labs(y = "PC4 = Contraste entre comprimento\nda face e volume do neurocrânio ") +
  #coord_cartesian(ylim=c(-2, 58), xlim = c(-10, 580)) +
  coord_fixed() +
  guides(size = FALSE, 
         text = FALSE, 
         colour = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ), 
         fill = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ),
         shape = guide_legend(title.position = "top", override.aes = list(alpha = 1, size = 5), nrow = 1) ) +
  theme(legend.position = "bottom", 
        axis.title= element_text(size=15),
        plot.title = element_text(lineheight=.8, size = 13))

pc_plot23 <- ggplot(plot.W, aes(PC2, -PC3)) +
  #geom_polygon(aes(hPC1, hPC4, fill = Family, color = Family, group= Especie, log = "x"), data = hulls, alpha=0.2) + 
  geom_point(data = points.info,
             aes (PC2, -PC3, group= Genero, shape = Regiao, color = Family, log = "x"), 
             alpha = 0.7, size = 2) + 
  #theme(legend.text = element_text(angle = 90, size = 10, hjust = 3, vjust = 3, face = 'italics')) +
  #scale_shape(name = "Localização e status", labels = c("Madagascar - viventes", "Madagascar - Extintos", "Fora de Madagascar - Viventes")) +
  geom_text(data = ddply(points.info, .(Genero, Family), numcolwise(mean)),
            aes(PC2, -PC3, label= Genero, color = Family),  alpha = 0.4, size = 6) +
  #        scale_size(name = "Size reference\nAverage PC1 score per genus", breaks = c(5.2, 5.6, 6.0),
  #                labels = c("Small", "Medium", "Large") ) +
  scale_colour_brewer( palette="Set1") +
  theme_bw() + 
  #labs(x = "PC1 = Tamanho") +
  #labs(y = "PC4 = Contraste entre comprimento\nda face e volume do neurocrânio ") +
  #coord_cartesian(ylim=c(-2, 58), xlim = c(-10, 580)) +
  coord_fixed() +
  guides(size = FALSE, 
         text = FALSE, 
         colour = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ), 
         fill = guide_legend(title.position = "top", title = "Família", size = 40, override.aes = list(alpha = 1) ),
         shape = guide_legend(title.position = "top", override.aes = list(alpha = 1, size = 5), nrow = 1) ) +
  theme(legend.position = "bottom", 
        axis.title= element_text(size=15),
        plot.title = element_text(lineheight=.8, size = 13))

pc_plot123 <- plot_grid(pc_plot12+ theme(legend.position = "none"), 
                    pc_plot13+ theme(legend.position = "none"), 
                    pc_plot23+ theme(legend.position = "none"), ncol =1, vjust = -2, rel_widths = c(1,3,1000))
  plot_grid(pc_plot123, pc_plot14, ncol = 1, rel_widths =  c(10,3))
