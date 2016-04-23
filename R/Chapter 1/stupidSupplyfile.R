Modulemurs$Plot$Oral <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Oral") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") + 
  theme(axis.title.x = element_blank())

Modulemurs$Plot$Nasal <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Nasal") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Zygomatic <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Zygomatic") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Orbit <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Orbit") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Base <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Base") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Vault <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Vault") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Face <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Face") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Neuro <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Neuro") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
Modulemurs$Plot$Full <- Modulemurs$test.modularity.dist %>%  # pot em 10 x 15
  ldply(function(x) x$plotaisso) %>% filter (sp != "Euoticus_elegantulus") %>% 
  filter (sp != "Prolemur_simus") %>% 
  filter (hypotesis == "Full Integration") %>% 
  ggplot (data = ., aes (group = sp) ) + 
  geom_linerange(aes(x = sp, ymin = min, ymax = max, group = interaction (hypotesis, sp)), size =4, alpha = 0.2) +
  geom_point(aes(x = sp, y = mean, group = interaction (hypotesis, sp)), size = 4, color = "darkgrey", alpha = 0.7) +
  geom_point(aes(x = sp, y = abs(AVG.ratio), group = interaction (hypotesis, sp))) +
  facet_wrap(~ hypotesis) + theme_bw() + coord_flip() + xlab("") + ylab("AVG ratio") +
theme(axis.title.x = element_blank())
