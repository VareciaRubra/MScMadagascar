if(!require(plyr)) {install.packages('plyr'); library(plyr)}
if(!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
if(!require(magrittr)) {install.packages('magrittr'); library(magrittr)}
if(!require(lme4)) {install.packages('lme4'); library(lme4)}
if(!require(ggplot2)) {install.packages('ggplot2'); library(ggplot2)}
if(!require(tidyr)) {install.packages('tidyr'); library(tidyr)}
#if(!require(MCMCglmm)) {install.packages('MCMCglmm'); library(MCMCglmm)}
if(!require(reshape2)) {install.packages('reshape2'); library(reshape2)}
if(!require(evolqg)) {install.packages('evolqg'); library(evolqg)}
#if(!require(evolqg)) {devtools::install_github('lem-usp/evolqg'); library(evolqg)}
if(!require(readr)) {devtools::install_github('hadley/readr'); library(readr)}
if(!require(doParallel)) {install.packages('doParallel'); library(doParallel)}
#Registrando o numero de cores : 3 em casa, 7 no lab e até 10 no darwin
#para descobrir rodar no terminal: nproc
#abrir no terminal htop para ver os cores trabalhando
registerDoParallel(cores = 15)
#abrir no terminal htop para ver os cores trabalhando

arquivo.bruto = "Data/1Master_Organized_Factors.csv"
arquivo.saida =  "Data/1Master_Organized_Factors.csv"
#read csv and create table dataframe
raw.regular <- read.csv(arquivo.bruto, head = T)
raw.data<- tbl_df(raw.regular)
raw.data$Especie %<>% gsub("\\.", "", .)
# forçando a ordem dos fatores  como sendo a ordem que aparece 
# no proprio dataframe-------> 
#isso evita que funçoes da classe apply reoordenem os resultados numérica/alfabeticamente.

raw.data$Tombo <- factor (raw.data$Tombo, levels = unique(raw.data$Tombo) )

################## Árvore James ##################
treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
species <- treefile$tip.label[treefile$tip.label %in% unique(raw.data$Especie)]
pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)
raw.data %<>% 
  filter(Tombo != "044001b") %>% #comentario de pq foi removido
  filter(Tombo != "2386(16b)") %>%
  filter(Tombo != "16351") %>%
  filter(Tombo != "043004d") %>%
  filter(Tombo != "11882") #%>%
  

#contando quantos tem por espécie e por Take (pra ter certeza de que se removeu um cara sairam as duas replicas)
raw.data %>% count(Especie) 
raw.data %>% count(Take)
#write.csv(raw.data, arquivo.saida , row.names= F)

#organizando databases por ordem taxonomica 
unique(raw.data$Regiao)
#raw.data$Regiao %<>% gsub("Tarsiidae", "Out_Madagascar", .)
All.raw.main.data<- dlply(raw.data, .(All), tbl_df)
Gen.raw.main.data<- dlply(raw.data, .(Genero), tbl_df)
Sp.raw.main.data<- dlply(raw.data, .(Especie), tbl_df)


#write.csv(raw.data, arquivo.saida , row.names= F)

