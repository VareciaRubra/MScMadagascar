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
registerDoParallel(cores = 30)
#abrir no terminal htop para ver os cores trabalhando

arquivo.bruto = "Data/2Master_Organized_Factors.csv"
arquivo.saida =  "Data/3MasterDB_following_JHerrera_PhylOrder.csv"
#read csv and create table dataframe
raw.regular <- read.csv(arquivo.bruto, head = T)
original.names<- unique(raw.regular$Especie)
raw.regular$Especie %<>% gsub("\\.", "", .)

########## Renomeando minhas especies de acordo com a arvore do James ############
recoderFunc <- function(data, oldvalue, newvalue) {
  # convert any factors to characters
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  data[data == 0] <- NA
  newvec <- data
  # put recoded values into the correct position in the return vector
  for (i in unique(oldvalue)) newvec[data == i ] <- newvalue[oldvalue == i ]
  newvec
}

rename.james <- read.csv("myspecies.csv", header = TRUE)
raw.regular$Especie<- recoderFunc(data = raw.regular$Especie, oldvalue = rename.james$Anna, newvalue = rename.james$James.Edited)

################################ Colocando a ordem dos fatores de espécie como sendo a ordem em que eles aparecem na árvore
################## Árvore James ##################
treefile = read.nexus(file = "~/ataches/fbd369agerange_gooddates.tre")
species <- treefile$tip.label[treefile$tip.label %in% unique(raw.regular$Especie)]
################# Arvore indexada pelas especies que eu tenho na amostra #############
pruned.tree<-drop.tip(treefile,treefile$tip.label[-match(species, treefile$tip.label)])
plot(pruned.tree)
tiporder<- treefile$edge[,2][ treefile$edge[, 2] %in% 1:length(treefile$tip.label) ]
ordered.species <- as.data.frame(as.character(treefile$tip.label[tiporder]) )

ordered.species<- read.csv("ordered_species_james_plus.csv", header = TRUE)
james.original.order<- as.data.frame(!is.na(ordered.species$.id))
james.edited.order <- as.data.frame(ordered.species$.editedid)
current.order <- james.edited.order
names(current.order)<- ".id"

my.species <- unique(raw.regular$Especie)
species <- current.order[current.order$.id %in% my.species,]
raw.regular$Especie <- factor (raw.regular$Especie, levels = species )
######### salvando no formato tbl_df
raw.data<- tbl_df(raw.regular)
raw.data %<>% arrange(. , Especie)
# forçando a ordem dos fatores  como sendo a ordem que aparece 
# no proprio dataframe------->  ordenado de acordo com a arvore do James ################
#isso evita que funçoes da classe apply reoordenem os resultados numérica/alfabeticamente.
raw.data$Especie <- factor (raw.data$Especie, levels = unique(raw.data$Especie) )
raw.data$Genero <- factor (raw.data$Genero, levels = unique(raw.data$Genero) )
raw.data$Familia <- factor (raw.data$Familia, levels = unique(raw.data$Familia) )
raw.data$Tombo <- factor (raw.data$Tombo, levels = unique(raw.data$Tombo) )

raw.data %<>% 
  filter(Tombo != "044001b") %>% #comentario de pq foi removido
  filter(Tombo != "2386(16b)") %>%
  filter(Tombo != "16351") %>%
  filter(Tombo != "043004d") %>%
  filter(Tombo != "11882") #%>%
  

#contando quantos tem por espécie e por Take (pra ter certeza de que se removeu um cara sairam as duas replicas)
raw.data %>% dplyr::count(Especie) 
raw.data %>% dplyr::count(Take)
#write.csv(raw.data, arquivo.saida , row.names= F)

#organizando databases por ordem taxonomica 
unique(raw.data$Regiao)
#raw.data$Regiao %<>% gsub("Tarsiidae", "Out_Madagascar", .)
All.raw.main.data<- dlply(raw.data, .(All), tbl_df)
Gen.raw.main.data<- dlply(raw.data, .(Genero), tbl_df)
Sp.raw.main.data<- dlply(raw.data, .(Especie), tbl_df)


write.csv(raw.data, arquivo.saida.2 , row.names= F)

