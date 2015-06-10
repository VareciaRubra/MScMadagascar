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
#Registrando o numero de cores : 3 em casa e 7 no lab
#para descobrir rodar no terminal: nproc
registerDoParallel(cores = 7)

arquivo.bruto = "Data/Hapalemur_Shower.csv"
arquivo.saida = "Data/Hapalemur_Clean"

#read csv and create table dataframe
raw.data <- tbl_df(read_csv(arquivo.bruto))
# forçando a ordem dos fatores como sendo a ordem que aparece 
# no proprio dataframe-------> 
#isso evita que funçoes da classe apply reoordenem os resultados numérica/alfabeticamente.
raw.data$Tombo <- factor (raw.data$Tombo, levels = unique(raw.data$Tombo) )

#substituindo os espaços nos momes de especies
raw.data$Especie %<>% gsub(" ", "_", .)
#removendo os ? nos casos de especies que tava em duvida 
#fazer depois de rodar o biplot PC1xPC2
raw.data$Especie %<>% gsub("[?]", "", .)
#Substituindo o . por _ nos labels 
names(raw.data) <- gsub(".", "_", names(raw.data), fixed = T)

#renomeando o BAOPI que escrevi errado lá no arquivo original de coletar dado
raw.data %<>% 
  rename(., BA_OPI = BAO_PI) 

#contando quantos tem por espécie e por Take (pra ter certeza de que se removeu um cara sairam as duas replicas)
raw.data %>% count(Especie) 
raw.data %>% count(Take)

#############################################
############## SEGUNDA PARTE ################
#############################################
#Removendo outliers (procedimento enquanto tou olhando os graficos, cada cara escroto que percebo add uma linha)
#Fazer isso junto com a rodada dos demais scripts
raw.data %<>% 
  filter(Tombo != "MO-1882-1559") %>% #medidas muito erradas
  filter(Tombo != "MO-1880-2510") %>% # BAOPI muito erradas
  filter(Tombo != "MO-1932-3354") %>% # BAOPI muito erradas
  filter(Tombo != "6291") %>% #BAOPI muito erradas
  filter(Tombo != "MO-1882-1557") %>% # BAOPI muito erradas
  filter(Tombo != "63357") %>% #LD_AS zuado
  filter(Tombo != "63355") %>% #APET_TS zuado
  filter(Tombo != "83664") #%>% #IS_PM zuado

#depois de limpar os zuados, salvar em um arquivo
#write.csv(raw.data, arquivo.saida, row.names= F)

#organizando databases por genero e outro por especie
raw.main.data <- dlply(raw.data, .(Especie), tbl_df)
all.raw.main.data<- dlply(raw.data, .(Genero), tbl_df)

raw.main.data$Hapalemur_griseus %>% count(Take)

