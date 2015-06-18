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
registerDoParallel(cores = 6)
#abrir no terminal htop para ver os cores trabalhando

arquivo.bruto = "Data/1Master_Organized_Factors.csv"
arquivo.saida =  "Data/1Master_Organized_Factors.csv"
raw.regular <- read.csv(arquivo.bruto, head = T)

#read csv and create table dataframe
raw.data <- tbl_df(read_csv(arquivo.bruto))
table(is.na(raw.data$Tombo))
raw.data<- tbl_df(raw.regular)
table(is.na(raw.data$Tombo))
# forçando a ordem dos fatores como sendo a ordem que aparece 
# no proprio dataframe-------> 
#isso evita que funçoes da classe apply reoordenem os resultados numérica/alfabeticamente.
raw.data$Tombo <- factor (raw.data$Tombo, levels = unique(raw.data$Tombo) )

#contando quantos tem por espécie e por Take (pra ter certeza de que se removeu um cara sairam as duas replicas)
raw.data %>% count(Especie) 
raw.data %>% count(Take)
#write.csv(raw.data, arquivo.saida , row.names= F)

#organizando databases por ordem taxonomica 
unique(raw.data$Regiao)
#raw.data$Regiao %<>% gsub("Tarsiidae", "Out_Madagascar", .)
All.raw.main.data<- dlply(raw.data, .(All), tbl_df)
raw.main.data<- dlply(raw.data, .(Genero), tbl_df)
raw.main.data<- dlply(raw.data, .(Especie), tbl_df)


write.csv(raw.data, arquivo.saida , row.names= F)

