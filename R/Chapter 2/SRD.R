#######################################################################################################################################
########################################################################################################################################
# SRD results exploration

mx.list.taxonomy <- cov.mx[mask]

mx.list.taxonomy$Tarsius <- ancestral.mx$'43'
mx.list.taxonomy$Microcebus <- ancestral.mx$'51'
mx.list.taxonomy$Mirza <- ancestral.mx$Mirza_coquereli
mx.list.taxonomy$Cheirogaleus <- ancestral.mx$'52'
mx.list.taxonomy$Phaner <- ancestral.mx$Phaner_furcifer
mx.list.taxonomy$W.Cheirogaleidae <- ancestral.mx$'48'
mx.list.taxonomy$Lepilemur <- ancestral.mx$'53'
mx.list.taxonomy$Avahi <- ancestral.mx$'58'
mx.list.taxonomy$Propithecus <- ancestral.mx$'59'
mx.list.taxonomy$Indri <- ancestral.mx$Indri_indri
mx.list.taxonomy$W.Indridae <- ancestral.mx$'56'
mx.list.taxonomy$Eulemur <- ancestral.mx$'65'
mx.list.taxonomy$Hapalemur <- ancestral.mx$Hapalemur_griseus
mx.list.taxonomy$Prolemur <- ancestral.mx$Prolemur_simus
mx.list.taxonomy$Lemur <- ancestral.mx$Lemur_catta
mx.list.taxonomy$Varecia <- ancestral.mx$'76'
mx.list.taxonomy$W.Lemuridae <- ancestral.mx$'63'
mx.list.taxonomy$Daubentonia <- ancestral.mx$Daubentonia_madagascariensis
mx.list.taxonomy$W.Madagascar <- ancestral.mx$'45'
mx.list.taxonomy$Perodicticus <- ancestral.mx$Perodicticus_potto
mx.list.taxonomy$Loris <- ancestral.mx$Loris_tardigradus
mx.list.taxonomy$W.Lorisidae <- ancestral.mx$'78'
mx.list.taxonomy$Nycticebus <- ancestral.mx$Nycticebus_coucang
mx.list.taxonomy$Euoticus <- ancestral.mx$Euoticus_elegantulus
mx.list.taxonomy$Otolemur <- ancestral.mx$Otolemur_crassicaudatus
mx.list.taxonomy$Galago <- ancestral.mx$Galago_senegalensis
mx.list.taxonomy$W.Galagidae <- ancestral.mx$'80'
mx.list.taxonomy$W.OutMadagascar <- ancestral.mx$'77'
mx.list.taxonomy$W.Strepsirrhini <- ancestral.mx$'44'
mx.list.taxonomy$W.Prosimian <- ancestral.mx$'42'
mx.list.taxonomy$Saguinus.P <- Saguinus_P.cov
mx.list.taxonomy$Saguinus.G <- Saguinus_G.cov



SRD(cov.x = mx.list.taxonomy)


srd.results.all <- SRD(cov.x = cov.mx[mask], parallel = TRUE)
