# format Rev treefile 

setwd("~/Desktop/relaxed_dating/data/")
library(ape)

tree <- read.tree("Bomarea_multilocus2_BS10_best.tre")

# drop tips from intraspecific sampling 
to_drop <- c(
  # b. edulis clade 1: everything but salta
  #only keep morelos
  "Bomarea_edulis_Veracruz_Tribble65",
  "Bomarea_edulis_Oaxaca_Tribble76",
  #"Bomarea_edulis_MorelosN_Tribble325",
  "Bomarea_acutifolia_Mexico_Tribble79",
  "Bomarea_edulis_Guatemala_Breckon2075",
  "Bomarea_edulis_Chiapas_Tribble90",
  
  "Bomarea_edulis_Argentina_Deginani1587",
  "Bomarea_edulis_Brazil_CaxiasdoSul_Kegler89",
  "Bomarea_edulis_Argentina_M_ulguradeRomero1876",
  "Bomarea_dolichocarpa_Peru_Barbour5069",
  "Bomarea_edulis_Brazil_Bahia_Thomas13631",
  "Bomarea_edulis_Brazil_MatoGrosso_Eiten9877",
  
  "Bomarea_edulis_Bolivia_Solomon7616",
  
  "Bomarea_edulis_MorelosS_Tribble377",
  "Bomarea_edulis_Veracruz_Tribble552",
  "Bomarea_edulis_Queretaro_Tribble701",
  "Bomarea_edulis_CR_J_F_Morales6635",
  "Bomarea_edulis_Nicaragua_Stevens21839",
  "Bomarea_edulis_Panama_Duke13712",
  
  "Bomarea_edulis_Suriname_Lindeman455",
  "Bomarea_edulis_Venezuela_AymardC_1314",

  "Bomarea_edulis_Cuba_Smith3264",
  "Bomarea_edulis_Haiti_Proctor10821",
  "Bomarea_edulis_DR_Jimenez3903",
 
  "Bomarea_edulis_Nayarit_Tribble113",
  "Bomarea_edulis_Honduras_Harmon3784",
  
  # B edulis clade 2 Argentina (salta)
  # keep in (1 tip)
  "Bomarea_edulis_Argentina_Hunziker12327",
  #"Bomarea_edulis_Argentina_Novara2378"
  
  # James Graham samples - keep 1
  #"Bomarea_sp__enanorojo_Peru_Graham12607",
  "Bomarea_sp__enanoverde_Peru_Graham12600",
  "Bomarea_sp__enanoverde_Peru_Graham12609",
  "Bomarea_sp__enanorojo_Peru_Graham12599",
  
  # synonyms, keep 1
  "Bomarea_lehmannii_AlzateS_N",
  #"Bomarea_straminea_Alzate3300",
  #"Bomarea_angustipetala_Alzate5116",    
  
  # "Bomarea_multiflora_CultivatedinCAfromCol_Greenhouse",
  "Bomarea_multiflora_AlzateS_N",
  
  "Bomarea_obovata_CR_Bonifacino6050"
  #"Bomarea_obovata_Ecuador_Clark4985",
)

tree <- treeio::drop.tip(tree, tip = to_drop)

# date tree with PL 

mycalibration <- data.frame(
  node = c( # Luzuriaga stem age fossil calibration
           getMRCA(tree, tip = c("Luzuriaga_polyphylla_CultivatedinCA_UCBG90_2401",
                                 "Drymophila_moorei_Australia_Copeland4560")),
            # Bomarea crown age secondary calibration 
           getMRCA(tree, tip = c("Bomarea_salsilla_Chile_Ackerman545",
                                 "Bomarea_multiflora_CultivatedinCAfromCol_Greenhouse")),
            # Alstroemeria crown age secondary calibration 
           getMRCA(tree, tip = c("Alstroemeria_apertiflora_Hatschbach17552", 
                                 "Alstroemeria_revoluta_Watson6608")),
            # root age secondary calibration 
           getMRCA(tree, tip = c("Luzuriaga_polyphylla_CultivatedinCA_UCBG90_2401", 
                                 "Bomarea_salsilla_Chile_Ackerman545"))),
  age.min = c(23.2, 14.2, 18.3, 50),
  age.max = c(33.2, 14.4, 18.5, 55),
  soft.bounds = c(FALSE, FALSE, FALSE, FALSE))


##### Use penalized likelihood to date the tree
tree_dated <- chronos(tree, lambda = 1, model = "relaxed",
                      calibration = mycalibration,
                      control = chronos.control() )

write.tree(tree_dated, file = "rev_starting.tre")
