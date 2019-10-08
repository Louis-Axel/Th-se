########################################
## Début du chargement des librairies ##
########################################


library(lubridate)
library(plyr)
library(nlme)
library(lme4)
library(lattice)
library(ggplot2)
library(sf)
library(predictmeans)
library(sme)
library(influence.ME)
library(corrplot)
library(rccdates)
library(Metrics)
library(agricolae)
library(ape)
library(predictmeans)
library(emmeans)
library(multcomp)
library(broom)
library(tidyverse)
library(ggpubr)
library(afex)
library(MLmetrics)
library(devtools)
library(ggbiplot)
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(klaR)
library(factoextra)
library(vegan)
library(usedist)
library(readr)
library(RColorBrewer)
library(ggthemes)
library(kml)
library(mclust)
library(fpc)
library(mgcv)
library(brms)

########################################
## Fin du chargement des librairies   ##
########################################

########################################
## Début du chargement des données    ##
########################################

setwd("C:/Users/Louis-Axel/Documents/Thèse/Données")
eausol = read.csv("sol_eau.csv", dec = ",", stringsAsFactors=FALSE)
compoflo= read.csv("composition_floristique.csv", dec = ",", stringsAsFactors=FALSE)
parcelle = read.csv("parcelle.csv")
setwd("C:/Users/Louis-Axel/Documents/")
sol_chimique <- read.csv2("sol_chimique_new.csv")

########################################
## Fin du chargement des données      ##
########################################

###############################
## Renommage des traitements ##
###############################

for (i in 1:length(parcelle[,1])){
  if(parcelle$Traitt[i] == 1){ parcelle$Traitt[i] = "0-Temoin"}
  if(parcelle$Traitt[i] == 2){ parcelle$Traitt[i] = "1- Min70"}
  if(parcelle$Traitt[i] == 3){ parcelle$Traitt[i] = "3-Lis40"}
  if(parcelle$Traitt[i] == 4){ parcelle$Traitt[i] = "8-Mix4070"}
  if(parcelle$Traitt[i] == 5){ parcelle$Traitt[i] = "4-Lis70"}
  if(parcelle$Traitt[i] == 6){ parcelle$Traitt[i] = "5-Com72"}
  if(parcelle$Traitt[i] == 7){ parcelle$Traitt[i] = "7-Mix72"}
  if(parcelle$Traitt[i] == 8){ parcelle$Traitt[i] = "6-Com12"}
  if(parcelle$Traitt[i] == 9){ parcelle$Traitt[i] = "9-Mix1248"}
  if(parcelle$Traitt[i] == 10){ parcelle$Traitt[i] = "2-Min120"}
}

########################################
## Modification nouvelle database     ##
########################################

sol_chimique$Code_site = sol_chimique$site
sol_chimique$Parcelle = sol_chimique$parcelle

############################################
## Création de la nouvelle table de C.org ##
############################################

sol_chimique$C.org1er =  as.numeric(as.character(sol_chimique$C.orga..g.100g.sec.))
sol_chimique$C.org2em = as.numeric(as.character(sol_chimique$C.orga..g.100g.MS105.)) *as.numeric(as.character(sol_chimique$MSr.105øC))
sol_chimique$C.org2em3 = as.numeric(as.character(sol_chimique$C.orga..g.100g.MS105.))

for (i in 1 : length(sol_chimique[,1])){
if ( is.na(sol_chimique$C.org1er[i])) {sol_chimique$C.org1er[i] = 0}
  if ( is.na(sol_chimique$C.org2em[i])) {sol_chimique$C.org2em[i] = 0}
  if ( is.na(sol_chimique$C.org2em3[i])) {sol_chimique$C.org2em3[i] = 0}
}

sol_chimique$C.orgNC = (sol_chimique$C.org1er + sol_chimique$C.org2em3)

# Dans C.orgNC on a les valeurs de carbone du sol sans modification 

sol_chimique$C.org = (sol_chimique$C.org1er + sol_chimique$C.org2em)

# Dans C.org on a les valeurs de carbone du sol avec modification

##################################################################################
##  Avec ce code on associe parclle et traitements dans la table de composition ##
##  floristique                                                                 ##
##################################################################################

for(j in 1:length(parcelle$Parcelle)){ 
  for(i in 1:length(compoflo$parcelle)){
    if (compoflo$parcelle[i] == parcelle$Parcelle[j] &
        compoflo$Code_site[i] == parcelle$Code_site[j]){
      compoflo$traitement[i] = parcelle$Traitt[j]
    }}}

#####################################################################################
## Avec ce code on associe parcelle et traitements dans la table de chimie du sol  ##
#####################################################################################

for(j in 1:length(parcelle$Parcelle)){ 
  for(i in 1:length(sol_chimique$Parcelle)){
    if (sol_chimique$Parcelle[i] == parcelle$Parcelle[j] &
        sol_chimique$Code_site[i] == parcelle$Code_site[j]){
      sol_chimique$traitement[i] = parcelle$Traitt[j]
    }}}

for(j in 1:length(parcelle$Parcelle)){ 
  for(i in 1:length(sol_chimique$Parcelle)){
    if (sol_chimique$Parcelle[i] == parcelle$Parcelle[j] &
        sol_chimique$Code_site[i] == parcelle$Code_site[j]){
      sol_chimique$bloc[i] = parcelle$Bloc[j]
    }}}

#########################
## Renommage des sites ## 
#########################

sol_chimique$Code_site = as.character(sol_chimique$Code_site)

for (i in 1:length(sol_chimique[,1])){
  if(sol_chimique$Code_site[i] == "LYC" | sol_chimique$Code_site[i] == "0-LYC" )
                                        { sol_chimique$Code_site[i] = "0-LYC"
                                          sol_chimique$SOL[i] = "arénosol"}
  if(sol_chimique$Code_site[i] == "SIC" | sol_chimique$Code_site[i] == "1-SIC" )
                                        { sol_chimique$Code_site[i] = "1-SIC"
                                          sol_chimique$SOL[i] = "andosol" }
  if(sol_chimique$Code_site[i] == "SED" | sol_chimique$Code_site[i] == "2-SED" )
                                        { sol_chimique$Code_site[i] = "2-SED" 
                                          sol_chimique$SOL[i] = "andosol"}
  if(sol_chimique$Code_site[i] == "MAR" | sol_chimique$Code_site[i] == "3-MAR" )
                                        { sol_chimique$Code_site[i] = "3-MAR" 
                                          sol_chimique$SOL[i] = "andosol"}
}

sol_chimique$Code_site = as.factor(sol_chimique$Code_site)

####################################################################################
## Dans la table chimie du sol on transforme la table pour avoir la Date correcte ##
####################################################################################

sol_chimique$Date_prelevement = dmy(sol_chimique$Date_prel)

#########################################
## On récupère une donnée année unique ##
#########################################

sol_chimique$Date2 = substr(as.character(sol_chimique$Date_prelevement), 1, 4)
sol_chimique$Date3 = as.Date(substr(as.character(sol_chimique$Date_prelevement), 1, 4), format = "%Y")
sol_chimique$date_test = as.numeric(sol_chimique$Date2) - 2004

#################################################################
## Création des tables initiales et finales (via les horizons) ##
#################################################################

Check_1 = F ; Check_2 = F ; Check_3 = F ; Check_4 = F 

for (i in 1 : length(sol_chimique[,1])){
  if (sol_chimique$Code_site[i] == levels(sol_chimique$Code_site)[1] & Check_1 == F){
    Horizon_1 = sol_chimique$horizon[i]
    Check_1 = TRUE
  }
  if (sol_chimique$Code_site[i] == levels(sol_chimique$Code_site)[2] & Check_2 == F){
    Horizon_2 = sol_chimique$horizon[i]
    Check_2 = TRUE
  }
  if (sol_chimique$Code_site[i] == levels(sol_chimique$Code_site)[3] & Check_3 == F){
    Horizon_3 = sol_chimique$horizon[i]
    Check_3 = TRUE
  }
  if (sol_chimique$Code_site[i] == levels(sol_chimique$Code_site)[4] & Check_4 == F){
    Horizon_4 = sol_chimique$horizon[i]
    Check_4 = TRUE
  }
}
Horizon = cbind(as.character(Horizon_1),
                as.character(Horizon_2),
                as.character(Horizon_3),
                as.character(Horizon_4))

# Nous venons de récupérer les horizons sur lesqueles ont été pratiqués les relevés annuels 
# pour les différents sites 

sol_chimique$horizon_num = as.numeric(sol_chimique$horizon)
sol_chimique$trai_num = as.numeric(substr(sol_chimique$traitement, 1, 1))
 
# Maintenans on stock la vieille base de donnée et on garde dans sol_chimique tous le reste

sol_chimique_old = sol_chimique

sol_chimique_new = sol_chimique

for (i in 1:length(sol_chimique[,1])){
  for (j in levels(sol_chimique$Code_site)){
  if (sol_chimique$Code_site[i] != j & Horizon[match(j, levels(sol_chimique$Code_site))] != as.character(sol_chimique$horizon[i])){
    sol_chimique_new[i,] = NA
  } 
  }
}

naa = list()

for (i in 1:length(sol_chimique_new[,1])){
if (is.na(sol_chimique_new$Code_site[i]))
  { 
  naa = c(naa, i)
}
}

sol_chimique_new = sol_chimique_new[-as.numeric(naa),]

#####
## On voit bien que la conversion ne peut pas fonctionner
#####

xyplot(C.org ~ Date3 | Code_site, data = sol_chimique_new[sol_chimique_new$traitement == "0-Temoin",], typ = c("p", "r"))
xyplot(C.orgNC ~ Date3 | Code_site, data = sol_chimique_new[sol_chimique_new$traitement == "0-Temoin",], typ = c("p", "r"))

sol_chimique = sol_chimique_new 

############################################
## Transformation des types de variables  ##
############################################

compoflo$Code_site = as.factor(compoflo$Code_site)
compoflo$parcelle = as.factor(compoflo$parcelle)
parcelle$Parcelle = as.factor(parcelle$Parcelle)
sol_chimique$traitement = as.factor(sol_chimique$traitement)
sol_chimique$CN = (sol_chimique$C.org) / (sol_chimique$N.tot)
sol_chimique_new = sol_chimique_new[-c(1011, 1015, 1020),]
sol_chimique_LYC = sol_chimique[sol_chimique$Code_site == "0-LYC",]
sol_chimique_SIC = sol_chimique[sol_chimique$Code_site == "1-SIC",]
sol_chimique_MAR = sol_chimique[sol_chimique$Code_site == "3-MAR",]
sol_chimique_SED = sol_chimique[sol_chimique$Code_site == "2-SED",]

sol_chimique_LYC_new = sol_chimique_new[sol_chimique_new$Code_site == "0-LYC",]
sol_chimique_SIC_new = sol_chimique_new[sol_chimique_new$Code_site == "1-SIC",]
sol_chimique_MAR_new = sol_chimique_new[sol_chimique_new$Code_site == "3-MAR",]
sol_chimique_SED_new = sol_chimique_new[sol_chimique_new$Code_site == "2-SED",]

######
## On va écrire sous forme de csv les différe
######

write.csv(sol_chimique, file = "DATAFRAME_sol_chimique.csv")
write.csv(sol_chimique_new, file = "DATAFRAME_sol_chimique_new.csv")

sol_chimique = read.csv("DATAFRAME_sol_chimique.csv", sep =",")
sol_chimique_new = read.csv("DATAFRAME_sol_chimique_new.csv", sep = ",")

##################
## Moyenne LYC  ##
##################

Moyenne_LYC = array(as.numeric(0), c(10,8,length(levels(as.factor(sol_chimique_LYC_new$Date2)))))

for (i in levels(as.factor(sol_chimique_LYC_new$Date2)))
{
  for (j in levels(as.factor(sol_chimique_LYC_new$traitement)))
  {
    for (l in levels(as.factor(sol_chimique_LYC_new$bloc)))
    {
      for (k in 1:length(sol_chimique_LYC_new$C.org))
      {
        if (sol_chimique_LYC_new$Date2[k] == i & as.factor(sol_chimique_LYC_new$traitement[k]) == j
            & as.factor(sol_chimique_LYC_new$bloc)[k] == l)
        {
          Moyenne_LYC[match(j, levels(as.factor(sol_chimique_LYC_new$traitement))),
                      match(l, levels(as.factor(sol_chimique_LYC_new$bloc))),
                      match(i, levels(as.factor(sol_chimique_LYC_new$Date2)))] = as.numeric(sol_chimique_LYC_new$C.orgNC[k])
          Moyenne_LYC[match(j, levels(as.factor(sol_chimique_LYC_new$traitement))), 
                      7,
                      match(i, levels(as.factor(sol_chimique_LYC_new$Date2)))] = as.numeric(sol_chimique_LYC_new$Date2)[k]
          Moyenne_LYC[match(j, levels(as.factor(sol_chimique_LYC_new$traitement))), 
                      8,
                      match(i, levels(as.factor(sol_chimique_LYC_new$Date2)))] = as.numeric(substr((sol_chimique_LYC_new$traitement)[k],1,1))
        
        }
      }
    }
  }
}

Moyenne_LYC

for (i in 1:10){
  for (j in 1:length(levels(as.factor(sol_chimique_LYC_new$Date2)))){
    Moyenne_LYC[i,4,j] = as.numeric(mean( as.numeric(as.character( Moyenne_LYC[i,1:3,j] ) ) , na.rm = T))
    Moyenne_LYC[i,5,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_LYC[i,1:3,j] ) ), na.rm = T))
    Moyenne_LYC[i,6,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_LYC[i,1:3,j] ) ), na.rm = T) / (sqrt(3)))
  }
}  


for (j in 1:length(data_frame_LYC[,1])){
  for (i in 0:9){
    if (data_frame_LYC$V8[j] == i)
      data_frame_LYC$initiale[j] = data_frame_LYC$V4[i+1]
    data_frame_LYC$err_initiale[j] = data_frame_LYC$V5[i+1]
  }
}

for(i in 1:length(data_frame_LYC[,1])){
  data_frame_LYC$err_initiale[j] = sd(data_frame_LYC$V4/data_frame_LYC$initiale)
  
  }

write.csv(data_frame_LYC, "data_frame_corgNC_LYC.csv")

ggplot(data_frame_LYC, aes(x = V7, y = V4, group = V8))+
  geom_smooth()+
  geom_errorbar( aes(x = V7, ymax = V4 + V5, ymin = V4 - V5))+
  scale_x_discrete(limits = seq(from = 2004, to = 2018, by = 2)) +
  geom_hline( aes(yintercept= initiale), color = "red", size=1)+
  geom_point()+facet_wrap(~ V8) + theme_bw()+ ggtitle("Evolution of carbon concentration for LYC") +
  xlab("Date") + ylab("Carbon")



ggplot(data_frame_LYC, aes(x = V7, y = (V4/initiale-1)*100))+
  scale_x_discrete(limits = seq(from = 2004, to = 2018, by = 2)) +
  geom_smooth()+
  geom_errorbar( aes(x = V7, ymax = ((V4/initiale + sqrt((V5/V4)^2+(err_initiale/initiale)^2))-1)*100,
                     ymin = ((V4/initiale - sqrt((V5/V4)^2+(err_initiale/initiale)^2))-1)*100))+
  geom_hline( aes(yintercept= 0), color = "red", size=1)+
  geom_point()+facet_wrap(~ V8) + theme_bw() + ggtitle("Evolution of carbon concentration in % for LYC") +
  xlab("Date") + ylab("Carbon / initiale")

#################
## Moyenne SIC ##
#################

Moyenne_SIC = array(as.integer(0), c(10,8,length(levels(as.factor(sol_chimique_SIC_new$Date2)))))

for (i in levels(as.factor(sol_chimique_SIC_new$Date2)))
{
  for (j in levels(as.factor(sol_chimique_SIC_new$traitement)))
  {
    for (l in levels(as.factor(sol_chimique_SIC_new$bloc)))
    {
      for (k in 1:length(sol_chimique_SIC_new$C.org))
      {
        if (sol_chimique_SIC_new$Date2[k] == i & as.factor(sol_chimique_SIC_new$traitement[k]) == j
            & as.factor(sol_chimique_SIC_new$bloc)[k] == l)
        {
          Moyenne_SIC[match(j, levels(as.factor(sol_chimique_SIC_new$traitement))),
                      match(l, levels(as.factor(sol_chimique_SIC_new$bloc))),
                      match(i, levels(as.factor(sol_chimique_SIC_new$Date2)))] = as.numeric(sol_chimique_SIC_new$C.orgNC[k])
          Moyenne_SIC[match(j, levels(as.factor(sol_chimique_SIC_new$traitement))), 
                      7,
                      match(i, levels(as.factor(sol_chimique_SIC_new$Date2)))] = as.numeric(sol_chimique_SIC_new$Date2)[k]
          Moyenne_SIC[match(j, levels(as.factor(sol_chimique_SIC_new$traitement))), 
                      8,
                      match(i, levels(as.factor(sol_chimique_SIC_new$Date2)))] = as.numeric(substr((sol_chimique_SIC_new$traitement)[k],1,1))
          
        }
      }
    }
  }
}

Moyenne_SIC

for (i in 1:10){
  for (j in 1:length(levels(as.factor(sol_chimique_SIC_new$Date2)))){
    Moyenne_SIC[i,4,j] = as.numeric(mean( as.numeric(as.character( Moyenne_SIC[i,1:3,j] ) ) , na.rm = T))
    Moyenne_SIC[i,5,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_SIC[i,1:3,j] ) ), na.rm = T))
    Moyenne_SIC[i,6,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_SIC[i,1:3,j] ) ), na.rm = T) / (sqrt(3)))
  }
}

data_frame_SIC = as.data.frame(rbind(Moyenne_SIC[,,1],Moyenne_SIC[,,2],Moyenne_SIC[,,3],Moyenne_SIC[,,4],Moyenne_SIC[,,5],Moyenne_SIC[,,6]
                                     ,Moyenne_SIC[,,7],Moyenne_SIC[,,8],Moyenne_SIC[,,9],Moyenne_SIC[,,9],Moyenne_SIC[,,10],
                                     Moyenne_SIC[,,11],Moyenne_SIC[,,12],Moyenne_SIC[,,13],Moyenne_SIC[,,14]))

for (j in 1:length(data_frame_SIC[,1])){
for (i in 0:9){
if (data_frame_SIC$V8[j] == i)
data_frame_SIC$initiale[j] = data_frame_SIC$V4[i+1]
data_frame_SIC$err_initiale[j] = data_frame_SIC$V5[i+1]
}
}

for(i in 1:length(data_frame_SIC[,1])){
data_frame_SIC$err_initiale[j] = sd(data_frame_SIC$V4/data_frame_SIC$initiale)
}

write.csv(data_frame_SIC, "data_frame_corgNC_SIC.csv")

ggplot(data_frame_SIC, aes(x = V7, y = V4 ,  group = V8))+
  geom_smooth()+
  geom_errorbar( aes(x = V7, ymax = V4 + V5, ymin = V4 - V5))+
  scale_x_discrete(limits = seq(from = 2004, to = 2018, by = 2)) +
  geom_hline( aes(yintercept= initiale), color = "red", size=1)+
  geom_point()+facet_wrap(~ V8, ncol = 5)+theme_bw() + ggtitle("Evolution of carbon concentration for SIC") +
  xlab("Date") + ylab("Carbon")

ggplot(data_frame_SIC, aes(x = V7, y = 100*(V4/initiale-1)))+
  scale_x_discrete(limits = seq(from = 2004, to = 2018, by = 2)) +
  geom_smooth() +
  geom_errorbar( aes(x = V7, ymax = 100*(V4/initiale + sqrt((V5/V4)^2+(err_initiale/initiale)^2)-1),
                     ymin = 100*(V4/initiale - sqrt((V5/V4)^2+(err_initiale/initiale)^2)-1)))+
  geom_point()+facet_wrap(~ V8, ncol = 5)+theme_bw()+ geom_hline(yintercept=1, 
                                                       color = "red", size=1) +
  ggtitle("Evolution of carbon concentration in % for SIC") +
  xlab("Date") + ylab("Carbon / initiale")

#####

####################
## Moyenne SED    ##
####################

Moyenne_SED = array(as.integer(0), c(10,8,length(levels(as.factor(sol_chimique_SED_new$Date2)))))

for (i in levels(as.factor(sol_chimique_SED_new$Date2)))
{
  for (j in levels(as.factor(sol_chimique_SED_new$traitement)))
  {
    for (l in levels(as.factor(sol_chimique_SED_new$bloc)))
    {
      for (k in 1:length(sol_chimique_SED_new$C.org))
      {
        if (sol_chimique_SED_new$Date2[k] == i & as.factor(sol_chimique_SED_new$traitement[k]) == j
            & as.factor(sol_chimique_SED_new$bloc)[k] == l)
        {
          Moyenne_SED[match(j, levels(as.factor(sol_chimique_SED_new$traitement))),
                      match(l, levels(as.factor(sol_chimique_SED_new$bloc))),
                      match(i, levels(as.factor(sol_chimique_SED_new$Date2)))] = as.numeric(sol_chimique_SED_new$C.orgNC[k])
          Moyenne_SED[match(j, levels(as.factor(sol_chimique_SED_new$traitement))), 
                      7,
                      match(i, levels(as.factor(sol_chimique_SED_new$Date2)))] = as.numeric(sol_chimique_SED_new$Date2)[k]
          Moyenne_SED[match(j, levels(as.factor(sol_chimique_SED_new$traitement))), 
                      8,
                      match(i, levels(as.factor(sol_chimique_SED_new$Date2)))] = as.numeric(substr((sol_chimique_SED_new$traitement)[k],1,1))
          
        }
      }
    }
  }
}

Moyenne_SED

for (i in 1:10){
  for (j in 1:length(levels(as.factor(sol_chimique_SED_new$Date2)))){
    Moyenne_SED[i,4,j] = as.numeric(mean( as.numeric(as.character( Moyenne_SED[i,1:3,j] ) ) , na.rm = T))
    Moyenne_SED[i,5,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_SED[i,1:3,j] ) ), na.rm = T))
    Moyenne_SED[i,6,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_SED[i,1:3,j] ) ), na.rm = T) / (sqrt(3)))
  }
}


data_frame_SED = as.data.frame(rbind(Moyenne_SED[,,1],Moyenne_SED[,,2],Moyenne_SED[,,3],Moyenne_SED[,,4],Moyenne_SED[,,5],Moyenne_SED[,,6]
                                     ,Moyenne_SED[,,7],Moyenne_SED[,,8],Moyenne_SED[,,9],Moyenne_SED[,,10]
                                     ,Moyenne_SED[,,11],Moyenne_SED[,,12],Moyenne_SED[,,13]))

for (j in 1:length(data_frame_SED[,1])){
  for (i in 0:9){
    if (data_frame_SED$V8[j] == i)
      data_frame_SED$initiale[j] = data_frame_SED$V4[i+1]
    data_frame_SED$err_initiale[j] = data_frame_SED$V5[i+1]
  }
}

for(i in 1:length(data_frame_SED[,1])){
  data_frame_SED$err_initiale[j] = sd(data_frame_SED$V4/data_frame_SED$initiale)
}

write.csv(data_frame_SED, "data_frame_corgNC_SED.csv")

ggplot(data_frame_SED, aes(x = V7, y = V4))+
  geom_smooth()+
  geom_errorbar( aes(x = V7, ymax = V4 + V5, ymin = V4 - V5))+
  scale_x_discrete(limits = seq(from = 2004, to = 2012, by = 2)) +
  geom_point()+facet_wrap(~ V8, ncol = 5)+theme_bw()+ geom_hline(aes(yintercept=initiale), 
                                                       color = "red", size=1)+
  ggtitle("Evolution of carbon concentration flat for SED") +
  xlab("Date") + ylab("Carbon")



ggplot(data_frame_SED, aes(x = V7, y = 100*(V4/initiale-1)))+
  scale_x_discrete(limits = seq(from = 2004, to = 2012, by = 2)) +
  geom_smooth() +
  geom_errorbar( aes(x = V7, ymax = 100*(V4/initiale + sqrt((V5/V4)^2+(err_initiale/initiale)^2)-1),
                     ymin = 100*(V4/initiale - sqrt((V5/V4)^2+(err_initiale/initiale)^2)-1)))+
  geom_point()+facet_wrap(~ V8, ncol = 5 )+theme_bw()+ geom_hline(yintercept=1, 
                                                       color = "red", size=1)+
  ggtitle("Evolution of carbon concentration in % for SED") +
  xlab("Date") + ylab("Carbon / initiale")


#####

#####################
## Moyenne MAR     ##
#####################

Moyenne_MAR = array(as.integer(0), c(10,8,length(levels(as.factor(sol_chimique_MAR_new$Date2)))))

for (i in levels(as.factor(sol_chimique_MAR_new$Date2)))
{
  for (j in levels(as.factor(sol_chimique_MAR_new$traitement)))
  {
    for (l in levels(as.factor(sol_chimique_MAR_new$bloc)))
    {
      for (k in 1:length(sol_chimique_MAR_new$C.org))
      {
        if (sol_chimique_MAR_new$Date2[k] == i & as.factor(sol_chimique_MAR_new$traitement[k]) == j
            & as.factor(sol_chimique_MAR_new$bloc)[k] == l)
        {
          Moyenne_MAR[match(j, levels(as.factor(sol_chimique_MAR_new$traitement))),
                      match(l, levels(as.factor(sol_chimique_MAR_new$bloc))),
                      match(i, levels(as.factor(sol_chimique_MAR_new$Date2)))] = as.numeric(sol_chimique_MAR_new$C.orgNC[k])
          Moyenne_MAR[match(j, levels(as.factor(sol_chimique_MAR_new$traitement))), 
                      7,
                      match(i, levels(as.factor(sol_chimique_MAR_new$Date2)))] = as.numeric(sol_chimique_MAR_new$Date2)[k]
          Moyenne_MAR[match(j, levels(as.factor(sol_chimique_MAR_new$traitement))), 
                      8,
                      match(i, levels(as.factor(sol_chimique_MAR_new$Date2)))] = as.numeric(substr((sol_chimique_MAR_new$traitement)[k],1,1))
          
        }
      }
    }
  }
}

Moyenne_MAR

for (i in 1:10){
  for (j in 1:length(levels(as.factor(sol_chimique_MAR_new$Date2)))){
    Moyenne_MAR[i,4,j] = as.numeric(mean( as.numeric(as.character( Moyenne_MAR[i,1:3,j] ) ) , na.rm = T))
    Moyenne_MAR[i,5,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_MAR[i,1:3,j] ) ), na.rm = T))
    Moyenne_MAR[i,6,j] = as.numeric(sd(  as.numeric(as.character( Moyenne_MAR[i,1:3,j] ) ), na.rm = T) / (sqrt(3)))
  }
}

data_frame_MAR = as.data.frame(rbind(Moyenne_MAR[,,1],Moyenne_MAR[,,2],Moyenne_MAR[,,3],Moyenne_MAR[,,4],Moyenne_MAR[,,5],Moyenne_MAR[,,6]
                                     ,Moyenne_MAR[,,7],Moyenne_MAR[,,8],Moyenne_MAR[,,9],Moyenne_MAR[,,10]))

for (j in 1:length(data_frame_MAR[,1])){
  for (i in 0:9){
    if (data_frame_MAR$V8[j] == i)
      data_frame_MAR$initiale[j] = data_frame_MAR$V4[i+1]
    data_frame_MAR$err_initiale[j] = data_frame_MAR$V5[i+1]
  }
}

for(i in 1:length(data_frame_MAR[,1])){
  data_frame_MAR$err_initiale[j] = sd(data_frame_MAR$V4/data_frame_MAR$initiale)
}

write.csv(data_frame_MAR, "data_frame_corgNC_MAR.csv")

ggplot(data_frame_MAR, aes(x = V7, y = V4))+
  geom_smooth()+
  geom_errorbar( aes(x = V7, ymax = V4 + V5, ymin = V4 - V5))+
  scale_x_discrete(limits = seq(from = 2004, to = 2012, by = 2)) +
  geom_point()+facet_wrap(~ V8, ncol = 5)+theme_bw()+geom_hline(aes(yintercept = initiale), color = 'red')



ggplot(data_frame_MAR, aes(x = V7, y = V4/initiale))+
  scale_x_discrete(limits = seq(from = 2004, to = 2012, by = 2)) +
  geom_smooth() +
  geom_errorbar( aes(x = V7, ymax = V4/initiale + sqrt((V5/V4)^2+(err_initiale/initiale)^2),
                     ymin = V4/initiale - sqrt((V5/V4)^2+(err_initiale/initiale)^2)))+
  geom_point()+facet_wrap(~ V8, ncol = 5)+geom_hline(aes(yintercept = 1), color = 'red')+theme_bw()

#####

#####################################################################################
##### Test sur les moyennes pour regarder quels traitements sont significatifs     ##
#####################################################################################

for (i in 1:10){
  a = t.test(Moyenne_LYC[1,1:3,], Moyenne_LYC[i,1:3,], paired = T)
  print(c( formatC(a$p.value, digits = 5, format = "e"), as.character(i)))
}


for (i in 1:10){
  a = t.test(Moyenne_SIC[1,1:3,], Moyenne_SIC[i,1:3,], paired = T)
  print(c( formatC(a$p.value, digits = 5, format = "e"), as.character(i)))
}

for (i in 1:10){
  a = t.test(Moyenne_SED[1,1:3,], Moyenne_SED[i,1:3,], paired = T)
  print(c( formatC(a$p.value, digits = 5, format = "e"), as.character(i)))
}

for (i in 1:10){
  a = t.test(Moyenne_MAR[1,1:3,], Moyenne_MAR[i,1:3,], paired = T)
  print(c( formatC(a$p.value, digits = 5, format = "e"), as.character(i)))
}

Moyenne_LYC

############################################
# Pour le moment les données en sont là et #
# sont chargeables directement             #
############################################


setwd("C:/Users/Louis-Axel/Documents/")

sol_chimique = read.csv("DATAFRAME_sol_chimique.csv", sep =",")
sol_chimique_new = read.csv("DATAFRAME_sol_chimique_new.csv", sep = ",")

sol_chimique_LYC = sol_chimique[sol_chimique$Code_site == "0-LYC",]
sol_chimique_SIC = sol_chimique[sol_chimique$Code_site == "1-SIC",]
sol_chimique_MAR = sol_chimique[sol_chimique$Code_site == "3-MAR",]
sol_chimique_SED = sol_chimique[sol_chimique$Code_site == "2-SED",]

sol_chimique_LYC_new = sol_chimique_new[sol_chimique_new$Code_site == "0-LYC",]
sol_chimique_SIC_new = sol_chimique_new[sol_chimique_new$Code_site == "1-SIC",]
sol_chimique_MAR_new = sol_chimique_new[sol_chimique_new$Code_site == "3-MAR",]
sol_chimique_SED_new = sol_chimique_new[sol_chimique_new$Code_site == "2-SED",]

data_frame_LYC = read.csv("data_frame_corgNC_LYC.csv")
data_frame_SIC = read.csv("data_frame_corgNC_SIC.csv")
data_frame_SED = read.csv("data_frame_corgNC_SED.csv")
data_frame_MAR = read.csv("data_frame_corgNC_MAR.csv")

#############################################
# Exemple type de réalisation de boxplot   ##
#############################################

xyplot(as.numeric(as.character(pH.eau)) ~ Date_prelevement | traitement, data = sol_chimique_LYC, typ = c("p", "r"))
xyplot(as.numeric(as.character(pH.eau)) ~ Date_prelevement | traitement, data = sol_chimique_SED, typ = c("p", "r"))
xyplot(as.numeric(as.character(pH.eau)) ~ Date_prelevement | traitement, data = sol_chimique_SIC, typ = c("p", "r"))
xyplot(as.numeric(as.character(pH.eau)) ~ Date_prelevement | traitement, data = sol_chimique_MAR, typ = c("p", "r"))

#####################################
## Ici nous commencons les modèles ##
#####################################

model_LYC = lme((sqrt(C.orgNC)) ~ traitement, random =~ Date3|traitement,data = sol_chimique_LYC_new,  method = "ML") 
summary(model_LYC)
plot(model_LYC)
hist(resid(model_LYC))
plot(model_LYC, Date3 ~ resid(.))
plot(model_LYC, sqrt((C.orgNC)) ~ fitted(.) | Date3, abline = c(0,1), id = 0.025, idResType = "pearson")
plot(model_LYC, sqrt((C.orgNC)) ~ fitted(.) | traitement, abline = c(0,1), id = 0.025, idResType = "pearson")
plot.lme(model_LYC, sqrt((C.orgNC)) ~ fitted(.) | Parcelle, abline = c(0,1), id = 0.025)
plot(model_LYC, (fitted(.)^2)~ as.numeric(Date3) | traitement)

## Note : Ici un simple modèle mixte sur la racine carré avec unb effet mixte sur le traitement en fonction de la data 

sol_chimique_SIC_new$traitement = as.factor(sol_chimique_SIC_new$traitement)
model_SIC = lme(sqrt(C.orgNC) ~ traitement , random =~ Date3|Parcelle, data = sol_chimique_SIC_new, method = "ML") 
r = glht(model_SIC, linfct = mcp(traitement = "Tukey"), test = adjusted("holm"))
tuk.cld <- cld(r, by = NULL, covar = F)
conf = confint(r)
old.par <- par(mai=c(1,1,1.5,1), no.readonly = F)
col = brewer.pal(n = 10, name = "Spectral")
plot(tuk.cld, col = col, type = "response")
rr = emmeans(model_SIC, list(pairwise ~ traitement), adjust = "tukey")
plot(rr)

## Ici on fait la même mais un test de Tukey supplémentaire est effectué

summary(model_SIC)
plot(model_SIC)
hist(resid(model_SIC))
plot(model_SIC, Date3 ~ fitted(.))
plot(model_SIC, sqrt((C.orgNC)) ~ fitted(.) | Date3, abline = c(0,1), id = 0.025, idResType = "pearson")
plot(model_SIC, sqrt((C.orgNC)) ~ fitted(.) | traitement, abline = c(0,1), id = 0.025, idResType = "pearson")
plot.lme(model_SIC, sqrt((C.orgNC)) ~ fitted(.) | Parcelle, abline = c(0,1), id = 0.025)
plot(model_SIC, (fitted(.)^2) ~ as.numeric(Date3) | traitement)

## Model total 
## Ici on va inclure les sites en tant qie facteurs dans notre analyse

sol_chimique_new$Date3 = as.factor(sol_chimique_new$Date3)
sol_chimique_new$Date2 = as.factor(sol_chimique_new$Date2)
sol_chimique_new$date_test = as.numeric(sol_chimique_new$date_test)

sol_chimique_SIC_new$Date2 = as.factor(sol_chimique_SIC_new$Date2)

model_Mtot = lme(log(C.orgNC)~ Date2 * traitement, 
                 random =~ 1|Parcelle, data = sol_chimique_SIC_new) 
AIC(model_Mtot) 

model_Mtot = lme(log(C.orgNC)~ traitement * Date2, 
                 random =~ 1|Code_site, data = sol_chimique_new) 
AIC(model_Mtot)

model_Mtot2 = lmer(log(C.orgNC) ~ 1 + traitement * (Code_site)  *
                     ((date_test)|Code_site),
                   data = sol_chimique_new)

AIC(model_Mtot2)

model_Mtot3 = lmer(log(C.orgNC) ~   traitement   +
                     (bs(date_test, 3) | Code_site),
                   data = sol_chimique_new)

model_Mtot4 = lmer(log(C.orgNC) ~     traitement * Code_site +
                     (bs(date_test, 3) | (Code_site)) *
                     (bs(date_test, 3) | (traitement)) ,
                   data = sol_chimique_new)
                   
## Utilisation du modèle mixte bayésien 
  
aa = brm(data = sol_chimique_new, log(C.orgNC) ~ traitement * Code_site +
             (date_test | traitement) +(date_test | traitement) , 
         core = 8, chain = 8, iter = 2000)

summary(aa)

newdata = data.frame(traitement = levels(sol_chimique_new$traitement),
                     Code_site = (sol_chimique_new$Code_site))
fit = fitted(aa,
             newdata = newdata,
             re_formula = NA,
              summary = T) 

colnames(fit) = c('fit', "se", "lwr", "upr")
df_plot = cbind(newdata, fit)
boxplot(df_plot$fit ~ df_plot$traitement)

fit1 = as.data.frame(fitted(aa,
             newdata = newdata,
             re_formula = NA,
             summary = F) )

colnames(fit1) = newdata$traitement

a_vs_b = fit1$`0-Temoin` - (
  fit1$`1- Min70`+
    fit1$`2-Min120`+
    fit1$`3-Lis40`+
    fit1$`4-Lis70`+
    fit1$`5-Com72`+
    fit1$`6-Com12`+
    fit1$`7-Mix72`+
    fit1$`8-Mix4070`+
    fit1$`9-Mix1248`
) 

hist(a_vs_b)
quantile(a_vs_b, probs = c(.5, .025, .975))

pp_check(aa)
marginal_effects(aa)
plot(aa, pars = c("traitement", "Code_site")) 
plot(marginal_effects(aa, effects = "Code_site:traitement"))

## Fin du modèle mixte bayésien

## Modèle mixtes avec spline cubique automatique : marche bien. Sans doute le même code que LePabic

gam.mod = gam(log(C.orgNC) ~  s(trai_num) + (Code_site),
                   data = sol_chimique_new)

summary(gam.mod)
anova(gam.mod)

plot(gam.mod$gam,pages=1)
sol_chimique_new$resid.gam.mod <- residuals(gam.mod, type = "pearson")
sol_chimique_new$fit.gam.mod <- residuals(gam.mod, type = "pearson")

ggplot(data = sol_chimique_new) +
  geom_point(aes(x = date_test, y = resid.gam.mod, group = Code_site, color = Code_site)) + 
  facet_wrap(~traitement)

gam.mod.1 <- gam(log(C.orgNC) ~ traitement * Code_site + s(date_test, by = Code_site) +
                   + s(date_test, by = traitement), data = sol_chimique_new)


sol_chimique_new$resid.gam.mod <- residuals(gam.mod.1, type = "pearson")
sol_chimique_new$fit.gam.mod <- predict(gam.mod.1, type = "link")
plot(sol_chimique_new$fit.gam.mod, sol_chimique_new$resid.gam.mod)

summary(gam.mod$lme)
summary(gam.mod$gam)
anova(gam.mod$gam) 
gam.check(gam.mod$gam) 

ggplot(data = sol_chimique_new) + geom_line(aes(x = date_test, y = fit.gam.mod, group = Code_site)) + 
  facet_wrap(~traitement)

plot(Variogram(gam.mod.1$lme, robust = TRUE, data = sol_chimique_new, form = ~date_test | 
                 traitement))


dat <- gamSim(sol_chimique_new$C.orgNC,n=200,scale=2)

## Fin du modèle mixte avec spline cubiques
  library(splines)


ab = anova(model_Mtot, model_Mtot2, model_Mtot3)
anova(model_Mtot3)
AIC(model_Mtot3)
summary(model_Mtot3)

anova.lme(model_Mtot2, model_Mtot3)
plot(model_Mtot3)
hist(resid(model_Mtot3))
plot(model_Mtot3, log(C.orgNC) ~ fitted(.) | Date3, abline = c(0,1), id = 0.0025, idResType = "pearson")
plot(model_Mtot3, (log(C.orgNC)) ~ fitted(.) | traitement, abline = c(0,1), id = 0.0025, idResType = "pearson")
plot(model_Mtot3, log(C.orgNC) ~ fitted(.) | Parcelle, abline = c(0,1), id = 0.025)


plot(model_Mtot3, exp(fitted(.)) ~ as.numeric(date_test) | traitement, 
     typ=c("r", "p"))


sol_chimique_new$fitted = fitted(gam.mod.1)

ggplot(sol_chimique_new, aes(x = Date2,
                             y = (fitted),
                             col = traitement,
                             fill = traitement,
                             group = traitement))+
  geom_smooth()+ facet_wrap(. ~ traitement, nrow = 2)

ggplot(sol_chimique_new, aes(x = Date2, y =  (exp(fitted)) , col = traitement,
                             fill = traitement, group = traitement))+
  geom_smooth()+ facet_wrap(. ~ Code_site, scales = "free")

a = summary(model_Mtot3)
sort(a$coefficients[,1])

xyplot(exp(fitted(model_Mtot3)) ~ as.numeric(sol_chimique_SIC_new$date_test) | 
         sol_chimique_SIC_new$traitement)


r = glht(model_Mtot3, linfct = mcp(traitement = "Tukey"), test = adjusted("holm"))

tuk.cld <- cld(r, by = NULL, covar = F)
plot(tuk.cld$x, tuk.cld$lp^3, col = col, ylim = c(0,30))
e = tuk.cld$mcletters$Letters
text(sort(unique(tuk.cld$x)), 25, labels = c(e[[1]],e[[2]],e[[3]],e[[4]],e[[5]],
                               e[[6]],e[[7]],e[[8]],e[[9]],e[[10]]))

rc = lsmeans(model_Mtot3, pairwise ~ traitement)

plot(rc)


mod <- lm(C.orgNC ~ traitement * Code_site + date_test, 
           data=sol_chimique_new)

sol_chimique_new = sol
fit <- sme(as.data.frame(sol_chimique_new[sol_chimique_new$traitement == "5-Com72",
                            c("C.orgNC", "traitement", "Code_site")]), lambda.mu=0,lambda.v=0)

sol_chimique_Corg$y = (sol_chimique_SIC_new$C.orgNC)
sol_chimique_Corg$z = (sol_chimique_SIC_new$date_test)
sol_chimique_Corg$code_site = (sol_chimique_SIC_new$Code_site)
sol_chimique_Corg$traitement = (sol_chimique_SIC_new$traitement)

sol_chimique_Corg = as.data.frame(sol_chimique_Corg)

sol_chimique_Corg$Zt <- smspline(~ z, data=sol_chimique_Corg)

fit1s <- lme(log(y) ~  traitement  , data=sol_chimique_Corg,random=list(z=pdIdent(~Zt - 1)))

plot(sol_chimique_Corg$z,sol_chimique_Corg$y,pch="o",type="n",main="Spline fits: lme(y ~ time, random=list(all=pdIdent(~Zt-1)))",xlab="time",ylab="y")
points(sol_chimique_Corg$z,sol_chimique_Corg$y,col=1)
points(sol_chimique_Corg$z, exp(fitted(fit1s)),col=2)
                                                                                                                                                                                                                                                
summary(fit1s)
smSplineEx1$Zt <- smspline(C.orgNC~ Date2, data=sol_chimique_new)
                           
summary(mod)
plot(mod)
anova(mod)
AIC(mod)
hist(resid(mod))
# r = glht(mod, linfct = mcp(traitement = "Tukey"), test = adjusted("holm"))
 tuk.cld <- cld(r, by = NULL, covar = F)
plot(tuk.cld$x, tuk.cld$lp, col = col, ylim = c(0,30)))
 e = tuk.cld$mcletters$Letters
 text(sort(unique(tuk.cld$x)), 25, labels = c(e[[1]],e[[2]],e[[3]],e[[4]],e[[5]],
                                             e[[6]],e[[7]],e[[8]],e[[9]],e[[10]]))

sol_chimique_new$fit = exp(fitted(mod))
xyplot(data = sol_chimique_new, fit - C.orgNC  ~ date_test | traitement, facet.by = traitement )
xyplot(data = sol_chimique_new, fit   ~ date_test | traitement, facet.by = traitement )
plot(fitted(mod))
plot(mod,(fitted(.)) - C.orgNC ~ Date2 | traitement)

## 

Moyenne_SIC

plot(predict(model_Mtot), residuals(model_Mtot))
a = Moyenne_SIC[,,1]

for (i in 2:length(Moyenne_SIC[1,1,])){
  b = (Moyenne_SIC[,,i])
  a= rbind(a,b)
}

a = as.data.frame(a)

mydata = as.data.frame(cbind(a$V7,
                             a$V8,
                             a$V4))

colnames(mydata) = c("trait", "année","CORG")

df.cluster <- reshape(mydata, timevar="année", idvar="trait", direction="wide")
names(df.cluster) <- c("id", paste("t", 1:(ncol(df.cluster)-1)))
df.cluster.cld <- clusterLongData(df.cluster)
kml(df.cluster.cld, nbClusters = 10, nbRedrawing=10000)
nouveau_choice= choice(df.cluster.cld)
mean.trajectories <- calculTrajMean(df.cluster.cld["traj"], df.cluster.cld['c2'][[1]]['clusters'])



sol_chimique_new$traitement = as.factor(sol_chimique_new$traitement)
sol_chimique_new$Code_site = as.factor(sol_chimique_new$Code_site)
model_all_tr = lme((sqrt(C.orgNC)) ~ traitement * Code_site  , random =~ Date2|Parcelle,
                   data = sol_chimique_new, method = "ML")

plot(model_all_tr)
lsm = lsmeans(model_all_tr, pairwise ~ traitement | Code_site)
plot(lsm)


tuk = glht(model_all_tr, mcp(traitement  = "Tukey"), test = adjusted("holm"))
tuk.cld <- cld(as.data.frame(tuk))
conf = confint(r)
old.par <- par(mai=c(1,1,1.5,1), no.readonly = F)
col = brewer.pal(n = 10, name = "Spectral")

plot(tuk.cld$x,  tuk.cld$lp)
kml(tuk.cld, nbClusters=2, nbRedrawing=1)

rr = cbind(tuk.cld$x, tuk.cld$y, tuk.cld$lp)

plot(rr[,2], rr[,3])

rr_LYC = rr[names(rr[,1]) == "0-LYC",]

boxplot(rr_LYC[,2] + rr_LYC[,3]  ~ rr_LYC[,1])

rr_SIC = rr[names(rr[,1]) == "1-SIC",]

boxplot(rr_SIC[,2] ~ rr_SIC[,1])

rr_SED = rr[names(rr[,1]) == "2-SED",]

boxplot(rr_SED[,2] ~ rr_SED[,1])

rr_MAR = rr[names(rr[,1]) == "3-MAR",]

boxplot(rr_MAR[,2] ~ rr_MAR[,1])


plot(tuk.cld, col = col, type = "response")





AIC(model_all_tr)
logLik(model_all_tr)
BIC(model_all_tr)

flm_all = lm(log(C.org)~ Code_site * traitement, data = sol_chimique)

AIC(flm_all)
logLik(flm_all)
BIC(flm_all)

par(mfrow= c(1,2))
plot( fitted(model_all_tr), log(sol_chimique$C.org))
plot( fitted(flm_all), log(sol_chimique$C.org))

par(mfrow= c(1,2))
plot( (fitted(model_all_tr) - log(sol_chimique$C.org))^2, ylim = c(0,2))
plot( (fitted(flm_all) - log(sol_chimique$C.org))^2, ylim = c(0,2))

rmse(fitted(model_all_tr) , log(sol_chimique$C.org))
rmse(fitted(flm_all) , log(sol_chimique$C.org))

sol_chimique$C.org.fitted = exp(fitted(model_all_tr))

plot(fitted(model_all_tr))

xyplot(C.org~ Date3 | traitement, data = sol_chimique[sol_chimique$Code_site == "3-MAR",], 
       fit = model_all_tr,reg.line = T,
       aspect = "fill", typ = c("p", "r"), strip = FALSE)


sol_chimique_ando = rbind(sol_chimique_MAR,sol_chimique_SIC, sol_chimique_SED)
for (i in 1:length(sol_chimique_ando[,1])) {
  if(sol_chimique_ando$C.org[i] < 7){sol_chimique_ando$C.org[i] = NA}}

sol_chimique$

spp = ggscatter(sol_chimique_ando, x = "date_test", y = "C.org", color = "traitement", 
               palette = col,  add = "reg.line",        # Color by groups "cyl"
               size = 1, alpha = 0.7, conf.int = TRUE   ,   facet.by = "traitement",                    # Change point shape by groups "cyl"
               ellipse = F, 
               fill = "traitement", cor.coef = F, ellipse.type = "norm"
) +
  stat_cor(method = "pearson", label.x = 2, label.y = 16) +
  stat_regline_equation(label.x = 2,label.y = 19)

spp

sp = ggscatter(sol_chimique_LYC, x = "date_test", y = "C.org", color = "traitement", 
               palette = col,  add = "reg.line",        # Color by groups "cyl"
               size = 1, alpha = 0.7, conf.int = TRUE   ,   facet.by = "traitement",                    # Change point shape by groups "cyl"
               ellipse = F, 
               fill = "traitement", cor.coef = F, ellipse.type = "norm"
) +
  stat_cor(method = "pearson", label.x = 2, label.y = ) +
  stat_regline_equation(label.x = 2,label.y = 7)

sp

sol_chimique

anova(model_all_tr, flm_all)

par(mfrow = c(1,1))

plot(fitted(flm_all), log(sol_chimique_cor$C.org))
points(fitted(model_all_tr), log(sol_chimique_cor$C.org), col = 'red')
abline(0,1)

par(mfrow = c(1,1))
q1 = qqnorm(fitted(flm_all) -log(sol_chimique_cor$C.org), plot.it = F)
q2 = qqnorm(fitted(model_all_tr) -log(sol_chimique_cor$C.org), col = 'red', plot.it = F)

plot(range(q1$x, q2$x), range(q1$y, q2$y), type = "n")
points(q1, cex = 0.2)
points(q2, col = "red", pch = 2, cex = 0.2)
qqline(fitted(model_all_tr) -log(sol_chimique_cor$C.org), col = 'red')
qqline(fitted(flm_all) -log(sol_chimique_cor$C.org))

points((q1$x-q2$x),(q1$y-q2$y))

m2 = mean(fitted(model_all_tr) -log(sol_chimique_cor$C.org))
sd2 = sd(fitted(model_all_tr) -log(sol_chimique_cor$C.org))
m1 = mean(fitted(flm_all) -log(sol_chimique_cor$C.org))
sd1 = sd(fitted(flm_all) -log(sol_chimique_cor$C.org))

hist(fitted(flm_all) -log(sol_chimique_cor$C.org), prob = T, ylim = c(0,4))
curve(dnorm(x, m1, sd1), add = T)
curve(dnorm(x, m2, sd2), col = 'red', add = T)
hist(fitted(model_all_tr) -log(sol_chimique_cor$C.org), prob = T, ylim = c(0,4))
curve(dnorm(x, m1, sd1), add = T)
curve(dnorm(x, m2, sd2), col = 'red', add = T)

m2 = mean(abs(fitted(model_all_tr) -log(sol_chimique_cor$C.org)))
sd2 = sd(abs(fitted(model_all_tr) -log(sol_chimique_cor$C.org)))
m1 = mean(abs(fitted(flm_all) -log(sol_chimique_cor$C.org)))
sd1 = sd(abs(fitted(flm_all) -log(sol_chimique_cor$C.org)))

hist(abs(fitted(flm_all) -log(sol_chimique_cor$C.org)), prob = T, ylim = c(0,10))
curve(dnorm(x, m1, sd1), add = T)
curve(dnorm(x, m2, sd2), col = 'red', add = T)
hist(abs(fitted(model_all_tr) -log(sol_chimique_cor$C.org)), prob = T, ylim = c(0,10))
curve(dnorm(x, m1, sd1), add = T)
curve(dnorm(x, m2, sd2), col = 'red', add = T)

d1 = (fitted(flm_all) -log(sol_chimique_cor$C.org))
d2 = (fitted(model_all_tr) -log(sol_chimique_cor$C.org))

m3 = mean((fitted(flm_all)-fitted(model_all_tr)))
sd3 = sd((fitted(flm_all)-fitted(model_all_tr)))
par(mfrow = c(1,1))
hist((fitted(flm_all)-fitted(model_all_tr)), prob = T, ylim = c(0,8))
curve(dnorm(x, m3, sd3), add = T)
                
qqnorm((fitted(flm_all)-fitted(model_all_tr)))
       
a = rmse(fitted(flm_all), log(sol_chimique_cor$C.org))
a1 = rmse(fitted(model_all_tr), log(sol_chimique_cor$C.org))
print( (1 - a1/a) * 100 )

## Le modèle mixte est bien meilleur que le modèle non mixte (on améliore de 12.33% )

model_all_tr = lme(log(C.org) ~ Code_site * traitement, random =~ Date3|Code_site,
                   data = sol_chimique_cor, method = "ML")

AIC(model_all_tr)
logLik(model_all_tr)

flm_all = glm(log(C.org)~ Code_site * traitement, data = sol_chimique)

AIC(flm_all)
logLik(flm_all)

## Le modèle mixte est bien meilleur que le modèle non mixte (on améliore de 19% avant correction )

summary(model_all_tr) ; anova(model_all_tr) ; intervals(model_all_tr) ; plot(model_all_tr)

shapiro.test(residuals(model_all_tr))
plot(predict(model_all_tr), residuals(model_all_tr))


model_all_tr = lme(log(C.org) ~ traitement * Code_site, random =~ Date3|Parcelle,
                   data = sol_chimique_cor, method = "ML")

plot(fitted(model_all_tr))

modlist <- lmer(log(C.org) ~ Code_site* traitement + (Date3|Code_site) ,
                data = sol_chimique_cor)


plot(fitted(modlist), log(sol_chimique_cor$C.org))
points(fitted(model_all_tr), log(sol_chimique_cor$C.org), col = 'red')

cooks2 <- CookD(model = model_all_tr)

plot((cooks2))

which.max(cooks2)

coodtr = CookD(model_all_tr, plot = F)
qqnorm(model_all_tr)
qqnorm(model_all_tr, ~ranef(.))

plot(model_all_tr, exp(fitted(.)) ~ (C.org) | Code_site, abline = c(0,1))

sol_chimique_PCA = sol_chimique[,c(4,5,6,7,9,10,11,12,13)]

res.pca = PCA(sol_chimique_PCA, scale.unit = TRUE, ncp = 5, graph = TRUE)

var <- get_pca_var(res.pca)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_ind(res.pca, col.ind="contrib") +
  scale_color_gradient2(low="blue", mid="red", midpoint=2, space ="Lab")

corrplot(var$cos2, is.corr=FALSE)

fviz_pca_ind(res.pca,
             geom.ind = c("point","text"), # Montre les points seulement (mais pas le "text"),
             col.ind = sol_chimique$traitement, # colorer by groups
              addEllipses = T, label = "var", invisible ="none",
             legend.title = "Groups",gradient.cols = redmono
)+
  ggpubr::fill_palette("simpsons")+      # Couleur des individus
  ggpubr::color_palette("simpsons")

sol_chimique$Date2 = substr(as.character(sol_chimique$Date_prelevement), 1, 4)

fviz_pca_biplot(res.pca, 
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = sol_chimique$traitement,
                addEllipses = T,
                legend.title = list(fill = "Traitement", color = "Traitement"),
                repel = TRUE        # Evite le chévauchement du texte
)+
   ggpubr::fill_palette("simpsons")+      # Couleur des individus
  ggpubr::color_palette("simpsons")   

sol_chimique_cor$ratio[i] = 0

## Il nous faut récupérer la veleur initiale pour chaque parcelle pour chaque traitement 
## pour chaque site

matrice_init = array(dim = c(30,4))

c(1:60, 91:150)

for (i in c(1:60, 91:150)){
  for (j in 1 : 30){
    for (k in levels(sol_chimique_cor$Code_site)){
      if (sol_chimique_cor$Parcelle[i] == j & sol_chimique_cor$Code_site[i] == k ){
        matrice_init[j, match(k, levels(sol_chimique_cor$Code_site))] = sol_chimique_cor$C.org[i]
      }
    }
  }
}

for (i in 1:length(sol_chimique$C.org)){ # Pour tous les relevés
  for (l in 1:30){
    for (k in levels(sol_chimique_cor$Code_site)){
    if (sol_chimique_cor$Parcelle[i] == l & sol_chimique_cor$Code_site[i] == k)
    {
      sol_chimique_cor$init[i] = matrice_init[l,match(k,levels(sol_chimique_cor$Code_site))]
    }
    }
  }}

sol_chimique_cor$diff = sol_chimique_cor$C.org - sol_chimique_cor$init

##

for (i in 1:length(sol_chimique_cor[,1])){

sol_chimique_cor$ratio[i] = (sol_chimique_cor$C.org[i] / sol_chimique_cor$init[i] ) }



sp = ggscatter(sol_chimique, x = "traitement", y = "C.org", color = "Code_site", 
               palette = "lancet",  add = "reg.line",        # Color by groups "cyl"
          size = 1, alpha = 0.7, conf.int = TRUE   ,   facet.by = "traitement",                    # Change point shape by groups "cyl"
           ellipse = F, 
          fill = "Code_site", cor.coef = T, ellipse.type = "norm"
)

sp
  stat_cor(aes(color = Code_site), label.x = 2)  

sp



sol_chimique$CN = sol_chimique$C.org / sol_chimique$N.tot
sp = ggscatter(sol_chimique, x = "pH", y = "N.tot",                       # Add regression line
               color = "traitement", palette = "aaas",  add = "reg.line",        # Color by groups "cyl"
               size = 1, alpha = 1, conf.int = F,                       # Change point shape by groups "cyl"
               fullrange = F, ellipse = T, repel = T, merge = T,
               fill = "lightgray"
)+
  stat_cor(aes(color = traitement), label.x = 2)  
xplot <- ggboxplot(sol_chimique_LYC, x = "traitement", y = "N.tot", 
                   color = "traitement", fill = "traitement", palette = "aaas",
                   alpha = 0.5, ggtheme = theme_bw())+
  rotate()

xplot
yplot <- ggboxplot(sol_chimique, x = "traitement", y = "N.tot",
                   color = "traitement", fill = "traitement", palette = c25,
                   alpha = 0.5, ggtheme = theme_bw())

sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")

plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "h", 
          rel_widths = c(1, 1), rel_heights = c(1, 7))



fm8 <- lmer(log(N.tot) ~  Code_site  * traitement  +
              (1  | Date_prelevement),
             sol_chimique, REML = T)

fm7= lme(log(N.tot) ~ Code_site * traitement,
          random =~ 1 | Date_prelevement, data = sol_chimique)


summary(fm7)
ggbiplot(fm7)
plot(PCA(fm7),alpha=0.05)

summary(fm8)
anova(fm8)

R = (random.effects(fm8))

cor(R)

plot(fm8)

coefficients(fm8)

rmse(fitted(fm8), sol_chimique$N.tot)
r.squaredGLMM(fm8)

summary(fm8)



cor = ranef(fm8)

MAPE(fitted(fm8), sol_chimique$N.tot)

rt = lmer(pH ~ traitement * Code_site +
       (Date_prelevement | Parcelle), data=sol_chimique)
r.squaredGLMM(fit1.lme)


xyplot(C.org  ~  factor(Date_prelevement) | traitement, 
       data = sol_chimique_SIC,grid = TRUE, type = c("p","smooth"), layout = c(1,1),
       index.cond = function(x, y) coef(lm(y ~ x))[1], auto.key =  T)
       
plot(rt)

confint(rt)

summary(rt)

resid(a)

a = sol_chimique_LYC$CEC - predict(fm8)
sol_chimique_LYC$diff = a
sol_chimique_LYC$pred = predict(fm8)

p <- ggplot(sol_chimique_SIC, aes(x= Date_prelevement, y=CEC, fill=traitement)) +
  geom_boxplot()+
 scale_size_manual(name="Predictions", values=c("Subjects"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 
p


q <- ggplot(sol_chimique, aes(x= Date_prelevement, y=pred, 
                                  fill=traitement)) +
  geom_line(aes( y = pred, group = traitement, 
                 col = traitement), lwd = 1.1)+
  geom_point(aes( y = pred, group = traitement))+
  scale_size_manual(name="Predictions", values=c("Subjects"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 
q

summary(fm8)

r <- ggplot(sol_chimique_LYC, aes(x= Date_prelevement, y=diff, 
                                  fill=traitement, group = traitement)) +
  geom_boxplot()+
  scale_size_manual(name="Predictions", values=c("Subjects"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 
r

anova(fm8)

distance = matrix(nrow = 10, ncol = 10)
factor = 4:8

for (i in levels(sol_chimique_LYC$traitement)){
  for (j in levels(sol_chimique_LYC$traitement)){
  a = mean(dist(rbind(sol_chimique_LYC[sol_chimique_LYC$traitement == i ,factor] , 
                    sol_chimique_LYC[sol_chimique_LYC$traitement == j ,factor] )))
  distance[match(i, levels(sol_chimique_LYC$traitement)),
           match(j, levels(sol_chimique_LYC$traitement))] = a 
  }
}

plot(hclust(as.dist(distance)))

for (i in levels(sol_chimique_SIC$traitement)){
  for (j in levels(sol_chimique_SIC$traitement)){
    a = mean(dist(rbind(sol_chimique_SIC[sol_chimique_SIC$traitement == i ,factor] , 
                   sol_chimique_SIC[sol_chimique_SIC$traitement == j ,factor] )))
    distance[match(i, levels(sol_chimique_SIC$traitement)),
             match(j, levels(sol_chimique_SIC$traitement))] = a 
  }
}

rownames(distance) = levels(sol_chimique$traitement)
colnames(distance) = levels(sol_chimique$traitement)
plot(hclust(as.dist(distance)))

for (i in levels(sol_chimique_SED$traitement)){
  for (j in levels(sol_chimique_SED$traitement)){
    a = dist(rbind(sol_chimique_SED[sol_chimique_SED$traitement == i ,4] , 
                   sol_chimique_SED[sol_chimique_SED$traitement == j ,4] ))
    distance[match(i, levels(sol_chimique_SED$traitement)),
             match(j, levels(sol_chimique_SED$traitement))] = a 
  }
}

plot(hclust(as.dist(distance)))

for (i in levels(sol_chimique_MAR$traitement)){
  for (j in levels(sol_chimique_MAR$traitement)){
    a = dist(rbind(sol_chimique_MAR[sol_chimique_MAR$traitement == i ,4] , 
                   sol_chimique_MAR[sol_chimique_MAR$traitement == j ,4] ))
    distance[match(i, levels(sol_chimique_MAR$traitement)),
             match(j, levels(sol_chimique_MAR$traitement))] = a 
  }
}

for (i in levels(sol_chimique$traitement)){
  for (j in levels(sol_chimique$traitement)){
    a = mean(dist(rbind(sol_chimique[sol_chimique$traitement == i ,factor] , 
                        sol_chimique[sol_chimique$traitement == j ,factor] )))
    distance[match(i, levels(sol_chimique$traitement)),
             match(j, levels(sol_chimique$traitement))] = a 
  }
}

hc = hclust(as.dist(distance))
plot(hc , hang = -1, cex = 0.6)
hcd <- as.dendrogram(hc)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
plot(hcd,  xlab = "Height", nodePar = nodePar, 
     edgePar = list(col = 2:3, lwd = 2:1))
plot(as.phylo(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)
plot(as.phylo(hc), type = "fan")
plot(as.phylo(hc), type = "cladogram", cex = 0.6,
     edge.color = "steelblue", edge.width = 2, edge.lty = 2,
     tip.color = "steelblue")

fviz_dend(hc, k = 4,                 # Cut in four groups
          cex = 0.8,                 # label size
          k_colors = "jco",
          color_labels_by_k = F,  # color labels by groups
          ggtheme = theme_bw(), rect = T     # Change theme
, horiz = T)

plot(hc)

hc %>% set("labels_col", value = c("green", "blue"), k=2) %>% 
  plot(main = "Color labels \nper cluster")

abline(h = 2, lty = 2)
for (i in 1:9){
  print(distance[i, (i+1):10])
}

distance2 = cbind(
            c(distance[2:10,1],
              distance[3:10,2],
              distance[4:10,3],
              distance[5:10,4],
              distance[6:10,5],
              distance[7:10,6],
              distance[8:10,7],
              distance[9:10,8],
              distance[10,9]),
            c("2","3","4","5","6","7","8","9","10",
              "3","4","5","6","7","8","9","10",
              "4","5","6","7","8","9","10",
              "5","6","7","8","9","10",
              "6","7","8","9","10",
              "7","8","9","10",
              "8","9","10",
              "9","10",
              "10"),
            c("1","1","1","1","1","1","1","1","1",
              "2","2","2","2","2","2","2","2",
              "3","3","3","3","3","3","3",
              "4","4","4","4","4","4",
              "5","5","5","5","5",
              "6","6","6","6",
              "7","7","7",
              "8","8",
              "9"))


c = list()

c$merge = matrix(c(distance2), nc = 10, byrow = T)
colnames(distance2) = c("Weight", "source", "target")

  a <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
a$merge <- matrix(c(-1, -2,
                    -3, -4,
                    1,  2), nc=2, byrow=TRUE ) 
a$height <- c(1, 1.5, 3)    # define merge heights
a$order <- 1:4              # order of leaves(trivial if hand-entered)
a$labels <- LETTERS[1:4]    # labels of leaves
class(a) <- "hclust"        # make it an hclust object
plot(a)                     # look at the result   

#convert to a dendrogram object if needed
ad <- as.dendrogram(a)

write.csv(file = "distance2.csv", distance2)

ME = dist(rbind(sol_chimique_MAR[,4] , sol_chimique_SED[,4]))
MI = dist(rbind(sol_chimique_MAR[,4] , sol_chimique_SIC[,4]))
ML = dist(rbind(sol_chimique_MAR[,4] , sol_chimique_LYC[,4]))
EI = dist(rbind(sol_chimique_SED[,4] , sol_chimique_SIC[,4]))
EL = dist(rbind(sol_chimique_SED[,4] , sol_chimique_LYC[,4]))
IL = dist(rbind(sol_chimique_SIC[,4] , sol_chimique_LYC[,4]))

link = c("ME", "MI", "ML", "EI", "EL", "IL")
weight = c(1/ME, 1/MI, 1/ML, 1/EI, 1/EL, 1/IL)
source = c("M", "M", "M", "E", "E", "I")
target = c("E", "I", "L", "I", "L", "L")

LINKS = as.data.frame(cbind( weight, source, target))

NODES = rbind("M", "I", "E", "L")
colnames(LINKS) = c("Weight", "SOURCE", "TARGET")

write.csv(file = "links.csv", LINKS)
write.csv(file = "nodes.csv", NODES)

sol_chimique_fin_SIC = sol_chimique[sol_chimique$Date2 == "2012" & 
                                      sol_chimique$Code_site == "1-SIC",]

compare_means(C.org ~ traitement, data = sol_chimique_fin_SIC)

p <- ggboxplot(sol_chimique_fin_SIC, x = "traitement", y = "C.org",
               color = "traitement", palette = "jco",
               add = "jitter")

# Change method
p + stat_compare_means()
p + stat_compare_means(method = "anova")

my_comparisons <- as.list( combn(levels(sol_chimique_fin_SIC$traitement), 1, simplify = F) )
# Pairwise comparison against reference
compare_means(C.org ~ traitement,  data = sol_chimique_fin_SIC, ref.group = "0-Temoin",
              method = "wilcox.test")
ggboxplot(sol_chimique_fin_SIC, x = "traitement", y = "C.org",
          color = "traitement", palette = "jco")+ 
  stat_compare_means(method = "anova", label.y = 18)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "0-Temoin")


ggboxplot(sol_chimique_fin, x = "traitement", y = "C.org", color = "traitement", 
          add = "jitter", legend = "none", facet.by = "Code_site", scales = "free_y") +
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova")   +   # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.", hide.ns = T)                      # Pairwise comparison against all

bwplot(C.org ~ Date2 | traitement,  data = sol_chimique_SED)
