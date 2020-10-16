library(vegan)


planilha <- read.table("Not Commit/planilha_completa.txt", header = T)
peixe <- read.table("Not Commit/Peixes.txt", header = T, stringsAsFactors = T)
#Constucting matrices------------------------------------------------------------------------------------------------
#Fish and Isolation variables####
peixe <- as.character(peixe$peixe)
traits <- read.csv("Not Commit/Traits.csv", header = T, stringsAsFactors = T)

peixe[which(peixe == "ausente")] <- "absent"
peixe[which(peixe == "presente")] <- "present"

peixe <- as.factor(peixe)

isolation <- c(rep("30",8),rep("120",8),rep("480",8))
isolation <- as.factor(isolation)
levels(isolation)
isolation <- relevel(isolation, ref="30")


peixe<- c(as.character(peixe),as.character(peixe),as.character(peixe))
isolation <- c(as.character(isolation),as.character(isolation),as.character(isolation))
peixe <- as.factor(peixe)
isolation <- as.factor(isolation)

SS <- as.factor(planilha$AM)

poca <- planilha$poca
ID <- as.factor(poca)




#Constructing Y matrix####
com <- planilha[,c(3:length(planilha))]

for(i in 1:dim(com)[2]){
  com[,i] <- as.numeric(as.character(com[,i]))
}

dim(com)

peixe<-peixe[which(is.na(com$B_Callibaetis) == F)]
isolation<-isolation[which(is.na(com$B_Callibaetis) == F)]
isolation <- relevel(isolation, ref="30")
SS<-SS[which(is.na(com$B_Callibaetis) == F)]
ID<-ID[which(is.na(com$B_Callibaetis) == F)]
com <- na.omit(com)




#com$Chironominae <- com$C_Apedilum + com$C_Asheum + com$C_Beardius + com$C_Caladomyia +
#com$C_Chironomus + com$C_Goeldchironomus + com$C_Parachironomus + com$C_Polypedilum + com$C_Tanytarsini

#com$Tanypodinae <- com$C_T_Ablabesmyia + com$C_T_Labrundinia + com$C_T_Larsia

#com<- com[,which(colnSSes(com)!= "C_Apedilum" &
#                     colnSSes(com)!= "C_Asheum"&
#                     colnSSes(com)!= "C_Beardius"&
#                     colnSSes(com)!= "C_Caladomyia"&
#                     colnSSes(com)!= "C_Chironomus"&
#                     colnSSes(com)!= "C_Goeldchironomus"&
#                     colnSSes(com)!= "C_Parachironomus"&
#                     colnSSes(com)!= "C_Polypedilum"&
#                     colnSSes(com)!= "C_Tanytarsini" &
#             colnSSes(com)!= "C_T_Ablabesmyia"&
#             colnSSes(com)!= "C_T_Labrundinia"&
#             colnSSes(com)!= "C_T_Larsia")]





####Removind some ponds for balanced design
ID
#C4
#C3
#B3
#A4

com_incomplete <- com[which(ID != "A4" & ID != "B3" & ID != "C3" & ID != "C4"),]
dim(com)
dim(com_incomplete)
com_incomplete_oc <- decostand(com_incomplete, method = "pa")
com_oc <- decostand(com, method = "pa")

com_incomplete <- com_incomplete[,which(colSums(com_incomplete_oc) > 0)]
com <- com[,which(colSums(com_oc) > 0)]

dim(com_incomplete)

isolation_incomplete <- isolation[which(ID != "A4" & ID != "B3" & ID != "C3" & ID != "C4")]
peixe_incomplete <- peixe[which(ID != "A4" & ID != "B3" & ID != "C3" & ID != "C4")]
SS_incomplete <- SS[which(ID != "A4" & ID != "B3" & ID != "C3" & ID != "C4")]
ID_incomplete <- ID[which(ID != "A4" & ID != "B3" & ID != "C3" & ID != "C4")]
ID_incomplete <- as.factor(as.character(ID_incomplete))


######

#####Constructing Y and environmental matrices for each sSSpling survey#####

isolation_SS1_Exp1 <- isolation[which(SS == "1")]
isolation_SS2_Exp1 <- isolation[which(SS == "2")]
isolation_SS3_Exp1 <- isolation[which(SS == "3")]
fish <- peixe
fish_SS1 <- peixe[which(SS == "1")]
fish_SS2 <- peixe[which(SS == "2")]
fish_SS3 <- peixe[which(SS == "3")]

com_SS1 <- com[which(SS == "1"), ]
com_SS1_oc <- decostand(com_SS1, method = "pa")
com_SS1_Exp1 <- com_SS1[,which(colSums(com_SS1_oc) > 0)]

com_SS2 <- com[which(SS == "2"), ]
com_SS2_oc <- decostand(com_SS2, method = "pa")
com_SS2_Exp1 <- com_SS2[,which(colSums(com_SS2_oc) > 0)]

com_SS3 <- com[which(SS == "3"), ]
com_SS3_oc <- decostand(com_SS3, method = "pa")
com_SS3_Exp1 <- com_SS3[,which(colSums(com_SS3_oc) > 0)]

com_SS2_SS3 <- com[which(SS == "2" | SS == "3"), ]
com_SS2_SS3_oc <- decostand(com_SS2_SS3, method = "pa")

SS2_SS3 <- SS[which(SS == "2" | SS == "3")]
com_SS2_SS3 <- com_SS2_SS3[,which(colSums(com_SS2_SS3_oc[which(SS2_SS3 == "2"),]) > 0 &
                                    colSums(com_SS2_SS3_oc[which(SS2_SS3 == "3"),]) > 0)]

Trait_SS2_SS3 <- traits[which(colSums(com_SS2_SS3_oc[which(SS2_SS3 == "2"),]) > 0 &
                                 colSums(com_SS2_SS3_oc[which(SS2_SS3 == "3"),]) > 0),]
Trait_SS2 <- traits[which(colSums(com_SS2_oc) > 0),]
Trait_SS3 <- traits[which(colSums(com_SS3_oc) > 0),]


set.seed(4)
first <- sample(row(com_SS3)[which(fish_SS3=="absent" & isolation_SS3_Exp1 == "30")], size = 1)
second <- sample(row(com_SS3)[which(fish_SS3=="absent" & isolation_SS3_Exp1 == "120")], size = 1)
com_SS3_Exp1_incomplete <- com_SS3_Exp1[-c(first,second),]
fish_incomplete_SS3 <- fish_SS3[-c(first,second)]
isolation_incomplete_SS3_Exp1 <- isolation_SS3_Exp1[-c(first,second)]


#####

isolation_SS2_SS3 <- isolation[which(SS == "2" | SS == "3")]
fish_SS2_SS3 <- peixe[which(SS == "2" | SS == "3")]
SS_SS2_SS3 <- SS[which(SS == "2" | SS == "3")]
ID_SS2_SS3 <- ID[which(SS == "2" | SS == "3")]


length(com_SS2_SS3_abundance)
length(isolation_SS2_SS3)
length(fish_SS2_SS3)
length(SS_SS2_SS3)

