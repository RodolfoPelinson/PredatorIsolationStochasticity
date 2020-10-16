library(lmPerm)
library(iNEXT)
library(lme4)
library(psycho)
library(tidyverse)
library(car)
library(RVAideMemoire)
source("Loading Data Exp 1.R")
source("beta_deviation.R")

isolation_SS1_Exp1 <- as.character(isolation_SS1_Exp1)
isolation_SS1_Exp1[which(isolation_SS1_Exp1 == "30")] <- "030"

isolation_SS2_Exp1 <- as.character(isolation_SS2_Exp1)
isolation_SS2_Exp1[which(isolation_SS2_Exp1 == "30")] <- "030"

isolation_SS3_Exp1 <- as.character(isolation_SS3_Exp1)
isolation_SS3_Exp1[which(isolation_SS3_Exp1 == "30")] <- "030"

isolation_Exp1 <- as.factor(c(isolation_SS1_Exp1, isolation_SS2_Exp1, isolation_SS3_Exp1))
fish <- as.factor(c(as.character(fish_SS1), as.character(fish_SS2), as.character(fish_SS3)))

isolation_Exp1 <- as.factor(c(isolation_SS2_Exp1, isolation_SS3_Exp1))
fish <- as.factor(c(as.character(fish_SS2), as.character(fish_SS3)))

fish_isolation_SS1 <- rep(NA, length(fish_SS1))
for(i in 1:length(fish_SS1)){
  fish_isolation_SS1[i] <- paste(isolation_SS1_Exp1[i],fish_SS1[i])
}

fish_isolation_SS2 <- rep(NA, length(fish_SS2))
for(i in 1:length(fish_SS2)){
  fish_isolation_SS2[i] <- paste(isolation_SS2_Exp1[i],fish_SS2[i])
}

fish_isolation_SS3 <- rep(NA, length(fish_SS3))
for(i in 1:length(fish_SS3)){
  fish_isolation_SS3[i] <- paste(isolation_SS3_Exp1[i],fish_SS3[i])
}

#fish_isolation <- as.factor(as.character(c(fish_isolation_SS1,fish_isolation_SS2,fish_isolation_SS3)))

fish_isolation <- as.factor(as.character(c(fish_isolation_SS2,fish_isolation_SS3)))

fish_isolation_SS1 <- rep(NA, length(fish_SS1))
for(i in 1:length(fish_SS1)){
    fish_isolation_SS1[i] <- paste("1",isolation_SS1_Exp1[i],fish_SS1[i])
}

fish_isolation_SS2 <- rep(NA, length(fish_SS2))
for(i in 1:length(fish_SS2)){
  fish_isolation_SS2[i] <- paste("2",isolation_SS2_Exp1[i],fish_SS2[i])
}

fish_isolation_SS3 <- rep(NA, length(fish_SS3))
for(i in 1:length(fish_SS3)){
  fish_isolation_SS3[i] <- paste("3",isolation_SS3_Exp1[i],fish_SS3[i])
}

AM_2_3 <- as.factor(as.character(AM[which(AM == "2" | AM == "3")]))
ID_2_3 <- as.factor(as.character(ID[which(AM == "2" | AM == "3")]))

#All <- c(as.character(fish_isolation_SS1),as.character(fish_isolation_SS2),as.character(fish_isolation_SS3))
All <- c(as.character(fish_isolation_SS2),as.character(fish_isolation_SS3))

All <- as.factor(All)

fish_isolation <- factor(fish_isolation, levels = c("030 absent","120 absent","480 absent","030 present","120 present","480 present"))
All <- factor(All, levels = c("2 030 absent","2 120 absent","2 480 absent","2 030 present","2 120 present","2 480 present", "3 030 absent","3 120 absent","3 480 absent","3 030 present","3 120 present","3 480 present"))



######################## droping samples where fish died ####
com_SS2_SS3 <- com_SS2_SS3[-which(ID_2_3 == "A2" #),]
                                 &AM_2_3 == "3"),]

fish_isolation <- fish_isolation[-which(ID_2_3 == "A2" #)]
                                       &AM_2_3 == "3")]

fish <- fish[-which(ID_2_3 == "A2" #)]
                   &AM_2_3 == "3")]

All <- All[-which(ID_2_3 == "A2" #)]
                 &AM_2_3 == "3")]

isolation_Exp1 <- isolation_Exp1[-which(ID_2_3 == "A2"#)]
                                       &AM_2_3 == "3")]


ID_2_3_2 <- ID_2_3[-which(ID_2_3 == "A2"#)]
                       &AM_2_3 == "3")]

AM_2_3 <- AM_2_3[-which(ID_2_3 == "A2"#)]
                        & AM_2_3 == "3")]

ID_2_3 <- ID_2_3_2







com_SS2_SS3 <- com_SS2_SS3[-which(ID_2_3 == "A7" #),]
                                  &AM_2_3 == "3"),]

fish_isolation <- fish_isolation[-which(ID_2_3 == "A7" #)]
                                        &AM_2_3 == "3")]

fish <- fish[-which(ID_2_3 == "A7" #)]
                    &AM_2_3 == "3")]

All <- All[-which(ID_2_3 == "A7" #)]
                  &AM_2_3 == "3")]

isolation_Exp1 <- isolation_Exp1[-which(ID_2_3 == "A7"#)]
                                        &AM_2_3 == "3")]


ID_2_3_2 <- ID_2_3[-which(ID_2_3 == "A7"#)]
                          &AM_2_3 == "3")]

AM_2_3 <- AM_2_3[-which(ID_2_3 == "A7"#)]
                        & AM_2_3 == "3")]

ID_2_3 <- ID_2_3_2








com_SS2_SS3 <- com_SS2_SS3[which(ID_2_3 != "A4" #),]
                                   &ID_2_3 != "A1"),]

fish_isolation <- fish_isolation[which(ID_2_3 != "A4" #)]
                                      & ID_2_3 != "A1")]

fish <- fish[which(ID_2_3 != "A4" #)]
                                    &   ID_2_3 != "A1")]

All <- All[which(ID_2_3 != "A4" #)]
                     & ID_2_3 != "A1")]

isolation_Exp1 <- isolation_Exp1[which(ID_2_3 != "A4"#)]
                   & ID_2_3 != "A1")]

AM_2_3 <- AM_2_3[which(ID_2_3 != "A4"#)]
                                        & ID_2_3 != "A1")]

ID_2_3 <- ID_2_3[which(ID_2_3 != "A4"#)]
                       & ID_2_3 != "A1")]
#######



####################################################### #Abundance ####
com_SS2_SS3_abundance_LOG <- log(rowSums(com_SS2_SS3))
com_SS2_SS3_abundance <- rowSums(com_SS2_SS3)

isolation_SS <-rep(NA, length(All))
for(i in 1:length(All)){
  if(AM_2_3[i] == "2" & isolation_Exp1[i] == "030"){isolation_SS[i] <- "2 030"}
  if(AM_2_3[i] == "2" & isolation_Exp1[i] == "120"){isolation_SS[i] <- "2 120"}
  if(AM_2_3[i] == "2" & isolation_Exp1[i] == "480"){isolation_SS[i] <- "2 480"}
  if(AM_2_3[i] == "3" & isolation_Exp1[i] == "030"){isolation_SS[i] <- "3 030"}
  if(AM_2_3[i] == "3" & isolation_Exp1[i] == "120"){isolation_SS[i] <- "3 120"}
  if(AM_2_3[i] == "3" & isolation_Exp1[i] == "480"){isolation_SS[i] <- "3 480"}
}
isolation_SS <- as.factor(isolation_SS)
isolation_SS <- factor(isolation_SS, levels = c("2 030", "3 030","2 120", "3 120","2 480", "3 480" ))
new_All <- factor(All, levels = c("2 030 absent","3 030 absent","2 120 absent","3 120 absent","2 480 absent","3 480 absent",
                                  "2 030 present","3 030 present","2 120 present","3 120 present","2 480 present","3 480 present"))


svg(filename = "Results/Abundance.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(com_SS2_SS3_abundance~isolation_SS, outline = F, ylab = "Abundance", xlab = "", at = c(1,2,4,5,7,8), lwd = 1.5, ylim = c(0,1400))
mylevels <- levels(new_All)
levelProportions <- summary(new_All)/length(dist_Exp1_deviation$observed_distances)
col <-  c(rep(c("sienna3","grey70"),3), rep(c("dodgerblue3","grey70"),3))
bg <- c(rep("sienna3",6), rep("dodgerblue3",6))
pch <- c(15,22,16,21,17,24,15,22,16,21,17,24)
for(i in 1:length(mylevels)){
  
  x <- c(1,2,4,5,7,8,1,2,4,5,7,8)[i]
  thislevel <- mylevels[i]
  thisvalues <- com_SS2_SS3_abundance[new_All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i], cex = 1.5, lwd = 3) 
  
}
boxplot(com_SS2_SS3_abundance~isolation_SS, add = T, col = "transparent", outline = F,at = c(1,2,4,5,7,8), lwd = 1.5)
box(lwd = 2.5)

dev.off()

mix_mod_ab <- lmer(com_SS2_SS3_abundance_LOG~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))

round(Anova(mix_mod_ab, test.statistic = "F"),3)
plot(mix_mod_ab)
qqnorm(resid(mix_mod_ab))
qqline(resid(mix_mod_ab))
boxplot(resid(mix_mod_ab, method = "normalized") ~ fish*isolation_Exp1*AM_2_3)


emmeans(mix_mod_ab, list(pairwise ~ isolation_Exp1|AM_2_3), adjust = "tukey")
emmeans(mix_mod_ab, list(pairwise ~ AM_2_3|isolation_Exp1), adjust = "tukey")


####################################################
mix_mod_ab_nb <- glmer.nb(com_SS2_SS3_abundance~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))

round(Anova(mix_mod_ab_nb, test.statistic = "Chisq"),3)
plot(mix_mod_ab_nb)
qqnorm(resid(mix_mod_ab_nb, method = "pearson"))
qqline(resid(mix_mod_ab_nb, method = "pearson"))
boxplot(resid(mix_mod_ab_nb, method = "pearson") ~ fish*isolation_Exp1*AM_2_3)

emmeans(mix_mod_ab_nb, list(pairwise ~ isolation_Exp1|AM_2_3), adjust = "tukey")
emmeans(mix_mod_ab_nb, list(pairwise ~ AM_2_3|isolation_Exp1), adjust = "tukey")




set.seed(2);pairwise_30_survey <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(isolation_Exp1 == "030")],AM_2_3[which(isolation_Exp1 == "030")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_120_survey <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(isolation_Exp1 == "120")],AM_2_3[which(isolation_Exp1 == "120")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_480_survey <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(isolation_Exp1 == "480")],AM_2_3[which(isolation_Exp1 == "480")], p.method = "none", nperm = 1000)

p.adjust(c(pairwise_30_survey$p.value,
           pairwise_120_survey$p.value,
           pairwise_480_survey$p.value), method = "fdr")



pairwise_SS2_isolation <- pairwise.t.test(com_SS2_SS3_abundance_LOG[which(AM_2_3 == "2")],isolation_Exp1[which(AM_2_3 == "2")], p.adjust.method = "none")
pairwise_SS3_isolation <- pairwise.t.test(com_SS2_SS3_abundance_LOG[which(AM_2_3 == "3")],isolation_Exp1[which(AM_2_3 == "3")], p.adjust.method = "none")


pairwise_SS2_isolation <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(AM_2_3 == "2")],isolation_Exp1[which(AM_2_3 == "2")], p.method = "none", nperm = 10000)
pairwise_SS3_isolation <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(AM_2_3 == "3")],isolation_Exp1[which(AM_2_3 == "3")], p.method = "none", nperm = 10000)



########################################33


#####


###################Richness####



X <- iNEXT(t(com_SS2_SS3), datatype = "abundance", q = 0,  knots = 100,se = T, conf = 0.95,nboot =100)

rare_rich <- X$AsyEst[which(X$AsyEst[,2] == "Species richness"),][,4]
rich <- X$AsyEst[which(X$AsyEst[,2] == "Species richness"),][,3]


svg(filename = "Results/Richness.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(rich~fish, outline = F, ylab = "Richness", xlab = "", at = c(1,2), lwd = 1.5)
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,1,1,2,2,2,1,1,1,2,2,2)[i]
  thislevel <- mylevels[i]
  thisvalues <- rich[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.3)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(rich~fish, add = T, col = "transparent", outline = F,at = c(1,2), lwd = 1.5)
box(lwd = 2.5)
dev.off()




fit_rich_rand <- lmer(rich~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
fit_rich <- lm(rich~fish*isolation_Exp1*AM_2_3)

rePCA(fit_rich_rand)
getME(fit_rich_rand, name ="ALL")

round(Anova(fit_rich_rand, test.statistic = "F"),3)
plot(fit_rich_rand)
qqnorm(resid(fit_rich_rand))
qqline(resid(fit_rich_rand))

round(Anova(fit_rich, test.statistic = "F"),3)
plot(fit_rich)
qqnorm(resid(fit_rich))
qqline(resid(fit_rich))




svg(filename = "Results/Rarefied Richness.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(rare_rich~fish_isolation, outline = F, ylab = "Richness", xlab = "", at = c(1,2,3,5,6,7), lwd = 1.5)
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(16,16,16,15,15,15,21,21,21,22,22,22)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- rare_rich[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(rare_rich~fish_isolation, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5)
box(lwd = 2.5)
dev.off()

fit_rare_rich <- lmer(rare_rich~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3), control = lmerControl(optimizer = "Nelder_Mead"))
round(Anova(fit_rare_rich),3)
estimated_means_rare_rich <- get_means(fit, "fish*isolation_Exp1*AM_2_3")
estimated_means_rare_rich




####################################################################################################################################
################################################ GAMMA TWO SAMPLES TOGETHER ########################################################  


## example for incidence frequencies based data (list of data.frame)
absent_30<-as.incfreq(t(decostand(com_SS2_SS3[which(fish_isolation == "030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS2_SS3[which(fish_isolation == "030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS2_SS3[which(fish_isolation == "120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS2_SS3[which(fish_isolation == "120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS2_SS3[which(fish_isolation == "480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS2_SS3[which(fish_isolation == "480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])



svg(filename = "Gamma SS2 SS3.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(10,30), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
arrows(y0 = c(Gamma$qD.LCL[1],Gamma$qD.LCL[2], Gamma$qD.LCL[3],Gamma$qD.LCL[4],Gamma$qD.LCL[5], Gamma$qD.LCL[6]),
       y1 = c(Gamma$qD.UCL[1],Gamma$qD.UCL[2], Gamma$qD.UCL[3],Gamma$qD.UCL[4],Gamma$qD.UCL[5], Gamma$qD.UCL[6]),
       x1 = c(1,2,3,5,6,7), x0 = c(1,2,3,5,6,7), 
       code = 3, angle = 90, length = 0.05, c("grey50"), lwd = 2)
points(y = c(Gamma$qD[1],Gamma$qD[2], Gamma$qD[3],Gamma$qD[4],Gamma$qD[5], Gamma$qD[6]), x = c(1,2,3,5,6,7), cex = 2, pch = c(15,16,17), col = c("sienna1","sienna3","sienna4","dodgerblue1","dodgerblue3", "dodgerblue4"))
axis(1, at = c(1,2,3,5,6,7), labels = c("30","120","480","30","120","480"), las = 1, cex.axis = 1)
box(lwd = 2.5)

dev.off()







####################################################################################################################################
################################################################# GAMMA SS1 ########################################################  


absent_30<-as.incfreq(t(decostand(com_SS1_Exp1[which(fish_isolation_SS1 == "1 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS1_Exp1[which(fish_isolation_SS1 == "1 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS1_Exp1[which(fish_isolation_SS1 == "1 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS1_Exp1[which(fish_isolation_SS1 == "1 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS1_Exp1[which(fish_isolation_SS1 == "1 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS1_Exp1[which(fish_isolation_SS1 == "1 480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])


svg(filename = "Gamma SS1.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(0,20), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
arrows(y0 = c(Gamma$qD.LCL[1],Gamma$qD.LCL[2], Gamma$qD.LCL[3],Gamma$qD.LCL[4],Gamma$qD.LCL[5], Gamma$qD.LCL[6]),
       y1 = c(Gamma$qD.UCL[1],Gamma$qD.UCL[2], Gamma$qD.UCL[3],Gamma$qD.UCL[4],Gamma$qD.UCL[5], Gamma$qD.UCL[6]),
       x1 = c(1,2,3,5,6,7), x0 = c(1,2,3,5,6,7), 
       code = 3, angle = 90, length = 0.05, c("grey50"), lwd = 2)
points(y = c(Gamma$qD[1],Gamma$qD[2], Gamma$qD[3],Gamma$qD[4],Gamma$qD[5], Gamma$qD[6]), x = c(1,2,3,5,6,7), cex = 2, pch = c(15,16,17), col = c("sienna1","sienna3","sienna4","dodgerblue1","dodgerblue3", "dodgerblue4"))
axis(1, at = c(1,2,3,5,6,7), labels = c("30","120","480","30","120","480"), las = 1, cex.axis = 1)
box(lwd = 2.5)

dev.off()




####################################################################################################################################
################################################################# GAMMA SS2 ########################################################  


absent_30<-as.incfreq(t(decostand(com_SS2_Exp1[which(fish_isolation_SS2 == "2 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS2_Exp1[which(fish_isolation_SS2 == "2 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS2_Exp1[which(fish_isolation_SS2 == "2 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS2_Exp1[which(fish_isolation_SS2 == "2 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS2_Exp1[which(fish_isolation_SS2 == "2 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS2_Exp1[which(fish_isolation_SS2 == "2 480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])



svg(filename = "Gamma SS2.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(10,30), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
arrows(y0 = c(Gamma$qD.LCL[1],Gamma$qD.LCL[2], Gamma$qD.LCL[3],Gamma$qD.LCL[4],Gamma$qD.LCL[5], Gamma$qD.LCL[6]),
       y1 = c(Gamma$qD.UCL[1],Gamma$qD.UCL[2], Gamma$qD.UCL[3],Gamma$qD.UCL[4],Gamma$qD.UCL[5], Gamma$qD.UCL[6]),
       x1 = c(1,2,3,5,6,7), x0 = c(1,2,3,5,6,7), 
       code = 3, angle = 90, length = 0.05, c("grey50"), lwd = 2)
points(y = c(Gamma$qD[1],Gamma$qD[2], Gamma$qD[3],Gamma$qD[4],Gamma$qD[5], Gamma$qD[6]), x = c(1,2,3,5,6,7), cex = 2, pch = c(15,16,17), col = c("sienna1","sienna3","sienna4","dodgerblue1","dodgerblue3", "dodgerblue4"))
axis(1, at = c(1,2,3,5,6,7), labels = c("30","120","480","30","120","480"), las = 1, cex.axis = 1)
box(lwd = 2.5)

dev.off()



####################################################################################################################################
################################################################# GAMMA SS3 ########################################################  


absent_30<-as.incfreq(t(decostand(com_SS3_Exp1[which(fish_isolation_SS3 == "3 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS3_Exp1[which(fish_isolation_SS3 == "3 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS3_Exp1[which(fish_isolation_SS3 == "3 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS3_Exp1[which(fish_isolation_SS3 == "3 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS3_Exp1[which(fish_isolation_SS3 == "3 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS3_Exp1[which(fish_isolation_SS3 == "3 480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])


svg(filename = "Gamma SS3.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(10,30), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
arrows(y0 = c(Gamma$qD.LCL[1],Gamma$qD.LCL[2], Gamma$qD.LCL[3],Gamma$qD.LCL[4],Gamma$qD.LCL[5], Gamma$qD.LCL[6]),
       y1 = c(Gamma$qD.UCL[1],Gamma$qD.UCL[2], Gamma$qD.UCL[3],Gamma$qD.UCL[4],Gamma$qD.UCL[5], Gamma$qD.UCL[6]),
       x1 = c(1,2,3,5,6,7), x0 = c(1,2,3,5,6,7), 
       code = 3, angle = 90, length = 0.05, c("grey50"), lwd = 2)
points(y = c(Gamma$qD[1],Gamma$qD[2], Gamma$qD[3],Gamma$qD[4],Gamma$qD[5], Gamma$qD[6]), x = c(1,2,3,5,6,7), cex = 2, pch = c(15,16,17), col = c("sienna1","sienna3","sienna4","dodgerblue1","dodgerblue3", "dodgerblue4"))
axis(1, at = c(1,2,3,5,6,7), labels = c("30","120","480","30","120","480"), las = 1, cex.axis = 1)
box(lwd = 2.5)

dev.off()






