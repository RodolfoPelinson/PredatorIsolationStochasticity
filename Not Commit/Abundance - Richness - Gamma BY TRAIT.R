library(lmPerm)
library(iNEXT)
library(lme4)
library(psycho)
library(tidyverse)
library(car)
library(emmeans)

####################################################### #Abundance ####
com_SS2_SS3_predators <- com_SS2_SS3[,which(Trait_SS2_SS3$trophic == "Pr")]

com_SS2_SS3_abundance_LOG <- log(rowSums(com_SS2_SS3_predators))
com_SS2_SS3_abundance <- rowSums(com_SS2_SS3_predators)

svg(filename = "Results/Predators/Abundance.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(com_SS2_SS3_abundance~fish_isolation, outline = F, ylab = "Abundance", xlab = "", at = c(1,2,3,5,6,7), lwd = 1.5, ylim = c(0,170))
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){

  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- com_SS2_SS3_abundance[All==thislevel]

  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3)

}
boxplot(com_SS2_SS3_abundance~fish_isolation, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5)
box(lwd = 2.5)
dev.off()


mix_mod_ab <- lmer(com_SS2_SS3_abundance_LOG~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
round(Anova(mix_mod_ab, test.statistic = "F"),3)
plot(mix_mod_ab)
qqnorm(resid(mix_mod_ab))
qqline(resid(mix_mod_ab))

set.seed(2);pairwise_isolation1 <- pairwise.perm.t.test(com_SS2_SS3_abundance,isolation_Exp1, p.method = "none", nperm = 1000)
set.seed(2);pairwise_isolation2 <- pairwise.perm.t.test(com_SS2_SS3_abundance,isolation_Exp1, p.method = "fdr", nperm = 1000)

set.seed(2);pairwise_SS_fishless <- pairwise.perm.t.test(com_SS2_SS3_abundance[which(fish == "absent")] ,AM_2_3[which(fish == "absent")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_SS_fish <- pairwise.perm.t.test(com_SS2_SS3_abundance[which(fish == "present")] ,AM_2_3[which(fish == "present")], p.method = "none", nperm = 1000)

p.adjust(c(pairwise_SS_fishless$p.value,
           pairwise_SS_fish$p.value), method = "fdr")



###################Richness####



X <- iNEXT(t(com_SS2_SS3_predators), datatype = "abundance", q = 0,  knots = 100,se = T, conf = 0.95,nboot =100)

rare_rich <- X$AsyEst[which(X$AsyEst[,2] == "Species richness"),][,4]
rich <- X$AsyEst[which(X$AsyEst[,2] == "Species richness"),][,3]


svg(filename = "Results/Predators/Richness.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(rich~fish_isolation, outline = F, ylab = "Richness", xlab = "", at = c(1,2,3,5,6,7), lwd = 1.5)
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){

  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- rich[All==thislevel]

  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3)

}
boxplot(rich~fish_isolation, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5)
box(lwd = 2.5)
dev.off()




fit_rich_rand <- lmer(rich~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
round(Anova(fit_rich_rand, test.statistic = "F"),3)
plot(fit_rich_rand)
qqnorm(resid(fit_rich_rand))
qqline(resid(fit_rich_rand))


set.seed(2);pairwise_isolation1 <- pairwise.perm.t.test(rich,isolation_Exp1, p.method = "none", nperm = 1000)
set.seed(2);pairwise_isolation2 <- pairwise.perm.t.test(rich,isolation_Exp1, p.method = "fdr", nperm = 1000)







####################################################################################################################################
################################################################# GAMMA SS2 ########################################################
com_SS2_predators <- com_SS2[,which(Trait_SS2$trophic == "Pr")]


absent_30<-as.incfreq(t(decostand(com_SS2_predators[which(fish_isolation_SS2 == "2 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS2_predators[which(fish_isolation_SS2 == "2 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS2_predators[which(fish_isolation_SS2 == "2 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS2_predators[which(fish_isolation_SS2 == "2 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS2_predators[which(fish_isolation_SS2 == "2 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS2_predators[which(fish_isolation_SS2 == "2 480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])



svg(filename = "Gamma SS2 pred.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(0,12), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
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
com_SS3_predators <- com_SS3[,which(Trait_SS3$trophic == "Pr")]


absent_30<-as.incfreq(t(decostand(com_SS3_predators[which(fish_isolation_SS3 == "3 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS3_predators[which(fish_isolation_SS3 == "3 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS3_predators[which(fish_isolation_SS3 == "3 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS3_predators[which(fish_isolation_SS3 == "3 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS3_predators[which(fish_isolation_SS3 == "3 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS3_predators[which(fish_isolation_SS3 == "3 480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])



svg(filename = "Gamma SS3 pred.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(0,12), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
arrows(y0 = c(Gamma$qD.LCL[1],Gamma$qD.LCL[2], Gamma$qD.LCL[3],Gamma$qD.LCL[4],Gamma$qD.LCL[5], Gamma$qD.LCL[6]),
       y1 = c(Gamma$qD.UCL[1],Gamma$qD.UCL[2], Gamma$qD.UCL[3],Gamma$qD.UCL[4],Gamma$qD.UCL[5], Gamma$qD.UCL[6]),
       x1 = c(1,2,3,5,6,7), x0 = c(1,2,3,5,6,7),
       code = 3, angle = 90, length = 0.05, c("grey50"), lwd = 2)
points(y = c(Gamma$qD[1],Gamma$qD[2], Gamma$qD[3],Gamma$qD[4],Gamma$qD[5], Gamma$qD[6]), x = c(1,2,3,5,6,7), cex = 2, pch = c(15,16,17), col = c("sienna1","sienna3","sienna4","dodgerblue1","dodgerblue3", "dodgerblue4"))
axis(1, at = c(1,2,3,5,6,7), labels = c("30","120","480","30","120","480"), las = 1, cex.axis = 1)
box(lwd = 2.5)

dev.off()












#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#Preys



####################################################### #Abundance ####
com_SS2_SS3_preys <- com_SS2_SS3[,which(Trait_SS2_SS3$trophic == "Non_Pred")]

com_SS2_SS3_non_predators_abundance <- rowSums(com_SS2_SS3_preys)


com_SS2_SS3_abundance_LOG <- log(rowSums(com_SS2_SS3_preys))


svg(filename = "Results/Preys/Abundance.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
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

mix_mod_ab <- lmer(com_SS2_SS3_abundance~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
round(Anova(mix_mod_ab, test.statistic = "F"),3)
plot(mix_mod_ab)
qqnorm(resid(mix_mod_ab))
qqline(resid(mix_mod_ab))

library(MASS)
mod_ab_consumers <- glmer.nb(com_SS2_SS3_abundance~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
Anova(mod_ab_consumers, test.statistic = "Chisq")

set.seed(2);pairwise_30_survey <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(isolation_Exp1 == "030")],AM_2_3[which(isolation_Exp1 == "030")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_120_survey <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(isolation_Exp1 == "120")],AM_2_3[which(isolation_Exp1 == "120")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_480_survey <- pairwise.perm.t.test(com_SS2_SS3_abundance_LOG[which(isolation_Exp1 == "480")],AM_2_3[which(isolation_Exp1 == "480")], p.method = "none", nperm = 1000)

p.adjust(c(pairwise_30_survey$p.value,
           pairwise_120_survey$p.value,
           pairwise_480_survey$p.value), method = "fdr")




###################Richness####



X <- iNEXT(t(com_SS2_SS3_preys), datatype = "abundance", q = 0,  knots = 100,se = T, conf = 0.95,nboot =100)

rare_rich <- X$AsyEst[which(X$AsyEst[,2] == "Species richness"),][,4]
rich <- X$AsyEst[which(X$AsyEst[,2] == "Species richness"),][,3]


svg(filename = "Results/Preys/Richness.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
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
round(Anova(fit_rich_rand, test.statistic = "F"),3)
plot(fit_rich_rand)
qqnorm(resid(fit_rich_rand))
qqline(resid(fit_rich_rand))















####################################################################################################################################
################################################################# GAMMA SS2 ########################################################
com_SS2_preys <- com_SS2[,which(Trait_SS2$trophic == "Non_Pred")]


absent_30<-as.incfreq(t(decostand(com_SS2_preys[which(fish_isolation_SS2 == "2 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS2_preys[which(fish_isolation_SS2 == "2 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS2_preys[which(fish_isolation_SS2 == "2 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS2_preys[which(fish_isolation_SS2 == "2 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS2_preys[which(fish_isolation_SS2 == "2 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS2_preys[which(fish_isolation_SS2 == "2 480 present"),], method = "pa")))

all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])



svg(filename = "Gamma SS2 prey.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(5,18), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
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
com_SS3_preys <- com_SS3[,which(Trait_SS3$trophic == "Non_Pred")]


absent_30<-as.incfreq(t(decostand(com_SS3_preys[which(fish_isolation_SS3 == "3 030 absent"),], method = "pa")))
present_30<-as.incfreq(t(decostand(com_SS3_preys[which(fish_isolation_SS3 == "3 030 present"),], method = "pa")))
absent_120<-as.incfreq(t(decostand(com_SS3_preys[which(fish_isolation_SS3 == "3 120 absent"),], method = "pa")))
present_120<-as.incfreq(t(decostand(com_SS3_preys[which(fish_isolation_SS3 == "3 120 present"),], method = "pa")))
absent_480<-as.incfreq(t(decostand(com_SS3_preys[which(fish_isolation_SS3 == "3 480 absent"),], method = "pa")))
present_480<-as.incfreq(t(decostand(com_SS3_preys[which(fish_isolation_SS3 == "3 480 present"),], method = "pa")))
all <- list(absent_30=absent_30 ,absent_120 =absent_120, absent_480=absent_480, present_30=present_30,present_120=present_120,present_480=present_480)

X_ALL <- iNEXT(all, datatype = "incidence_freq", q = 0, knots = 40,se = T, conf = 0.95,nboot =1000, size = c(1:4))


Gamma <- rbind(absent_30 = X_ALL$iNextEst$absent_30[4,c(4,5,6)],
               absent_120 = X_ALL$iNextEst$absent_120[4,c(4,5,6)],
               absent_480 = X_ALL$iNextEst$absent_480[4,c(4,5,6)],
               present_30 = X_ALL$iNextEst$present_30[4,c(4,5,6)],
               present_120 = X_ALL$iNextEst$present_120[4,c(4,5,6)],
               present_480 = X_ALL$iNextEst$present_480[4,c(4,5,6)])



svg(filename = "Gamma SS3 prey.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
plot(c(NA,NA,NA,NA,NA,NA),ylim = c(5,18), xlim = c(0.5,7.5),type = "p", xaxt = "n",ylab = "Gamma",xlab = "Treatment", cex.lab = 1, cex.axis = 1)
arrows(y0 = c(Gamma$qD.LCL[1],Gamma$qD.LCL[2], Gamma$qD.LCL[3],Gamma$qD.LCL[4],Gamma$qD.LCL[5], Gamma$qD.LCL[6]),
       y1 = c(Gamma$qD.UCL[1],Gamma$qD.UCL[2], Gamma$qD.UCL[3],Gamma$qD.UCL[4],Gamma$qD.UCL[5], Gamma$qD.UCL[6]),
       x1 = c(1,2,3,5,6,7), x0 = c(1,2,3,5,6,7),
       code = 3, angle = 90, length = 0.05, c("grey50"), lwd = 2)
points(y = c(Gamma$qD[1],Gamma$qD[2], Gamma$qD[3],Gamma$qD[4],Gamma$qD[5], Gamma$qD[6]), x = c(1,2,3,5,6,7), cex = 2, pch = c(15,16,17), col = c("sienna1","sienna3","sienna4","dodgerblue1","dodgerblue3", "dodgerblue4"))
axis(1, at = c(1,2,3,5,6,7), labels = c("30","120","480","30","120","480"), las = 1, cex.axis = 1)
box(lwd = 2.5)

dev.off()
