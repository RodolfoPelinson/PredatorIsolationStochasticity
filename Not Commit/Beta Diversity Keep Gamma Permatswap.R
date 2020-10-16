library(AICcmodavg)
library(lme4)
library(car)
library(RVAideMemoire)
library(emmeans)


###################################Diversidade Beta


dist_Exp1_deviation <- beta_deviation(com_SS2_SS3, strata = All, times = 1000,
                                      transform = NULL, dist = "bray", fixedmar="both",
                                      shuffle = "both", method = "quasiswap", seed = 2, group = All) #Function that calculates the dissimilarity, distances to centroid and deviations


####################################################### Observed #################################################### 
fit_observed <- lmer(dist_Exp1_deviation$observed_distances~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
round(Anova(fit_observed, test.statistic = "F", type = "II"),3)
plot(fit_observed)
qqnorm(resid(fit_observed))
qqline(resid(fit_observed))


emmeans(fit_observed, list(pairwise ~ isolation_Exp1|fish), adjust = "tukey")


#Isolation comparison
set.seed(2);pairwise_observed_isolation_absent <- pairwise.perm.t.test(dist_Exp1_deviation$observed_distances[which(fish == "absent")],
                                                      isolation_Exp1[which(fish == "absent")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_observed_isolation_present <- pairwise.perm.t.test(dist_Exp1_deviation$observed_distances[which(fish == "present")],
                                                       isolation_Exp1[which(fish == "present")], p.method = "none", nperm = 1000)

p.adjust(na.omit(c(pairwise_observed_isolation_absent$p.value,pairwise_observed_isolation_present$p.value )), method = "fdr")


####################################################### Expected #################################################### 
fit_expected_rand <- lmer(dist_Exp1_deviation$expected_distances~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))
fit_expected <- lm(dist_Exp1_deviation$expected_distances~fish*isolation_Exp1*AM_2_3)

round(Anova(fit_expected_rand, test.statistic = "F"),3)
plot(fit_expected_rand)
qqnorm(resid(fit_expected_rand))
qqline(resid(fit_expected_rand))

#Isolation comparison
set.seed(3);pairwise_expected_isolation_absent <- pairwise.perm.t.test(dist_Exp1_deviation$expected_distances[which(fish == "absent")],
                                                                       isolation_Exp1[which(fish == "absent")], p.method = "none", nperm = 1000)
set.seed(3);pairwise_expected_isolation_present <- pairwise.perm.t.test(dist_Exp1_deviation$expected_distances[which(fish == "present")],
                                                                        isolation_Exp1[which(fish == "present")], p.method = "none", nperm = 1000)

p.adjust(na.omit(c(pairwise_expected_isolation_absent$p.value,pairwise_expected_isolation_present$p.value )), method = "fdr")




# Sampling by isolation
pairwise_SS_expected_30 <- pairwise.perm.t.test(dist_Exp1_deviation$expected_distances[which(isolation_Exp1 == "030")],
                                        AM_2_3[which(isolation_Exp1 == "030")], p.method = "none", nperm = 1000)
pairwise_SS_expected_120 <- pairwise.perm.t.test(dist_Exp1_deviation$expected_distances[which(isolation_Exp1 == "120")],
                                         AM_2_3[which(isolation_Exp1 == "120")], p.method = "none", nperm = 1000)
pairwise_SS_expected_480 <- pairwise.perm.t.test(dist_Exp1_deviation$expected_distances[which(isolation_Exp1 == "480")],
                                         AM_2_3[which(isolation_Exp1 == "480")], p.method = "none", nperm = 1000)

p.adjust(c(pairwise_expected_30$p.value, pairwise_expected_120$p.value, pairwise_expected_480$p.value), method = "fdr")




####################################################### Deviation #################################################### 

fit_deviation_rand <- lmer(dist_Exp1_deviation$deviation_distances~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3),control = lmerControl(optimizer = c("bobyqa")))
fit_deviation_rand <- lmer(dist_Exp1_deviation$deviation_distances~fish*isolation_Exp1*AM_2_3 + (1|ID_2_3))

fit_deviation <- lm(dist_Exp1_deviation$deviation_distances~fish*isolation_Exp1*AM_2_3)
plot(fit_expected_rand)
qqnorm(resid(fit_expected_rand))
qqline(resid(fit_expected_rand))

round(Anova(fit_deviation_rand, test.statistic = "F"),3)
round(Anova(fit_deviation, test.statistic = "F"),3)

emmeans(fit_deviation_rand, list(pairwise ~ isolation_Exp1|fish), adjust = "tukey")


#Isolation comparison
set.seed(2);pairwise_deviation_isolation_absent <- pairwise.perm.t.test(dist_Exp1_deviation$deviation_distances[which(fish == "absent")],
                                                      isolation_Exp1[which(fish == "absent")], p.method = "none", nperm = 1000)
set.seed(2);pairwise_deviation_isolation_present <- pairwise.perm.t.test(dist_Exp1_deviation$deviation_distances[which(fish == "present")],
                                                       isolation_Exp1[which(fish == "present")], p.method = "none", nperm = 1000)

p.adjust(na.omit(c(pairwise_deviation_isolation_absent$p.value,pairwise_deviation_isolation_present$p.value )), method = "fdr")





# Basic boxplot
svg(filename = "Results Keep Gamma permatswap/Beta Observed.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(dist_Exp1_deviation$observed_distances~fish_isolation, outline = F, ylab = "Distance to Centroid (Observed)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5)
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- dist_Exp1_deviation$observed_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(dist_Exp1_deviation$observed_distances~fish_isolation, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5)
box(lwd = 2.5)
dev.off()


svg(filename = "Results Keep Gamma permatswap/Beta Expected.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(dist_Exp1_deviation$expected_distances~fish_isolation, outline = F, ylab = "Distance to Centroid (Expected)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5)
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$expected_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- dist_Exp1_deviation$expected_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(dist_Exp1_deviation$expected_distances~fish_isolation, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5)
box(lwd = 2.5)
dev.off()


svg(filename = "Results Keep Gamma permatswap/Beta Deviation.svg", width = 5, height = 5, pointsize = 12, bg = "transparent")
boxplot(dist_Exp1_deviation$deviation_distances~fish_isolation, outline = F, ylab = "Distance to Centroid (Deviation)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(-2,10), lwd = 1.5)
mylevels <- levels(All)
levelProportions <- summary(All)/length(dist_Exp1_deviation$deviation_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- dist_Exp1_deviation$deviation_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(dist_Exp1_deviation$deviation_distances~fish_isolation, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5)
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
box(lwd = 2.5)
dev.off()









se <- function(x){
  sd(x)/sqrt(length(x))
}



means_observed <- tapply(dist_Exp1_deviation$observed_distances, fish_isolation, mean)
means_observed <- means_observed[c(1,3,5,2,4,6)]
se_observed <- tapply(dist_Exp1_deviation$observed_distances, fish_isolation, se)
se_observed <- se_observed[c(1,3,5,2,4,6)]

plot(c(1,1.5,2,3,3.5,4), means_observed,xlim =c(0.5,4.5), ylim = c(0.2,0.6))
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
arrows(x0 = c(1,1.5,2,3,3.5,4), x1 = c(1,1.5,2,3,3.5,4), y0 = (means_observed + se_observed), y1 = (means_observed - se_observed), length = 0, lwd = 2)
points(c(1,1.5,2,3,3.5,4), means_observed, col = c(rep("sienna3",3),rep("dodgerblue3",3)), pch = c(rep(15,3),rep(16,3)), cex = 2)
box(lwd = 2.5)




means_expected <- tapply(dist_Exp1_deviation$expected_distances, fish_isolation, mean)
means_expected <- means_expected[c(1,3,5,2,4,6)]
se_expected <- tapply(dist_Exp1_deviation$expected_distances, fish_isolation, se)
se_expected <- se_expected[c(1,3,5,2,4,6)]

plot(c(1,1.5,2,3,3.5,4), means_expected,xlim =c(0.5,4.5), ylim = c(0.2,0.6))
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
arrows(x0 = c(1,1.5,2,3,3.5,4), x1 = c(1,1.5,2,3,3.5,4), y0 = (means_expected + se_expected), y1 = (means_expected - se_expected), length = 0, lwd = 2)
points(c(1,1.5,2,3,3.5,4), means_expected, col = c(rep("sienna3",3),rep("dodgerblue3",3)), pch = c(rep(15,3),rep(16,3)), cex = 2)
box(lwd = 2.5)



means_deviation <- tapply(dist_Exp1_deviation$deviation_distances, fish_isolation, mean)
means_deviation <- means_deviation[c(1,3,5,2,4,6)]
se_deviation <- tapply(dist_Exp1_deviation$deviation_distances, fish_isolation, se)
se_deviation <- se_deviation[c(1,3,5,2,4,6)]

plot(c(1,1.5,2,3,3.5,4), means_deviation,xlim =c(0.5,4.5), ylim = c(-1,4))
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
arrows(x0 = c(1,1.5,2,3,3.5,4), x1 = c(1,1.5,2,3,3.5,4), y0 = (means_deviation + se_deviation), y1 = (means_deviation - se_deviation), length = 0, lwd = 2)
points(c(1,1.5,2,3,3.5,4), means_deviation, col = c(rep("sienna3",3),rep("dodgerblue3",3)), pch = c(rep(15,3),rep(16,3)), cex = 2)
box(lwd = 2.5)













