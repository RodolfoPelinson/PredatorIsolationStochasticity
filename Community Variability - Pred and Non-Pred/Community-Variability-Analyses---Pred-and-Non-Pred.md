Community Variability - Predatory and Non-Predatory Insects
================
Rodolfo Pelinson
20/10/2020

This is the same community variability analyses presented in the main
paper, but separating the communities into predatory and non-predatory
insects.

If you haven’t, install the package:

``` r
install.packages("devtools")
devtools::install_github("RodolfoPelinson/Pelinson.et.al.2020B")
```

These are the packages you will need to run this code:

``` r
library(Pelinson.et.al.2020B)
library(lme4) # Version 1.1-23
library(car) # Version 3.0-7
library(emmeans) # Version 1.4.8
library(vegan) # Version 2.5-6
```

## Community Variability

### Considering only Predatory insects for the last two surveys.

First loading data

``` r
data(com_SS2_SS3_predators, All, fish_SS2_SS3, isolation_SS2_SS3, SS_SS2_SS3, ID_SS2_SS3, fish_isolation_SS2_SS3)
```

    ## Warning in data(com_SS2_SS3_predators, All, fish_SS2_SS3, isolation_SS2_SS3, :
    ## data set 'fish_isolation_SS2_SS3' not found

Computing observed and expected distances to centroid, and
beta-deviation.

``` r
beta_deviation_SS2_SS3_predators <- beta_deviation(com_SS2_SS3_predators, strata = All, times = 10000,
                                      transform = NULL, dist = "bray", fixedmar="both",
                                      shuffle = "both", method = "quasiswap", seed = 2, group = All) 
```

#### Observed Community Variability

Running ANOVA table for observed distances to group centroids, or
observed beta-diversity/community variability in each treatment.

``` r
fit_observed_SS2_SS3 <- lmer(beta_deviation_SS2_SS3_predators$observed_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
round(Anova(fit_observed_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3_predators$observed_distances
    ##                                           Chisq Df Pr(>Chisq)   
    ## fish_SS2_SS3                              1.199  1      0.273   
    ## isolation_SS2_SS3                         2.129  2      0.345   
    ## SS_SS2_SS3                                5.533  1      0.019 * 
    ## fish_SS2_SS3:isolation_SS2_SS3            6.370  2      0.041 * 
    ## fish_SS2_SS3:SS_SS2_SS3                   9.009  1      0.003 **
    ## isolation_SS2_SS3:SS_SS2_SS3              7.528  2      0.023 * 
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 3.697  2      0.157   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_observed_SS2_SS3, list(pairwise ~ fish_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of fish_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  fish_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  absent        0.329 0.0339 32    0.249    0.409
    ##  present       0.273 0.0339 32    0.193    0.352
    ## 
    ## SS_SS2_SS3 = 3:
    ##  fish_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  absent        0.315 0.0360 32    0.231    0.400
    ##  present       0.473 0.0399 32    0.379    0.566
    ## 
    ## Results are averaged over the levels of: isolation_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of fish_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast         estimate     SE df t.ratio p.value
    ##  absent - present   0.0562 0.0480 32  1.171  0.2501 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast         estimate     SE df t.ratio p.value
    ##  absent - present  -0.1573 0.0538 32 -2.927  0.0063 
    ## 
    ## Results are averaged over the levels of: isolation_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger

``` r
emmeans(fit_observed_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.323 0.0416 32    0.218    0.428
    ##  120                0.313 0.0416 32    0.208    0.418
    ##  480                0.266 0.0416 32    0.162    0.371
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.287 0.0454 32    0.172    0.401
    ##  120                0.419 0.0454 32    0.305    0.533
    ##  480                0.477 0.0489 32    0.353    0.600
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120    0.0100 0.0588 32  0.171  0.9976 
    ##  30 - 480    0.0566 0.0588 32  0.962  0.7165 
    ##  120 - 480   0.0465 0.0588 32  0.792  0.8191 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   -0.1324 0.0641 32 -2.064  0.1350 
    ##  30 - 480   -0.1900 0.0667 32 -2.850  0.0226 
    ##  120 - 480  -0.0576 0.0667 32 -0.863  0.7778 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

Distance to centroid for predators is greater in pods with fish and in
higher isolation treatments, but only for the last survey

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3_predators$observed_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Distance to Centroid (Observed)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n", main = "Predatory Insects")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3_predators$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3_predators$observed_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3_predators$observed_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses---Pred-and-Non-Pred_files/figure-gfm/7-1.png)<!-- -->

#### Expected Community Variability

Running ANOVA table for expected distances to group centroids, or
expected beta-diversity/community variability in each treatment.

``` r
fit_expected_SS2_SS3 <- lmer(beta_deviation_SS2_SS3_predators$expected_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_expected_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3_predators$expected_distances
    ##                                            Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                              12.353  1     <2e-16 ***
    ## isolation_SS2_SS3                          3.958  2      0.138    
    ## SS_SS2_SS3                                 1.930  1      0.165    
    ## fish_SS2_SS3:isolation_SS2_SS3             3.541  2      0.170    
    ## fish_SS2_SS3:SS_SS2_SS3                    3.975  1      0.046 *  
    ## isolation_SS2_SS3:SS_SS2_SS3               6.463  2      0.039 *  
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  6.679  2      0.035 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Patterns are similar to those observed for the observed distances to
centroid.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3_predators$expected_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Distance to Centroid (Expected)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n", main = "Predatory Insects")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3_predators$expected_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3_predators$expected_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3_predators$expected_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses---Pred-and-Non-Pred_files/figure-gfm/9-1.png)<!-- -->

#### Beta-Deviation

Running ANOVA table for the deviations of expected distances to group
centroids from observed distances.

``` r
fit_deviation_SS2_SS3 <- lmer(beta_deviation_SS2_SS3_predators$deviation_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_deviation_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3_predators$deviation_distances
    ##                                           Chisq Df Pr(>Chisq)   
    ## fish_SS2_SS3                              4.384  1      0.036 * 
    ## isolation_SS2_SS3                         3.746  2      0.154   
    ## SS_SS2_SS3                                0.912  1      0.340   
    ## fish_SS2_SS3:isolation_SS2_SS3            2.067  2      0.356   
    ## fish_SS2_SS3:SS_SS2_SS3                   9.372  1      0.002 **
    ## isolation_SS2_SS3:SS_SS2_SS3              0.703  2      0.704   
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 4.893  2      0.087 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_deviation_SS2_SS3, list(pairwise ~ fish_SS2_SS3|isolation_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of fish_SS2_SS3 | isolation_SS2_SS3`
    ## isolation_SS2_SS3 = 30:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent        1.113 0.490 14.9   -0.105     2.33
    ##  present       0.259 0.535 17.4   -1.049     1.57
    ## 
    ## isolation_SS2_SS3 = 120:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent        1.869 0.490 14.9    0.651     3.09
    ##  present       0.538 0.535 17.4   -0.770     1.85
    ## 
    ## isolation_SS2_SS3 = 480:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent        0.262 0.535 17.4   -1.046     1.57
    ##  present       0.143 0.535 17.4   -1.165     1.45
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of fish_SS2_SS3 | isolation_SS2_SS3`
    ## isolation_SS2_SS3 = 30:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    0.854 0.725 16.2 1.177   0.2562 
    ## 
    ## isolation_SS2_SS3 = 120:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    1.331 0.725 16.2 1.835   0.0849 
    ## 
    ## isolation_SS2_SS3 = 480:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    0.119 0.756 17.4 0.157   0.8769 
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger

Beta deviation seems to increase with isolation, but only in fishless
ponds.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3_predators$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Beta-Deviation", xlab = "", at = c(1,2,3,5,6,7),ylim = c(-2,6), lwd = 1.5, col = "transparent", xaxt="n", main = "Predatory Insects")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3_predators$deviation_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3_predators$deviation_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3_predators$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
box(lwd = 2.5)
```

![](Community-Variability-Analyses---Pred-and-Non-Pred_files/figure-gfm/11-1.png)<!-- -->

### Considering only Non-Predatory insects for the last two surveys.

First loading data

``` r
data(com_SS2_SS3_non_predators)
```

Computing observed and expected distances to centroid, and
beta-deviation.

``` r
beta_deviation_SS2_SS3_non_predators <- beta_deviation(com_SS2_SS3_non_predators, strata = All, times = 10000,
                                      transform = NULL, dist = "bray", fixedmar="both",
                                      shuffle = "both", method = "quasiswap", seed = 2, group = All) 
```

#### Observed Community Variability

Running ANOVA table for observed distances to group centroids, or
observed beta-diversity/community variability in each treatment.

``` r
fit_observed_SS2_SS3 <- lmer(beta_deviation_SS2_SS3_non_predators$observed_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
round(Anova(fit_observed_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3_non_predators$observed_distances
    ##                                            Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                               3.887  1      0.049 *  
    ## isolation_SS2_SS3                          3.178  2      0.204    
    ## SS_SS2_SS3                                 3.265  1      0.071 .  
    ## fish_SS2_SS3:isolation_SS2_SS3            13.724  2      0.001 ***
    ## fish_SS2_SS3:SS_SS2_SS3                    0.292  1      0.589    
    ## isolation_SS2_SS3:SS_SS2_SS3              11.718  2      0.003 ** 
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  5.014  2      0.082 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_observed_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.283 0.0580 31.8   0.1368    0.429
    ##  120                0.450 0.0580 31.8   0.3038    0.596
    ##  480                0.472 0.0580 31.8   0.3257    0.618
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.413 0.0580 31.8   0.2663    0.559
    ##  120                0.422 0.0580 31.8   0.2756    0.568
    ##  480                0.197 0.0580 31.8   0.0503    0.343
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.508 0.0580 31.8   0.3613    0.654
    ##  120                0.475 0.0580 31.8   0.3290    0.622
    ##  480                0.441 0.0681 32.0   0.2690    0.612
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.596 0.0681 32.0   0.4241    0.767
    ##  120                0.234 0.0681 32.0   0.0622    0.405
    ##  480                0.325 0.0681 32.0   0.1536    0.497
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  30 - 120   -0.1670 0.0821 31.8 -2.034  0.1436 
    ##  30 - 480   -0.1888 0.0821 31.8 -2.301  0.0820 
    ##  120 - 480  -0.0219 0.0821 31.8 -0.267  0.9909 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  30 - 120   -0.0093 0.0821 31.8 -0.113  0.9993 
    ##  30 - 480    0.2160 0.0821 31.8  2.632  0.0385 
    ##  120 - 480   0.2253 0.0821 31.8  2.745  0.0293 
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  30 - 120    0.0323 0.0821 31.8  0.393  0.9722 
    ##  30 - 480    0.0669 0.0895 31.9  0.748  0.8425 
    ##  120 - 480   0.0347 0.0895 31.9  0.387  0.9733 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE   df t.ratio p.value
    ##  30 - 120    0.3619 0.0963 32.0  3.757  0.0021 
    ##  30 - 480    0.2705 0.0963 32.0  2.808  0.0251 
    ##  120 - 480  -0.0914 0.0963 32.0 -0.949  0.7253 
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

Similar to when we analysed the whole community, it seems that the
effect of isolation is dependent on the presence or absence of fish.
When fish is absent, there is no effect of isolation. When it is
present, there is a negative effect of isolation. This effect of
isolation is stronger from low and intermediate isolation to high
isolation in the second survey, and stronger from low to intermediate
and high isolation in the last survey.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3_non_predators$observed_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Distance to Centroid (Observed)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n", main = "Non-Predators")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3_non_predators$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3_non_predators$observed_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3_non_predators$observed_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses---Pred-and-Non-Pred_files/figure-gfm/15-1.png)<!-- -->

#### Expected Community Variability

Running ANOVA table for expected distances to group centroids, or
expected beta-diversity/community variability in each treatment.

``` r
fit_expected_SS2_SS3 <- lmer(beta_deviation_SS2_SS3_non_predators$expected_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_expected_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3_non_predators$expected_distances
    ##                                            Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                               0.821  1      0.365    
    ## isolation_SS2_SS3                         10.389  2      0.006 ** 
    ## SS_SS2_SS3                                 4.924  1      0.026 *  
    ## fish_SS2_SS3:isolation_SS2_SS3             4.296  2      0.117    
    ## fish_SS2_SS3:SS_SS2_SS3                    0.874  1      0.350    
    ## isolation_SS2_SS3:SS_SS2_SS3              15.799  2     <2e-16 ***
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  1.834  2      0.400    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.324 0.0405 32    0.222    0.426
    ##  120                0.429 0.0405 32    0.327    0.531
    ##  480                0.255 0.0405 32    0.153    0.356
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.564 0.0441 32    0.453    0.675
    ##  120                0.327 0.0441 32    0.216    0.439
    ##  480                0.343 0.0475 32    0.223    0.463
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   -0.1052 0.0572 32 -1.839  0.2091 
    ##  30 - 480    0.0692 0.0572 32  1.210  0.5529 
    ##  120 - 480   0.1744 0.0572 32  3.048  0.0137 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120    0.2366 0.0624 32  3.790  0.0019 
    ##  30 - 480    0.2209 0.0649 32  3.405  0.0054 
    ##  120 - 480  -0.0157 0.0649 32 -0.242  0.9931 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

There is a decrease in the expected distance to centroid with isolation,
and this decrease is stronger in for the last survey\!

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3_non_predators$expected_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Distance to Centroid (Expected)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n", main = "Non-Predators")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3_non_predators$expected_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3_non_predators$expected_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3_non_predators$expected_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses---Pred-and-Non-Pred_files/figure-gfm/17-1.png)<!-- -->

#### Beta-Deviation

Running ANOVA table for the deviations of expected distances to group
centroids from observed distances.

``` r
fit_deviation_SS2_SS3 <- lmer(beta_deviation_SS2_SS3_non_predators$deviation_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_deviation_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3_non_predators$deviation_distances
    ##                                            Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                               4.977  1      0.026 *  
    ## isolation_SS2_SS3                          6.650  2      0.036 *  
    ## SS_SS2_SS3                                 0.074  1      0.785    
    ## fish_SS2_SS3:isolation_SS2_SS3             7.499  2      0.024 *  
    ## fish_SS2_SS3:SS_SS2_SS3                    0.197  1      0.657    
    ## isolation_SS2_SS3:SS_SS2_SS3               3.644  2      0.162    
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 13.242  2      0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_deviation_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                -0.202 0.667 32   -1.882     1.48
    ##  120                0.580 0.667 32   -1.100     2.26
    ##  480                3.216 0.667 32    1.536     4.90
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                 0.856 0.667 32   -0.824     2.54
    ##  120               -0.253 0.667 32   -1.933     1.43
    ##  480               -0.230 0.667 32   -1.910     1.45
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                -0.260 0.667 32   -1.940     1.42
    ##  120                2.812 0.667 32    1.132     4.49
    ##  480                0.517 0.784 32   -1.457     2.49
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                -0.143 0.784 32   -2.117     1.83
    ##  120               -0.309 0.784 32   -2.283     1.67
    ##  480                1.684 0.784 32   -0.290     3.66
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  contrast  estimate    SE df t.ratio p.value
    ##  30 - 120   -0.7822 0.943 32 -0.829  0.7977 
    ##  30 - 480   -3.4182 0.943 32 -3.625  0.0030 
    ##  120 - 480  -2.6359 0.943 32 -2.795  0.0259 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  contrast  estimate    SE df t.ratio p.value
    ##  30 - 120    1.1086 0.943 32  1.176  0.5755 
    ##  30 - 480    1.0863 0.943 32  1.152  0.5913 
    ##  120 - 480  -0.0222 0.943 32 -0.024  1.0000 
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  contrast  estimate    SE df t.ratio p.value
    ##  30 - 120   -3.0721 0.943 32 -3.258  0.0080 
    ##  30 - 480   -0.7769 1.029 32 -0.755  0.8388 
    ##  120 - 480   2.2953 1.029 32  2.231  0.0953 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  contrast  estimate    SE df t.ratio p.value
    ##  30 - 120    0.1661 1.108 32  0.150  0.9983 
    ##  30 - 480   -1.8266 1.108 32 -1.648  0.2928 
    ##  120 - 480  -1.9927 1.108 32 -1.798  0.2253 
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

Similarly to patterns for the whole community, there is an increase in
beta deviation with isolation, but only in fishless ponds. This effect
seems stronger in the second survey.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3_non_predators$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Beta-Deviation", xlab = "", at = c(1,2,3,5,6,7),ylim = c(-3,8), lwd = 1.5, col = "transparent", xaxt="n", main = "Non-Predators")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3_non_predators$deviation_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3_non_predators$deviation_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3_non_predators$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
box(lwd = 2.5)
```

![](Community-Variability-Analyses---Pred-and-Non-Pred_files/figure-gfm/19-1.png)<!-- -->
