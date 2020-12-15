Community Variability - Varying Gamma
================
Rodolfo Pelinson
20/10/2020

``` r
library(Pelinson.et.al.2020B)
```

``` r
library(lme4)
library(car)
library(emmeans)
library(vegan)
```

## Community Variability Varying Gamma

### Whole community for the last two surveys.

First loading data

``` r
data(com_SS2_SS3, All, fish_SS2_SS3, isolation_SS2_SS3, SS_SS2_SS3, ID_SS2_SS3, fish_isolation_SS2_SS3)
```

    ## Warning in data(com_SS2_SS3, All, fish_SS2_SS3, isolation_SS2_SS3, SS_SS2_SS3, :
    ## data set 'fish_isolation_SS2_SS3' not found

Computing observed and expected distances to centroid, and
beta-deviation.

``` r
beta_deviation_SS2_SS3 <- beta_deviation(com_SS2_SS3, strata = NULL, times = 10000,
                                      transform = NULL, dist = "bray", fixedmar="both",
                                      shuffle = "both", method = "quasiswap", seed = 2, group = All) 
```

#### Expected Community Variability

Running ANOVA table for expected distances to group centroids, or
expected beta-diversity/community variability in each treatment.

``` r
fit_expected_SS2_SS3 <- lmer(beta_deviation_SS2_SS3$expected_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_expected_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3$expected_distances
    ##                                           Chisq Df Pr(>Chisq)   
    ## fish_SS2_SS3                              0.222  1      0.637   
    ## isolation_SS2_SS3                         1.895  2      0.388   
    ## SS_SS2_SS3                                0.708  1      0.400   
    ## fish_SS2_SS3:isolation_SS2_SS3            6.263  2      0.044 * 
    ## fish_SS2_SS3:SS_SS2_SS3                   0.614  1      0.433   
    ## isolation_SS2_SS3:SS_SS2_SS3              9.222  2      0.010 **
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 8.756  2      0.013 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.336 0.0322 14.9    0.249    0.422
    ##  120                0.377 0.0322 14.9    0.290    0.463
    ##  480                0.377 0.0351 17.4    0.285    0.470
    ## 
    ## fish_SS2_SS3 = present:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.444 0.0351 17.4    0.351    0.536
    ##  120                0.383 0.0351 17.4    0.290    0.476
    ##  480                0.301 0.0351 17.4    0.208    0.394
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  contrast   estimate     SE   df t.ratio p.value
    ##  30 - 120  -0.040761 0.0455 14.9 -0.895  0.7671 
    ##  30 - 480  -0.041628 0.0476 16.2 -0.874  0.7785 
    ##  120 - 480 -0.000867 0.0476 16.2 -0.018  1.0000 
    ## 
    ## fish_SS2_SS3 = present:
    ##  contrast   estimate     SE   df t.ratio p.value
    ##  30 - 120   0.060696 0.0497 17.4  1.222  0.5576 
    ##  30 - 480   0.142622 0.0497 17.4  2.871  0.0309 
    ##  120 - 480  0.081926 0.0497 17.4  1.649  0.3115 
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.316 0.0322 32    0.235    0.397
    ##  120                0.410 0.0322 32    0.329    0.491
    ##  480                0.346 0.0322 32    0.265    0.427
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.463 0.0351 32    0.375    0.552
    ##  120                0.350 0.0351 32    0.261    0.438
    ##  480                0.332 0.0378 32    0.237    0.428
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   -0.0937 0.0455 32 -2.059  0.1364 
    ##  30 - 480   -0.0298 0.0455 32 -0.655  0.8873 
    ##  120 - 480   0.0639 0.0455 32  1.404  0.4283 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120    0.1137 0.0497 32  2.288  0.0841 
    ##  30 - 480    0.1308 0.0516 32  2.534  0.0483 
    ##  120 - 480   0.0172 0.0516 32  0.332  0.9828 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.329 0.0455 32    0.215    0.444
    ##  120                0.381 0.0455 32    0.266    0.495
    ##  480                0.374 0.0455 32    0.260    0.489
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.303 0.0455 32    0.188    0.418
    ##  120                0.439 0.0455 32    0.324    0.554
    ##  480                0.318 0.0455 32    0.203    0.432
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.342 0.0455 32    0.228    0.457
    ##  120                0.372 0.0455 32    0.258    0.487
    ##  480                0.380 0.0535 32    0.246    0.515
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.584 0.0535 32    0.449    0.719
    ##  120                0.327 0.0535 32    0.192    0.461
    ##  480                0.284 0.0535 32    0.149    0.419
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120  -0.05137 0.0644 32 -0.798  0.8155 
    ##  30 - 480  -0.04512 0.0644 32 -0.701  0.8661 
    ##  120 - 480  0.00625 0.0644 32  0.097  0.9995 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120  -0.13608 0.0644 32 -2.114  0.1219 
    ##  30 - 480  -0.01454 0.0644 32 -0.226  0.9944 
    ##  120 - 480  0.12155 0.0644 32  1.888  0.1908 
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120  -0.03015 0.0644 32 -0.468  0.9544 
    ##  30 - 480  -0.03813 0.0702 32 -0.543  0.9316 
    ##  120 - 480 -0.00799 0.0702 32 -0.114  0.9993 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   0.25747 0.0757 32  3.403  0.0054 
    ##  30 - 480   0.29978 0.0757 32  3.963  0.0012 
    ##  120 - 480  0.04231 0.0757 32  0.559  0.9259 
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

Here the increase in expected distance to centroid increases with
isolation only in ponds with fish, in the last survey.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3$expected_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Distance to Centroid (Expected)", xlab = "", at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3$expected_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3$expected_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3$expected_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability---Varying-Gamma_files/figure-gfm/6-1.png)<!-- -->

#### Beta-Deviation

Running ANOVA table for the deviations of expected distances to group
centroids from observed distances.

``` r
fit_deviation_SS2_SS3 <- lmer(beta_deviation_SS2_SS3$deviation_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_deviation_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3$deviation_distances
    ##                                           Chisq Df Pr(>Chisq)  
    ## fish_SS2_SS3                              4.215  1      0.040 *
    ## isolation_SS2_SS3                         0.247  2      0.884  
    ## SS_SS2_SS3                                3.392  1      0.066 .
    ## fish_SS2_SS3:isolation_SS2_SS3            8.105  2      0.017 *
    ## fish_SS2_SS3:SS_SS2_SS3                   2.620  1      0.106  
    ## isolation_SS2_SS3:SS_SS2_SS3              4.941  2      0.085 .
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 1.813  2      0.404  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_deviation_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  isolation_SS2_SS3   emmean    SE   df lower.CL upper.CL
    ##  30                -0.00175 0.568 15.1  -1.5253    1.522
    ##  120                1.71806 0.568 15.1   0.1945    3.242
    ##  480                1.53239 0.615 17.6  -0.0898    3.155
    ## 
    ## fish_SS2_SS3 = present:
    ##  isolation_SS2_SS3   emmean    SE   df lower.CL upper.CL
    ##  30                 0.77421 0.615 17.6  -0.8480    2.396
    ##  120               -0.72849 0.615 17.6  -2.3507    0.894
    ##  480                0.03071 0.615 17.6  -1.5915    1.653
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  contrast  estimate    SE   df t.ratio p.value
    ##  30 - 120    -1.720 0.803 15.1 -2.141  0.1398 
    ##  30 - 480    -1.534 0.837 16.4 -1.832  0.2342 
    ##  120 - 480    0.186 0.837 16.4  0.222  0.9948 
    ## 
    ## fish_SS2_SS3 = present:
    ##  contrast  estimate    SE   df t.ratio p.value
    ##  30 - 120     1.503 0.870 17.6  1.727  0.2750 
    ##  30 - 480     0.743 0.870 17.6  0.855  0.7886 
    ##  120 - 480   -0.759 0.870 17.6 -0.873  0.7781 
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

No significant patterns.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Beta-Deviation", xlab = "", at = c(1,2,3,5,6,7),ylim = c(-2,10), lwd = 1.5, col = "transparent", xaxt="n")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3$deviation_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3$deviation_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
box(lwd = 2.5)
```

![](Community-Variability---Varying-Gamma_files/figure-gfm/8-1.png)<!-- -->
