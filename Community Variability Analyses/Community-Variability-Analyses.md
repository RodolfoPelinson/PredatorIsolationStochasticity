Community Variability Analyses
================
Rodolfo Pelinson
20/10/2020

This is the community variability analyses presented in the main paper.

These are the packages you will need to run this code:

``` r
library(lme4) # Version 1.1-23
library(car) # Version 3.0-7
library(emmeans) # Version 1.4.8
library(vegan) # Version 2.5-6
```

## Community Variability

### Whole community for the last two surveys.

First loading data

``` r
data(com_SS2_SS3, All, fish_SS2_SS3, isolation_SS2_SS3, SS_SS2_SS3, ID_SS2_SS3,
     fish_isolation_SS2_SS3)
```

    ## Warning in data(com_SS2_SS3, All, fish_SS2_SS3, isolation_SS2_SS3, SS_SS2_SS3, :
    ## data set 'fish_isolation_SS2_SS3' not found

Computing observed and expected distances to centroid, and
beta-deviation.

``` r
beta_deviation_SS2_SS3 <- beta_deviation(com_SS2_SS3, strata = All, times = 10000,
                                      transform = NULL, dist = "bray", fixedmar="both",
                                      shuffle = "both", method = "quasiswap", seed = 2,
                                      group = All) 
```

Looking at residual plots for observed, expected distances to centroids
and deviations.

``` r
fit_observed_SS2_SS3_G <- lmer(beta_deviation_SS2_SS3$observed_distances~fish_SS2_SS3*
                                 isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
plot(fit_observed_SS2_SS3_G)
```

![](Community-Variability-Analyses_files/figure-gfm/checking_distribution-1.png)<!-- -->

``` r
qqnorm(resid(fit_observed_SS2_SS3_G, type = "pearson"))
qqline(resid(fit_observed_SS2_SS3_G, type = "pearson"))
```

![](Community-Variability-Analyses_files/figure-gfm/checking_distribution-2.png)<!-- -->

``` r
fit_expected_SS2_SS3_G <- lmer(beta_deviation_SS2_SS3$expected_distances~fish_SS2_SS3*
                                 isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), REML = F)
```

    ## boundary (singular) fit: see ?isSingular

``` r
plot(fit_expected_SS2_SS3_G)
```

![](Community-Variability-Analyses_files/figure-gfm/checking_distribution-3.png)<!-- -->

``` r
qqnorm(resid(fit_expected_SS2_SS3_G, type = "pearson"))
qqline(resid(fit_expected_SS2_SS3_G, type = "pearson"))
```

![](Community-Variability-Analyses_files/figure-gfm/checking_distribution-4.png)<!-- -->

``` r
fit_deviation_SS2_SS3_G <- lmer(beta_deviation_SS2_SS3$deviation_distances~fish_SS2_SS3*
                                  isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), REML = F)
plot(fit_deviation_SS2_SS3_G)
```

![](Community-Variability-Analyses_files/figure-gfm/checking_distribution-5.png)<!-- -->

``` r
qqnorm(resid(fit_deviation_SS2_SS3_G, type = "pearson"))
qqline(resid(fit_deviation_SS2_SS3_G, type = "pearson"))
```

![](Community-Variability-Analyses_files/figure-gfm/checking_distribution-6.png)<!-- -->
Residual plots are not perfect, but they also does not seem too bad.

#### Observed Community Variability

Running ANOVA table for observed distances to group centroids, or
observed beta-diversity/community variability in each treatment.

``` r
fit_observed_SS2_SS3 <- lmer(beta_deviation_SS2_SS3$observed_distances~fish_SS2_SS3*
                               isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
round(Anova(fit_observed_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3$observed_distances
    ##                                            Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                               1.442  1      0.230    
    ## isolation_SS2_SS3                          2.071  2      0.355    
    ## SS_SS2_SS3                                 2.136  1      0.144    
    ## fish_SS2_SS3:isolation_SS2_SS3            17.217  2     <2e-16 ***
    ## fish_SS2_SS3:SS_SS2_SS3                    0.086  1      0.769    
    ## isolation_SS2_SS3:SS_SS2_SS3               7.317  2      0.026 *  
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  5.463  2      0.065 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_observed_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3),
        adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.322 0.0395 15.3    0.238    0.406
    ##  120                0.450 0.0395 15.3    0.366    0.534
    ##  480                0.425 0.0424 17.9    0.336    0.514
    ## 
    ## fish_SS2_SS3 = present:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.479 0.0424 17.9    0.390    0.568
    ##  120                0.329 0.0424 17.9    0.240    0.418
    ##  480                0.269 0.0424 17.9    0.179    0.358
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  2         estimate     SE   df t.ratio p.value
    ##  30 - 120   -0.1280 0.0559 15.3  -2.291  0.1056
    ##  30 - 480   -0.1029 0.0580 16.6  -1.776  0.2565
    ##  120 - 480   0.0250 0.0580 16.6   0.432  0.9645
    ## 
    ## fish_SS2_SS3 = present:
    ##  2         estimate     SE   df t.ratio p.value
    ##  30 - 120    0.1499 0.0600 17.9   2.499  0.0659
    ##  30 - 480    0.2105 0.0600 17.9   3.509  0.0076
    ##  120 - 480   0.0606 0.0600 17.9   1.010  0.6938
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_observed_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3),
        adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.340 0.0362 31.0    0.266    0.413
    ##  120                0.424 0.0362 31.0    0.350    0.498
    ##  480                0.311 0.0362 31.0    0.237    0.385
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.461 0.0394 31.5    0.381    0.542
    ##  120                0.355 0.0394 31.5    0.275    0.435
    ##  480                0.382 0.0423 31.8    0.296    0.468
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  2         estimate     SE   df t.ratio p.value
    ##  30 - 120   -0.0843 0.0512 31.0  -1.646  0.2948
    ##  30 - 480    0.0284 0.0512 31.0   0.554  0.9277
    ##  120 - 480   0.1127 0.0512 31.0   2.200  0.1024
    ## 
    ## SS_SS2_SS3 = 3:
    ##  2         estimate     SE   df t.ratio p.value
    ##  30 - 120    0.1062 0.0557 31.5   1.908  0.1841
    ##  30 - 480    0.0792 0.0578 31.7   1.370  0.4490
    ##  120 - 480  -0.0270 0.0578 31.7  -0.468  0.9546
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_observed_SS2_SS3, list(pairwise ~ SS_SS2_SS3|isolation_SS2_SS3),
        adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of SS_SS2_SS3 | isolation_SS2_SS3`
    ## isolation_SS2_SS3 = 30:
    ##  SS_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  2           0.340 0.0362 31.0    0.266    0.413
    ##  3           0.461 0.0394 31.5    0.381    0.542
    ## 
    ## isolation_SS2_SS3 = 120:
    ##  SS_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  2           0.424 0.0362 31.0    0.350    0.498
    ##  3           0.355 0.0394 31.5    0.275    0.435
    ## 
    ## isolation_SS2_SS3 = 480:
    ##  SS_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  2           0.311 0.0362 31.0    0.237    0.385
    ##  3           0.382 0.0423 31.8    0.296    0.468
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of SS_SS2_SS3 | isolation_SS2_SS3`
    ## isolation_SS2_SS3 = 30:
    ##  2     estimate     SE   df t.ratio p.value
    ##  2 - 3  -0.1218 0.0486 15.8  -2.505  0.0236
    ## 
    ## isolation_SS2_SS3 = 120:
    ##  2     estimate     SE   df t.ratio p.value
    ##  2 - 3   0.0687 0.0486 15.8   1.414  0.1769
    ## 
    ## isolation_SS2_SS3 = 480:
    ##  2     estimate     SE   df t.ratio p.value
    ##  2 - 3  -0.0710 0.0510 16.8  -1.391  0.1823
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger

It seems that the effect of isolation is dependent on the presence or
absence of fish. When fish is absent, there is no effect of isolation.
When it is present, there is a negative effect of isolation.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3$observed_distances~isolation_SS2_SS3*fish_SS2_SS3,
        outline = F, ylab = "Distance to Centroid (Observed)", xlab = "",
        at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n")
mylevels <- levels(All)
levelProportions <- summary(All)/length(beta_deviation_SS2_SS3$observed_distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- beta_deviation_SS2_SS3$observed_distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(beta_deviation_SS2_SS3$observed_distances~isolation_SS2_SS3*fish_SS2_SS3,
        add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7),
        lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8,
     at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses_files/figure-gfm/plot_observed-1.png)<!-- -->

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
    ##                                            Chisq Df Pr(>Chisq)   
    ## fish_SS2_SS3                               0.011  1      0.916   
    ## isolation_SS2_SS3                          7.063  2      0.029 * 
    ## SS_SS2_SS3                                 3.639  1      0.056 . 
    ## fish_SS2_SS3:isolation_SS2_SS3             8.425  2      0.015 * 
    ## fish_SS2_SS3:SS_SS2_SS3                    0.170  1      0.680   
    ## isolation_SS2_SS3:SS_SS2_SS3              12.700  2      0.002 **
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  6.840  2      0.033 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3),
        adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3`
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.374 0.0255 16.2    0.320    0.428
    ##  120                0.347 0.0255 16.2    0.293    0.401
    ##  480                0.278 0.0266 17.4    0.222    0.334
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3, SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3`
    ##  1         estimate     SE   df t.ratio p.value
    ##  30 - 120    0.0272 0.0361 16.2   0.755  0.8434
    ##  30 - 480    0.0962 0.0369 16.8   2.611  0.0541
    ##  120 - 480   0.0690 0.0369 16.8   1.872  0.2181
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3, SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3),
        adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.316 0.0345 14.9    0.243    0.390
    ##  120                0.353 0.0345 14.9    0.280    0.427
    ##  480                0.329 0.0376 17.4    0.250    0.409
    ## 
    ## fish_SS2_SS3 = present:
    ##  isolation_SS2_SS3 emmean     SE   df lower.CL upper.CL
    ##  30                 0.432 0.0376 17.4    0.352    0.511
    ##  120                0.340 0.0376 17.4    0.261    0.420
    ##  480                0.226 0.0376 17.4    0.147    0.305
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  2         estimate     SE   df t.ratio p.value
    ##  30 - 120   -0.0368 0.0488 14.9  -0.755  0.8444
    ##  30 - 480   -0.0130 0.0510 16.2  -0.256  0.9922
    ##  120 - 480   0.0238 0.0510 16.2   0.465  0.9563
    ## 
    ## fish_SS2_SS3 = present:
    ##  2         estimate     SE   df t.ratio p.value
    ##  30 - 120    0.0913 0.0532 17.4   1.716  0.2805
    ##  30 - 480    0.2055 0.0532 17.4   3.863  0.0036
    ##  120 - 480   0.1142 0.0532 17.4   2.147  0.1322
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3),
        adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.287 0.0345 32    0.217    0.357
    ##  120                0.392 0.0345 32    0.322    0.462
    ##  480                0.235 0.0345 32    0.165    0.305
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.461 0.0376 32    0.385    0.538
    ##  120                0.302 0.0376 32    0.225    0.378
    ##  480                0.321 0.0405 32    0.238    0.403
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  2         estimate     SE df t.ratio p.value
    ##  30 - 120   -0.1051 0.0488 32  -2.155  0.1120
    ##  30 - 480    0.0520 0.0488 32   1.067  0.6482
    ##  120 - 480   0.1571 0.0488 32   3.221  0.0088
    ## 
    ## SS_SS2_SS3 = 3:
    ##  2         estimate     SE df t.ratio p.value
    ##  30 - 120    0.1595 0.0532 32   2.999  0.0155
    ##  30 - 480    0.1405 0.0553 32   2.541  0.0476
    ##  120 - 480  -0.0191 0.0553 32  -0.345  0.9808
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3|SS_SS2_SS3),
        adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.271 0.0488 32   0.1717    0.370
    ##  120                0.346 0.0488 32   0.2470    0.446
    ##  480                0.278 0.0488 32   0.1786    0.377
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.303 0.0488 32   0.2032    0.402
    ##  120                0.437 0.0488 32   0.3380    0.537
    ##  480                0.192 0.0488 32   0.0924    0.291
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.362 0.0488 32   0.2623    0.461
    ##  120                0.360 0.0488 32   0.2606    0.459
    ##  480                0.381 0.0573 32   0.2641    0.498
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.561 0.0573 32   0.4440    0.677
    ##  120                0.243 0.0573 32   0.1267    0.360
    ##  480                0.261 0.0573 32   0.1439    0.377
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  3         estimate     SE df t.ratio p.value
    ##  30 - 120  -0.07533 0.0690 32  -1.092  0.6311
    ##  30 - 480  -0.00687 0.0690 32  -0.100  0.9995
    ##  120 - 480  0.06846 0.0690 32   0.993  0.6968
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  3         estimate     SE df t.ratio p.value
    ##  30 - 120  -0.13478 0.0690 32  -1.955  0.1679
    ##  30 - 480   0.11089 0.0690 32   1.608  0.3130
    ##  120 - 480  0.24566 0.0690 32   3.563  0.0035
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  3         estimate     SE df t.ratio p.value
    ##  30 - 120   0.00173 0.0690 32   0.025  1.0000
    ##  30 - 480  -0.01923 0.0752 32  -0.256  0.9920
    ##  120 - 480 -0.02096 0.0752 32  -0.279  0.9897
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  3         estimate     SE df t.ratio p.value
    ##  30 - 120   0.31735 0.0810 32   3.916  0.0013
    ##  30 - 480   0.30015 0.0810 32   3.704  0.0024
    ##  120 - 480 -0.01720 0.0810 32  -0.212  0.9954
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

Patterns are similar to those observed for the observed distances to
centroid.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3$expected_distances~isolation_SS2_SS3*fish_SS2_SS3,
        outline = F, ylab = "Distance to Centroid (Expected)", xlab = "",
        at = c(1,2,3,5,6,7),ylim = c(0,1), lwd = 1.5, col = "transparent", xaxt="n")
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
boxplot(beta_deviation_SS2_SS3$expected_distances~isolation_SS2_SS3*fish_SS2_SS3,
        add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7),
        lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8,
     at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses_files/figure-gfm/plot_expected-1.png)<!-- -->

#### Beta-Deviation

Running ANOVA table for the deviations of expected distances to group
centroids from observed distances.

``` r
fit_deviation_SS2_SS3 <- lmer(beta_deviation_SS2_SS3$deviation_distances~fish_SS2_SS3*
                                isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3),
                              control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_deviation_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3$deviation_distances
    ##                                           Chisq Df Pr(>Chisq)   
    ## fish_SS2_SS3                              4.129  1      0.042 * 
    ## isolation_SS2_SS3                         4.739  2      0.094 . 
    ## SS_SS2_SS3                                1.048  1      0.306   
    ## fish_SS2_SS3:isolation_SS2_SS3            9.345  2      0.009 **
    ## fish_SS2_SS3:SS_SS2_SS3                   0.096  1      0.757   
    ## isolation_SS2_SS3:SS_SS2_SS3              2.445  2      0.295   
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 2.595  2      0.273   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_deviation_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3),
        adjust = "sidak")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  isolation_SS2_SS3  emmean    SE   df lower.CL upper.CL
    ##  30                 0.0434 0.556 14.9  -1.1410     1.23
    ##  120                2.5286 0.556 14.9   1.3442     3.71
    ##  480                2.3009 0.606 17.4   1.0249     3.58
    ## 
    ## fish_SS2_SS3 = present:
    ##  isolation_SS2_SS3  emmean    SE   df lower.CL upper.CL
    ##  30                 0.9768 0.606 17.4  -0.2991     2.25
    ##  120               -0.0801 0.606 17.4  -1.3561     1.20
    ##  480                1.1787 0.606 17.4  -0.0973     2.45
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  2         estimate    SE   df t.ratio p.value
    ##  30 - 120    -2.485 0.786 14.9  -3.163  0.0192
    ##  30 - 480    -2.257 0.822 16.2  -2.746  0.0420
    ##  120 - 480    0.228 0.822 16.2   0.277  0.9901
    ## 
    ## fish_SS2_SS3 = present:
    ##  2         estimate    SE   df t.ratio p.value
    ##  30 - 120     1.057 0.857 17.4   1.234  0.5502
    ##  30 - 480    -0.202 0.857 17.4  -0.236  0.9938
    ##  120 - 480   -1.259 0.857 17.4  -1.469  0.4066
    ## 
    ## Results are averaged over the levels of: SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

Beta deviation seems to increase with isolation, but only in fishless
ponds.

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3,
        outline = F, ylab = "Beta-Deviation", xlab = "", at = c(1,2,3,5,6,7),
        ylim = c(-2,10), lwd = 1.5, col = "transparent", xaxt="n")
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
boxplot(beta_deviation_SS2_SS3$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3,
        add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7),
        lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8,
     at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
abline(h = 0, lty = 2, lwd = 2, col = "grey50")
box(lwd = 2.5)
```

![](Community-Variability-Analyses_files/figure-gfm/plot_deviation-1.png)<!-- -->

#### Observed Environmental Variability

\#\#\#\#Multivariate effect of Environmental Variability

``` r
env_data_SS2_SS3_st <- decostand(env_data_SS2_SS3, method = "stand", na.rm = TRUE)

env_data_SS2_SS3_dist <- vegdist(env_data_SS2_SS3_st, method = "euclidean", na.rm = TRUE)

betadisper_SS2_SS3 <- betadisper(env_data_SS2_SS3_dist, group = All)


fit_observed_SS2_SS3 <- lmer(betadisper_SS2_SS3$distances~fish_SS2_SS3*
                               isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))

round(Anova(fit_observed_SS2_SS3, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: betadisper_SS2_SS3$distances
    ##                                           Chisq Df Pr(>Chisq)
    ## fish_SS2_SS3                              1.823  1      0.177
    ## isolation_SS2_SS3                         0.091  2      0.956
    ## SS_SS2_SS3                                1.812  1      0.178
    ## fish_SS2_SS3:isolation_SS2_SS3            0.306  2      0.858
    ## fish_SS2_SS3:SS_SS2_SS3                   0.618  1      0.432
    ## isolation_SS2_SS3:SS_SS2_SS3              2.555  2      0.279
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 0.203  2      0.904

``` r
fit_observed_SS2_SS3_env <- lmer(beta_deviation_SS2_SS3$observed_distances ~ betadisper_SS2_SS3$distances + (1|ID_SS2_SS3))
round(Anova(fit_observed_SS2_SS3_env, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3$observed_distances
    ##                              Chisq Df Pr(>Chisq)
    ## betadisper_SS2_SS3$distances 1.621  1      0.203

``` r
fit_deviation_SS2_SS3_env <- lmer(beta_deviation_SS2_SS3$deviation_distances ~ betadisper_SS2_SS3$distances + (1|ID_SS2_SS3))
round(Anova(fit_observed_SS2_SS3_env, test.statistic = "Chisq"),3)
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: beta_deviation_SS2_SS3$observed_distances
    ##                              Chisq Df Pr(>Chisq)
    ## betadisper_SS2_SS3$distances 1.621  1      0.203

``` r
boxplot(betadisper_SS2_SS3$distances~isolation_SS2_SS3*fish_SS2_SS3,
        outline = F, ylab = "Distance to Centroid (Environmental Variability)", xlab = "",
        at = c(1,2,3,5,6,7),ylim = c(0,5), lwd = 1.5, col = "transparent", xaxt="n")
mylevels <- levels(All)
levelProportions <- summary(All)/length(betadisper_SS2_SS3$distances)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- betadisper_SS2_SS3$distances[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(betadisper_SS2_SS3$distances~isolation_SS2_SS3*fish_SS2_SS3,
        add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7),
        lwd = 1.5, xaxt="n")
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8,
     at =c(1,2,3,5,6,7))
axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Community-Variability-Analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
par(mfrow = c(1,2))
plot(beta_deviation_SS2_SS3$observed_distances ~ betadisper_SS2_SS3$distances, pch = 16, ylab = "Distance to Centroid (Observed Comm. Variability)", xlab = "Distance to Centroid (Environmental Variability)")

plot(beta_deviation_SS2_SS3$deviation_distances ~ betadisper_SS2_SS3$distances, pch = 16, ylab = "Beta Deviation", xlab = "Distance to Centroid (Environmental Variability)")
abline(h = 0, lty = 2, col = "grey50", lwd = 2)
```

![](Community-Variability-Analyses_files/figure-gfm/plot_env-1.png)<!-- -->

\#\#\#\#Univariate effect of Environmental Variability

First, we have to built a matrix with all univariate distances, that is,
how much each observation of each variable differs from its treatment
mean.

``` r
distances_env_uni <- data.frame(matrix(NA, ncol = ncol(env_data_SS2_SS3_st), nrow = nrow(env_data_SS2_SS3_st)))
for(i in 1:ncol(env_data_SS2_SS3_st)){
  if(anyNA(env_data_SS2_SS3_st[,i])){
     na_position <- match(NA, env_data_SS2_SS3_st[,i]) 
     distances_env_uni_dist <- vegdist(env_data_SS2_SS3_st[-na_position,i], method = "euclidean", na.rm = FALSE)
     betadisper_env_uni <- betadisper(distances_env_uni_dist, group = All[-na_position])
     distances <- betadisper_env_uni$distances
     distances <- append(distances, NA, after=na_position)
     distances_env_uni[,i] <- distances
  }else{
    distances_env_uni_dist <- vegdist(env_data_SS2_SS3_st[,i], method = "euclidean", na.rm = FALSE)
    betadisper_env_uni <- betadisper(distances_env_uni_dist, group = All)
    distances_env_uni[,i] <- betadisper_env_uni$distances
  }
}
colnames(distances_env_uni) <- colnames(env_data_SS2_SS3_st)
```

Now we can run ANOVAS for each environmental variable

``` r
uni_anova_env <- data.frame(matrix(NA, ncol = ncol(distances_env_uni), nrow = 7))
for(i in 1:ncol(distances_env_uni)){
  fit <- lmer(distances_env_uni[,i]~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
  uni_anova_env[,i] <- round(Anova(fit, test.statistic = "Chisq"),3)$`Pr(>Chisq)`
}
```

    ## boundary (singular) fit: see ?isSingular
    ## boundary (singular) fit: see ?isSingular
    ## boundary (singular) fit: see ?isSingular

``` r
rownames(uni_anova_env) <- c("fish", "isolation", "survey", "fish:isolation", "fish:survey", "isolation:survey","fish:isolation:survey")
colnames(uni_anova_env) <- colnames(distances_env_uni)
uni_anova_env
```

    ##                       depth conductivity    ph    do turbidity clorophylla_a
    ## fish                  0.723        0.393 0.129 0.221     0.054         0.915
    ## isolation             0.389        0.381 0.125 0.447     0.262         0.267
    ## survey                0.179        0.313 0.020 0.609     0.303         0.075
    ## fish:isolation        0.857        0.434 0.011 0.978     0.800         0.875
    ## fish:survey           0.484        0.333 0.523 0.180     0.242         0.552
    ## isolation:survey      0.349        0.426 0.058 0.923     0.678         0.280
    ## fish:isolation:survey 0.615        0.492 0.519 0.040     0.787         0.659
    ##                       phycocyanin canopy
    ## fish                        0.235  0.435
    ## isolation                   0.265  0.206
    ## survey                      0.134  0.458
    ## fish:isolation              0.940  0.060
    ## fish:survey                 0.974  0.831
    ## isolation:survey            0.518  0.935
    ## fish:isolation:survey       0.845  0.681

Because we are blindly looking for an effect without previous hypothesis
for specific variables, it is important to correct p values of each main
and interactive effect for multiple comparisons.

``` r
uni_anova_env_adjusted_p <- uni_anova_env
for(i in 1:nrow(uni_anova_env)){
  uni_anova_env_adjusted_p[i,] <- p.adjust(uni_anova_env_adjusted_p[i,], method = "fdr") 
}
uni_anova_env_adjusted_p
```

    ##                           depth conductivity     ph    do turbidity
    ## fish                  0.8262857    0.5800000 0.4700 0.470 0.4320000
    ## isolation             0.4445714    0.4445714 0.4272 0.447 0.4272000
    ## survey                0.3580000    0.4173333 0.1600 0.609 0.4173333
    ## fish:isolation        0.9780000    0.9780000 0.0880 0.978 0.9780000
    ## fish:survey           0.7360000    0.7360000 0.7360 0.736 0.7360000
    ## isolation:survey      0.8288000    0.8288000 0.4640 0.935 0.9040000
    ## fish:isolation:survey 0.8450000    0.8450000 0.8450 0.320 0.8450000
    ##                       clorophylla_a phycocyanin    canopy
    ## fish                         0.9150   0.4700000 0.5800000
    ## isolation                    0.4272   0.4272000 0.4272000
    ## survey                       0.3000   0.3573333 0.5234286
    ## fish:isolation               0.9780   0.9780000 0.2400000
    ## fish:survey                  0.7360   0.9740000 0.9497143
    ## isolation:survey             0.8288   0.8288000 0.9350000
    ## fish:isolation:survey        0.8450   0.8450000 0.8450000

Plotting it.

``` r
par(mfrow = c(4,2))
for(j in 1:ncol(distances_env_uni)){
  boxplot(distances_env_uni[,j]~isolation_SS2_SS3*fish_SS2_SS3,
          outline = F, ylab = "Variability", xlab = "",
          at = c(1,2,3,5,6,7), ylim = c(0,max(na.omit(distances_env_uni[,j]))*1.1), lwd = 1, col = "transparent", xaxt="n", main = paste(colnames(distances_env_uni)[j]))
  mylevels <- levels(All)
  levelProportions <- summary(All)/length(distances_env_uni[,j])
  col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
  bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
  pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
  for(i in 1:length(mylevels)){
    
    x<- c(1,2,3,5,6,7,1,2,3,5,6,7)[i]
    thislevel <- mylevels[i]
    thisvalues <- distances_env_uni[,j][All==thislevel]
    
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
    points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
    
  }
  boxplot(distances_env_uni[,j]~isolation_SS2_SS3*fish_SS2_SS3,
          add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7),
          lwd = 1, xaxt="n")
  axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8,
       at =c(1,2,3,5,6,7))
  axis(1,labels = c("Fishless","Fish"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
}
```

![](Community-Variability-Analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Lets also check if any if the variability in any of those variables have
an effect on observed community variability…

``` r
uni_anova_env_com <- data.frame(matrix(NA, nrow = ncol(distances_env_uni), ncol = 3))
for(i in 1:ncol(distances_env_uni)){
  fit <- lmer(beta_deviation_SS2_SS3$observed_distances ~ distances_env_uni[,i] + (1|ID_SS2_SS3))
  uni_anova_env_com[i,1] <- fit@beta[2]
  uni_anova_env_com[i,2] <- round(Anova(fit, test.statistic = "Chisq"),3)$`Pr(>Chisq)`
}
uni_anova_env_com[,3] <- p.adjust(uni_anova_env_com[,2], method = "fdr") 
rownames(uni_anova_env_com) <- colnames(env_data_SS2_SS3_st)
colnames(uni_anova_env_com) <- c("estimate","p","adjust.p")
uni_anova_env_com
```

    ##                   estimate     p  adjust.p
    ## depth          0.061860693 0.015 0.1200000
    ## conductivity   0.001619571 0.935 0.9350000
    ## ph             0.022181534 0.598 0.6834286
    ## do             0.023335011 0.519 0.6834286
    ## turbidity      0.023246876 0.432 0.6834286
    ## clorophylla_a  0.023786955 0.273 0.6834286
    ## phycocyanin    0.032000975 0.215 0.6834286
    ## canopy        -0.016881286 0.574 0.6834286

and beta deviation

``` r
uni_anova_env_dev <- data.frame(matrix(NA, nrow = ncol(distances_env_uni), ncol = 3))
for(i in 1:ncol(distances_env_uni)){
  fit <- lmer(beta_deviation_SS2_SS3$deviation_distances ~ distances_env_uni[,i] + (1|ID_SS2_SS3))
  uni_anova_env_dev[i,1] <- fit@beta[2]
  uni_anova_env_dev[i,2] <- round(Anova(fit, test.statistic = "Chisq"),3)$`Pr(>Chisq)`
}
uni_anova_env_dev[,3] <- p.adjust(uni_anova_env_dev[,2], method = "fdr") 
rownames(uni_anova_env_dev) <- colnames(env_data_SS2_SS3_st)
colnames(uni_anova_env_dev) <- c("estimate","p","adjust.p")
uni_anova_env_dev
```

    ##                  estimate     p  adjust.p
    ## depth          0.35781642 0.358 0.7160000
    ## conductivity   0.12249765 0.672 0.7680000
    ## ph             0.33688171 0.577 0.7680000
    ## do             0.96495138 0.057 0.2960000
    ## turbidity     -0.05039154 0.905 0.9050000
    ## clorophylla_a  0.55212446 0.074 0.2960000
    ## phycocyanin    0.25382971 0.498 0.7680000
    ## canopy        -0.58733323 0.130 0.3466667

We can plot it.

``` r
par(mfrow = c(4,2))
for(i in 1:ncol(distances_env_uni)){
  plot(beta_deviation_SS2_SS3$observed_distances ~ distances_env_uni[,i],
       pch = 16, ylab = "Observed Comm. Variability",
       xlab = "Environmental Variability", main = paste(colnames(distances_env_uni)[i]))
  
}
```

![](Community-Variability-Analyses_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
par(mfrow = c(4,2))
for(i in 1:ncol(distances_env_uni)){
  plot(beta_deviation_SS2_SS3$deviation_distances ~ distances_env_uni[,i],
       pch = 16, ylab = "Beta Deviation",
       xlab = "Environmental Variability", main = paste(colnames(distances_env_uni)[i]))
  
}
```

![](Community-Variability-Analyses_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->
