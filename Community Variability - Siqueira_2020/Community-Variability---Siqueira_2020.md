Community Variability - Siqueira et al 2020
================
Rodolfo Pelinson
20/10/2020

This is the same community variability analyses presented in the main
paper, but using a different null model. The same one as used in
[Siqueira et al. 2020](https://doi.org/10.1002/ecy.3014)

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

### Whole community for the last two surveys.

First loading data

``` r
data(com_SS2_SS3, All, fish_SS2_SS3, isolation_SS2_SS3, SS_SS2_SS3, ID_SS2_SS3)
```

Computing observed and expected distances to centroid, and
beta-deviation.

``` r
beta_deviation_SS2_SS3 <- beta_deviation_siqueira_et_al_2019(com_SS2_SS3, times = 10000,
                                                          transform = NULL, dist = "bray", seed = 2, group = All, keep.gamma = T) 
```

Looking at residual plots for observed, expected distances to centroids
and deviations.

``` r
fit_expected_SS2_SS3_G <- lmer(beta_deviation_SS2_SS3$expected_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
```

    ## boundary (singular) fit: see ?isSingular

``` r
plot(fit_expected_SS2_SS3_G)
```

![](Community-Variability---Siqueira_2020_files/figure-gfm/7-1.png)<!-- -->

``` r
qqnorm(resid(fit_expected_SS2_SS3_G, type = "pearson"))
qqline(resid(fit_expected_SS2_SS3_G, type = "pearson"))
```

![](Community-Variability---Siqueira_2020_files/figure-gfm/7-2.png)<!-- -->

``` r
fit_deviation_SS2_SS3_G <- lmer(beta_deviation_SS2_SS3$deviation_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
plot(fit_deviation_SS2_SS3_G)
```

![](Community-Variability---Siqueira_2020_files/figure-gfm/7-3.png)<!-- -->

``` r
qqnorm(resid(fit_deviation_SS2_SS3_G, type = "pearson"))
qqline(resid(fit_deviation_SS2_SS3_G, type = "pearson"))
```

![](Community-Variability---Siqueira_2020_files/figure-gfm/7-4.png)<!-- -->

The fit of the model for beta deviation here is pretty bad. If we were
to make inferences based on these results, we would have to look for a
better statistical model, maybe considering a different statistical
distribution, which is hard given that we have overdispersion and
negative values (i.e. excludes some distributions and data
transformations). But, since we are not interpreting results based on
this model, we will still use it.

#### Expected Community Variability

Running ANOVA table for expected distances to group centroids, or
expected beta-diversity/community variability in each treatment.

``` r
fit_expected_SS2_SS3 <- lmer(beta_deviation_SS2_SS3$expected_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_expected_SS2_SS3, test.statistic = "F"),3)
```

    ## Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
    ## 
    ## Response: beta_deviation_SS2_SS3$expected_distances
    ##                                               F Df Df.res Pr(>F)  
    ## fish_SS2_SS3                              0.142  1 16.626  0.712  
    ## isolation_SS2_SS3                         1.062  2 16.546  0.368  
    ## SS_SS2_SS3                                7.412  1 16.421  0.015 *
    ## fish_SS2_SS3:isolation_SS2_SS3            4.025  2 16.651  0.038 *
    ## fish_SS2_SS3:SS_SS2_SS3                   0.010  1 16.574  0.921  
    ## isolation_SS2_SS3:SS_SS2_SS3              3.271  2 16.483  0.064 .
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 3.853  2 16.604  0.042 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_expected_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|fish_SS2_SS3|SS_SS2_SS3), adjust = "tukey")
```

    ## $`emmeans of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.182 0.0575 32 0.037047    0.327
    ##  120                0.263 0.0575 32 0.118382    0.408
    ##  480                0.218 0.0575 32 0.073533    0.363
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.204 0.0575 32 0.058754    0.348
    ##  120                0.338 0.0575 32 0.192908    0.482
    ##  480                0.145 0.0575 32 0.000365    0.290
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.256 0.0575 32 0.111312    0.401
    ##  120                0.333 0.0575 32 0.187847    0.477
    ##  480                0.358 0.0675 32 0.187503    0.528
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean     SE df lower.CL upper.CL
    ##  30                 0.544 0.0675 32 0.374068    0.714
    ##  120                0.213 0.0675 32 0.042544    0.383
    ##  480                0.229 0.0675 32 0.058732    0.399
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | fish_SS2_SS3, SS_SS2_SS3`
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   -0.0813 0.0813 32 -1.001  0.5817 
    ##  30 - 480   -0.0365 0.0813 32 -0.449  0.8952 
    ##  120 - 480   0.0448 0.0813 32  0.552  0.8463 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 2:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   -0.1342 0.0813 32 -1.651  0.2397 
    ##  30 - 480    0.0584 0.0813 32  0.718  0.7544 
    ##  120 - 480   0.1925 0.0813 32  2.369  0.0606 
    ## 
    ## fish_SS2_SS3 = absent, SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120   -0.0765 0.0813 32 -0.942  0.6183 
    ##  30 - 480   -0.1015 0.0887 32 -1.145  0.4940 
    ##  120 - 480  -0.0250 0.0887 32 -0.282  0.9571 
    ## 
    ## fish_SS2_SS3 = present, SS_SS2_SS3 = 3:
    ##  contrast  estimate     SE df t.ratio p.value
    ##  30 - 120    0.3315 0.0955 32  3.471  0.0042 
    ##  30 - 480    0.3153 0.0955 32  3.302  0.0065 
    ##  120 - 480  -0.0162 0.0955 32 -0.170  0.9843 
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

Patterns are similar to those observed for the observed distances using
the other null model.

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

![](Community-Variability---Siqueira_2020_files/figure-gfm/9-1.png)<!-- -->

#### Beta-Deviation

Running ANOVA table for the deviations of expected distances to group
centroids from observed distances.

``` r
fit_deviation_SS2_SS3 <- lmer(beta_deviation_SS2_SS3$deviation_distances~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "nlminbwrap"))
round(Anova(fit_deviation_SS2_SS3, test.statistic = "F"),3)
```

    ## Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)
    ## 
    ## Response: beta_deviation_SS2_SS3$deviation_distances
    ##                                               F Df Df.res Pr(>F)  
    ## fish_SS2_SS3                              3.780  1 16.626  0.069 .
    ## isolation_SS2_SS3                         1.117  2 16.546  0.351  
    ## SS_SS2_SS3                                6.966  1 16.421  0.018 *
    ## fish_SS2_SS3:isolation_SS2_SS3            2.413  2 16.651  0.120  
    ## fish_SS2_SS3:SS_SS2_SS3                   5.365  1 16.574  0.034 *
    ## isolation_SS2_SS3:SS_SS2_SS3              3.958  2 16.483  0.039 *
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 0.054  2 16.604  0.947  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
emmeans(fit_deviation_SS2_SS3, list(pairwise ~ fish_SS2_SS3|SS_SS2_SS3), adjust = "tukey")
```

    ## $`emmeans of fish_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  fish_SS2_SS3 emmean   SE df lower.CL upper.CL
    ##  absent         13.2 6.99 32    -3.23     29.6
    ##  present        14.0 6.99 32    -2.44     30.4
    ## 
    ## SS_SS2_SS3 = 3:
    ##  fish_SS2_SS3 emmean   SE df lower.CL upper.CL
    ##  absent         49.5 7.42 32    32.07     66.9
    ##  present        16.1 8.21 32    -3.20     35.3
    ## 
    ## Results are averaged over the levels of: isolation_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of fish_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast         estimate    SE df t.ratio p.value
    ##  absent - present   -0.789  9.88 32 -0.080  0.9369 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast         estimate    SE df t.ratio p.value
    ##  absent - present   33.400 11.07 32  3.018  0.0050 
    ## 
    ## Results are averaged over the levels of: isolation_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger

``` r
emmeans(fit_deviation_SS2_SS3, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3), adjust = "tukey")
```

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                 18.41  8.56 32    -3.15     40.0
    ##  120                13.71  8.56 32    -7.85     35.3
    ##  480                 8.56  8.56 32   -13.00     30.1
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                 13.04  9.34 32   -10.49     36.6
    ##  120                30.75  9.34 32     7.22     54.3
    ##  480                54.53 10.06 32    29.20     79.9
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast  estimate   SE df t.ratio p.value
    ##  30 - 120      4.70 12.1 32  0.388  0.9205 
    ##  30 - 480      9.85 12.1 32  0.813  0.6976 
    ##  120 - 480     5.15 12.1 32  0.425  0.9054 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast  estimate   SE df t.ratio p.value
    ##  30 - 120    -17.71 13.2 32 -1.341  0.3834 
    ##  30 - 480    -41.50 13.7 32 -3.024  0.0132 
    ##  120 - 480   -23.79 13.7 32 -1.733  0.2086 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: tukey method for comparing a family of 3 estimates

Beta-deviation seems to be greater in fishless ponds, and in the higher
isolation treatment. **But, as we know by inspecting the residuals, this
model have a very poor fit.**

Plotting it:

``` r
boxplot(beta_deviation_SS2_SS3$deviation_distances~isolation_SS2_SS3*fish_SS2_SS3, outline = F, ylab = "Beta-Deviation", xlab = "", at = c(1,2,3,5,6,7),ylim = c(-20,180), lwd = 1.5, col = "transparent", xaxt="n")
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

![](Community-Variability---Siqueira_2020_files/figure-gfm/11-1.png)<!-- -->
