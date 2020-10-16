Alpha Diversity Analyses
================
Rodolfo Pelinson
16/10/2020

``` r
install.packages("devtools")
devtools::install_github("RodolfoPelinson/Pelinson.et.al.2020B")
```

``` r
library(Pelinson.et.al.2020B)
```

## Analysis of Alpha Diversity

Lets load the necessary packages:

``` r
library(lme4)
library(emmeans)
library(car)
```

### Whole Community

First, lets load the necessary data:

``` r
data(com_SS2_SS3_richness,
     isolation_SS2_SS3,
     fish_SS2_SS3,
     SS_SS2_SS3,
     ID_SS2_SS3,
     All)
```

Now, lets check what probability distribution should we choose using the
most complex model we have:

``` r
mix_model_G <- lmer(com_SS2_SS3_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), REML = F, control = lmerControl(optimizer = "bobyqa"))
mix_model_P <- glmer(com_SS2_SS3_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), family = "poisson", control = glmerControl(optimizer = "bobyqa"))
mix_model_NB <- glmer.nb(com_SS2_SS3_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = glmerControl(optimizer = "bobyqa"))

plot(mix_model_G)
```

![](Alpha-Diversity-Analyses_files/figure-gfm/checking_distribution-1.png)<!-- -->

``` r
plot(mix_model_P)
```

![](Alpha-Diversity-Analyses_files/figure-gfm/checking_distribution-2.png)<!-- -->

``` r
plot(mix_model_NB)
```

![](Alpha-Diversity-Analyses_files/figure-gfm/checking_distribution-3.png)<!-- -->

``` r
AIC(mix_model_G,mix_model_P,mix_model_NB)
```

    ##              df      AIC
    ## mix_model_G  14 190.7859
    ## mix_model_P  13 222.2224
    ## mix_model_NB 14 224.2226

Analysing the data:

``` r
mix_model_richness_NB <- lmer(com_SS2_SS3_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "bobyqa"))
Anova(mix_model_richness_NB, test.statistic = "Chisq")
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: com_SS2_SS3_richness
    ##                                             Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                              12.7360  1  0.0003587 ***
    ## isolation_SS2_SS3                          0.3934  2  0.8214269    
    ## SS_SS2_SS3                                 0.3189  1  0.5722735    
    ## fish_SS2_SS3:isolation_SS2_SS3             2.9020  2  0.2343391    
    ## fish_SS2_SS3:SS_SS2_SS3                    0.0119  1  0.9129797    
    ## isolation_SS2_SS3:SS_SS2_SS3               4.6086  2  0.0998285 .  
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  4.2128  2  0.1216776    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

It seems like there is only an effect of the presence of fish on species
richness.

Lets plot it:

``` r
boxplot(com_SS2_SS3_richness~fish_SS2_SS3, outline = F, ylab = "Richness", xlab = "", at = c(1,2), lwd = 1.5, col = "transparent", xaxt="n")
mylevels <- levels(All)
levelProportions <- summary(All)/length(com_SS2_SS3_richness)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,1,1,2,2,2,1,1,1,2,2,2)[i]
  thislevel <- mylevels[i]
  thisvalues <- com_SS2_SS3_richness[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.3)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(com_SS2_SS3_richness~fish_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2), lwd = 1.5, xaxt="n")
axis(1,labels = c("Fishless","Fish"), cex.axis = 1.5, at =c(1,2))

box(lwd = 2.5)
```

![](Alpha-Diversity-Analyses_files/figure-gfm/Plot_rch-1.png)<!-- -->

### Only Predatory Insects Community

First, lets load the necessary data:

``` r
data(com_SS2_SS3_predators_richness)
```

Analysing the data:

``` r
mix_model_predators_NB <- lmer(com_SS2_SS3_predators_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "bobyqa"))
Anova(mix_model_predators_NB, test.statistic = "Chisq")
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: com_SS2_SS3_predators_richness
    ##                                             Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                              20.2584  1  6.765e-06 ***
    ## isolation_SS2_SS3                         13.6283  2   0.001098 ** 
    ## SS_SS2_SS3                                 0.6010  1   0.438201    
    ## fish_SS2_SS3:isolation_SS2_SS3             3.3606  2   0.186321    
    ## fish_SS2_SS3:SS_SS2_SS3                    1.4155  1   0.234153    
    ## isolation_SS2_SS3:SS_SS2_SS3               3.2956  2   0.192477    
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 11.1988  2   0.003700 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Now pairwise differences:

``` r
emmeans(mix_model_predators_NB, list(pairwise ~ isolation_SS2_SS3), adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of isolation_SS2_SS3`
    ##  isolation_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  30                  4.72 0.287 16.9     3.96     5.48
    ##  120                 3.21 0.287 16.9     2.45     3.97
    ##  480                 3.53 0.295 18.1     2.75     4.30
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3, SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3`
    ##  contrast  estimate    SE   df t.ratio p.value
    ##  30 - 120     1.513 0.406 16.9  3.729  0.0050 
    ##  30 - 480     1.193 0.412 17.5  2.899  0.0290 
    ##  120 - 480   -0.319 0.412 17.5 -0.775  0.8323 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3, SS_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

``` r
emmeans(mix_model_predators_NB, list(pairwise ~ fish_SS2_SS3|isolation_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## $`emmeans of fish_SS2_SS3 | isolation_SS2_SS3, SS_SS2_SS3`
    ## isolation_SS2_SS3 = 30, SS_SS2_SS3 = 2:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent         5.25 0.484 29.4    4.109     6.39
    ##  present        3.50 0.484 29.4    2.359     4.64
    ## 
    ## isolation_SS2_SS3 = 120, SS_SS2_SS3 = 2:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent         4.25 0.484 29.4    3.109     5.39
    ##  present        1.75 0.484 29.4    0.609     2.89
    ## 
    ## isolation_SS2_SS3 = 480, SS_SS2_SS3 = 2:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent         4.25 0.484 29.4    3.109     5.39
    ##  present        3.25 0.484 29.4    2.109     4.39
    ## 
    ## isolation_SS2_SS3 = 30, SS_SS2_SS3 = 3:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent         4.75 0.484 29.4    3.609     5.89
    ##  present        5.39 0.560 31.4    4.072     6.70
    ## 
    ## isolation_SS2_SS3 = 120, SS_SS2_SS3 = 3:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent         4.25 0.484 29.4    3.109     5.39
    ##  present        2.59 0.560 31.4    1.272     3.90
    ## 
    ## isolation_SS2_SS3 = 480, SS_SS2_SS3 = 3:
    ##  fish_SS2_SS3 emmean    SE   df lower.CL upper.CL
    ##  absent         4.64 0.560 31.4    3.325     5.95
    ##  present        1.97 0.560 31.4    0.659     3.29
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of fish_SS2_SS3 | isolation_SS2_SS3, SS_SS2_SS3`
    ## isolation_SS2_SS3 = 30, SS_SS2_SS3 = 2:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    1.750 0.685 29.4  2.556  0.0160 
    ## 
    ## isolation_SS2_SS3 = 120, SS_SS2_SS3 = 2:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    2.500 0.685 29.4  3.652  0.0010 
    ## 
    ## isolation_SS2_SS3 = 480, SS_SS2_SS3 = 2:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    1.000 0.685 29.4  1.461  0.1547 
    ## 
    ## isolation_SS2_SS3 = 30, SS_SS2_SS3 = 3:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present   -0.637 0.740 30.7 -0.860  0.3963 
    ## 
    ## isolation_SS2_SS3 = 120, SS_SS2_SS3 = 3:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    1.663 0.740 30.7  2.248  0.0319 
    ## 
    ## isolation_SS2_SS3 = 480, SS_SS2_SS3 = 3:
    ##  contrast         estimate    SE   df t.ratio p.value
    ##  absent - present    2.667 0.792 31.4  3.369  0.0020 
    ## 
    ## Degrees-of-freedom method: kenward-roger

It seems like there is an effect of fish and isolation on predatory
insect richness. Also, the effect of fish stronger in low isolation
treatments in the second survey, and stronger in high isolation
treatments in the third survey.

Lets plot it:

``` r
boxplot(com_SS2_SS3_predators_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3, outline = F, ylab = "Richness", xlab = "", lwd = 1.5, ylim = c(0,8), col = "transparent", main = "Predators", at = c(1,2,4,5,7,8,10,11,13,14,16,17), xaxt="n")

All2 <- factor(All, levels = c("2 030 absent","2 030 present","2 120 absent","2 120 present","2 480 absent","2 480 present",
                               "3 030 absent","3 030 present","3 120 absent","3 120 present","3 480 absent","3 480 present"))

mylevels <- levels(All2)
levelProportions <- summary(All2)/length(com_SS2_SS3_predators_richness)
col <- c(rep(c("sienna3","dodgerblue3"),3), rep("grey70",6))
bg <- c(rep(c("sienna3","dodgerblue3"),3),rep(c("sienna3","dodgerblue3"),3))
pch <- c(15,15,16,16,17,17,22,22,21,21,24,24)
for(i in 1:length(mylevels)){

  x<- c(1,2,4,5,7,8,10,11,13,14,16,17)[i]
  thislevel <- mylevels[i]
  thisvalues <- com_SS2_SS3_predators_richness[All==thislevel]

  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.8)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3)

}
boxplot(com_SS2_SS3_predators_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3, add = T, outline = F,at = c(1,2,4,5,7,8,10,11,13,14,16,17), lwd = 1.5, col = "transparent", xaxt="n")
axis(1,labels = c("Abs.","Pre.", "Abs.","Pre.","Abs.","Pre.","Abs.","Pre.", "Abs.","Pre.","Abs.","Pre."), cex.axis = 0.8, at =c(1,2,4,5,7,8,10,11,13,14,16,17))
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 1, at =c(1.5,4.5,7.5,10.5,13.5,16.5), line = 1.5, tick = F )
axis(1,labels = c("Second","Third"), cex.axis = 1.2, at =c(4.5,13.5), line = 3, tick = F )
box(lwd = 2.5)
```

![](Alpha-Diversity-Analyses_files/figure-gfm/Plotting_rch_pr-1.png)<!-- -->

### Only Non-Predatory Insects Community

First, lets load the necessary data:

``` r
data(com_SS2_SS3_non_predators_richness)
```

Analysing the data:

``` r
mix_model_non_predators <- lmer(com_SS2_SS3_non_predators_richness~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), control = lmerControl(optimizer = "bobyqa"))
Anova(mix_model_non_predators, test.statistic = "Chisq")
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: com_SS2_SS3_non_predators_richness
    ##                                            Chisq Df Pr(>Chisq)  
    ## fish_SS2_SS3                              0.9749  1    0.32347  
    ## isolation_SS2_SS3                         5.5420  2    0.06260 .
    ## SS_SS2_SS3                                0.0438  1    0.83423  
    ## fish_SS2_SS3:isolation_SS2_SS3            0.7567  2    0.68499  
    ## fish_SS2_SS3:SS_SS2_SS3                   0.2384  1    0.62538  
    ## isolation_SS2_SS3:SS_SS2_SS3              6.5742  2    0.03736 *
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3 4.8928  2    0.08661 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Now pairwise differences:

``` r
emmeans(mix_model_non_predators, list(pairwise ~ isolation_SS2_SS3|SS_SS2_SS3), adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                  7.00 0.549 32     5.62     8.38
    ##  120                 6.75 0.549 32     5.37     8.13
    ##  480                 7.38 0.549 32     5.99     8.76
    ## 
    ## SS_SS2_SS3 = 3:
    ##  isolation_SS2_SS3 emmean    SE df lower.CL upper.CL
    ##  30                  5.58 0.599 32     4.07     7.09
    ##  120                 8.08 0.599 32     6.57     9.59
    ##  480                 7.67 0.646 32     6.04     9.29
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 3 estimates 
    ## 
    ## $`pairwise differences of isolation_SS2_SS3 | SS_SS2_SS3`
    ## SS_SS2_SS3 = 2:
    ##  contrast  estimate    SE df t.ratio p.value
    ##  30 - 120     0.250 0.777 32  0.322  0.9843 
    ##  30 - 480    -0.375 0.777 32 -0.483  0.9504 
    ##  120 - 480   -0.625 0.777 32 -0.805  0.8119 
    ## 
    ## SS_SS2_SS3 = 3:
    ##  contrast  estimate    SE df t.ratio p.value
    ##  30 - 120    -2.500 0.848 32 -2.949  0.0176 
    ##  30 - 480    -2.083 0.881 32 -2.365  0.0710 
    ##  120 - 480    0.417 0.881 32  0.473  0.9531 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: sidak method for 3 tests

It seems like there is an increase in richness from the low to
intermediate isolation, but not to high isolation. However, it only
happened in the last survey.

Lets plot it:

``` r
boxplot(com_SS2_SS3_non_predators_richness~isolation_SS2_SS3*SS_SS2_SS3, outline = F, ylab = "Richness", xlab = "", at = c(1,2,3,5,6,7), lwd = 1.5, col = "transparent", xaxt="n", main = "Non Predators")
mylevels <- levels(All)
levelProportions <- summary(All)/length(com_SS2_SS3_non_predators_richness)
col <- c(rep("sienna3",3), rep("dodgerblue3",3), rep("grey70",6))
bg <- c(rep("sienna3",3), rep("dodgerblue3",3),rep("sienna3",3), rep("dodgerblue3",3))
pch <- c(15,16,17,15,16,17,22,21,24,22,21,24)
for(i in 1:length(mylevels)){
  
  x<- c(1,2,3,1,2,3,5,6,7,5,6,7)[i]
  thislevel <- mylevels[i]
  thisvalues <- com_SS2_SS3_non_predators_richness[All==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(x, length(thisvalues)), amount=levelProportions[i]/0.3)
  points(myjitter, thisvalues, pch=pch[i], col=col[i], bg = bg[i] , cex = 1.5, lwd = 3) 
  
}
boxplot(com_SS2_SS3_non_predators_richness~isolation_SS2_SS3*SS_SS2_SS3, add = T, col = "transparent", outline = F,at = c(1,2,3,5,6,7), lwd = 1.5,  xaxt="n") 
axis(1,labels = c("30 m","120 m", "480 m","30 m","120 m", "480 m"), cex.axis = 0.8, at =c(1,2,3,5,6,7))
axis(1,labels = c("Second","Third"), cex.axis = 1, at =c(2,6), line = 1.5, tick = F )
box(lwd = 2.5)
```

![](Alpha-Diversity-Analyses_files/figure-gfm/Plot_rch_non_pr-1.png)<!-- -->
