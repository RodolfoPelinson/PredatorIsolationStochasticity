Community Size Analyses
================
Rodolfo Pelinson
16/10/2020

``` r
devtools::install_github("RodolfoPelinson/Pelinson.et.al.2020B")
```

``` r
library(Pelinson.et.al.2020B)
```

## Analysis of Community Size

First, lets load the necessary packages:

``` r
library(lme4)
```

    ## Loading required package: Matrix

``` r
library(emmeans)
library(car)
```

    ## Loading required package: carData

    ## Registered S3 methods overwritten by 'car':
    ##   method                          from
    ##   influence.merMod                lme4
    ##   cooks.distance.influence.merMod lme4
    ##   dfbeta.influence.merMod         lme4
    ##   dfbetas.influence.merMod        lme4

Now, lets load the necessary data:

``` r
data(com_SS2_SS3_abundance,
     isolation_SS2_SS3,
     fish_SS2_SS3,
     SS_SS2_SS3,
     ID_SS2_SS3,
     All)
```

Now, lets check what probability distribution should we choose using the
most complex model we have:

``` r
mix_model_G <- lmer(com_SS2_SS3_abundance~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), REML = F)
mix_model_P <- glmer(com_SS2_SS3_abundance~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3), family = "poisson")
mix_model_NB <- glmer.nb(com_SS2_SS3_abundance~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
```

    ## Warning in glmer.nb(com_SS2_SS3_abundance ~ fish_SS2_SS3 * isolation_SS2_SS3 * :
    ## no 'data = *' in glmer.nb() call ... Not much is guaranteed

``` r
plot(mix_model_G)
```

![](Community-Size-Analyses_files/figure-gfm/checking%20distribution-1.png)<!-- -->

``` r
plot(mix_model_P)
```

![](Community-Size-Analyses_files/figure-gfm/checking%20distribution-2.png)<!-- -->

``` r
plot(mix_model_NB)
```

![](Community-Size-Analyses_files/figure-gfm/checking%20distribution-3.png)<!-- -->

``` r
AIC(mix_model_G,mix_model_P,mix_model_NB)
```

    ##              df       AIC
    ## mix_model_G  14  623.8594
    ## mix_model_P  13 1313.5127
    ## mix_model_NB 14  580.8158

Looking at AIC values and the spread of residuals, it seems like
negative binomial is the best option here.

So we can go further analysing the data:

``` r
mix_model_NB <- glmer.nb(com_SS2_SS3_abundance~fish_SS2_SS3*isolation_SS2_SS3*SS_SS2_SS3 + (1|ID_SS2_SS3))
```

    ## Warning in glmer.nb(com_SS2_SS3_abundance ~ fish_SS2_SS3 * isolation_SS2_SS3 * :
    ## no 'data = *' in glmer.nb() call ... Not much is guaranteed

``` r
Anova(mix_model_NB, test.statistic = "Chisq")
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: com_SS2_SS3_abundance
    ##                                             Chisq Df Pr(>Chisq)    
    ## fish_SS2_SS3                               0.0131  1   0.908778    
    ## isolation_SS2_SS3                          0.3226  2   0.851039    
    ## SS_SS2_SS3                                37.2252  1  1.052e-09 ***
    ## fish_SS2_SS3:isolation_SS2_SS3             1.6093  2   0.447243    
    ## fish_SS2_SS3:SS_SS2_SS3                    4.3397  1   0.037233 *  
    ## isolation_SS2_SS3:SS_SS2_SS3              12.6888  2   0.001757 ** 
    ## fish_SS2_SS3:isolation_SS2_SS3:SS_SS2_SS3  0.5970  2   0.741948    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Now pairwise differences:

``` r
emmeans(mix_model_NB, list(pairwise ~ SS_SS2_SS3), adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of SS_SS2_SS3`
    ##  SS_SS2_SS3 emmean    SE  df asymp.LCL asymp.UCL
    ##  2            5.06 0.122 Inf      4.79      5.34
    ##  3            5.96 0.138 Inf      5.65      6.27
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3, isolation_SS2_SS3 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of SS_SS2_SS3`
    ##  contrast estimate   SE  df z.ratio p.value
    ##  2 - 3      -0.892 0.15 Inf -5.942  <.0001 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3, isolation_SS2_SS3 
    ## Results are given on the log (not the response) scale.

``` r
emmeans(mix_model_NB, list(pairwise ~ SS_SS2_SS3|fish_SS2_SS3), adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of SS_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  SS_SS2_SS3 emmean    SE  df asymp.LCL asymp.UCL
    ##  2            4.94 0.172 Inf      4.56      5.33
    ##  3            6.14 0.183 Inf      5.73      6.55
    ## 
    ## fish_SS2_SS3 = present:
    ##  SS_SS2_SS3 emmean    SE  df asymp.LCL asymp.UCL
    ##  2            5.18 0.173 Inf      4.80      5.57
    ##  3            5.77 0.202 Inf      5.32      6.22
    ## 
    ## Results are averaged over the levels of: isolation_SS2_SS3 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of SS_SS2_SS3 | fish_SS2_SS3`
    ## fish_SS2_SS3 = absent:
    ##  contrast estimate    SE  df z.ratio p.value
    ##  2 - 3       -1.20 0.202 Inf -5.930  <.0001 
    ## 
    ## fish_SS2_SS3 = present:
    ##  contrast estimate    SE  df z.ratio p.value
    ##  2 - 3       -0.59 0.220 Inf -2.685  0.0073 
    ## 
    ## Results are averaged over the levels of: isolation_SS2_SS3 
    ## Results are given on the log (not the response) scale.

``` r
emmeans(mix_model_NB, list(pairwise ~ SS_SS2_SS3|isolation_SS2_SS3), adjust = "sidak")
```

    ## NOTE: Results may be misleading due to involvement in interactions

    ## $`emmeans of SS_SS2_SS3 | isolation_SS2_SS3`
    ## isolation_SS2_SS3 = 30:
    ##  SS_SS2_SS3 emmean    SE  df asymp.LCL asymp.UCL
    ##  2            5.25 0.210 Inf      4.78      5.72
    ##  3            5.65 0.238 Inf      5.12      6.18
    ## 
    ## isolation_SS2_SS3 = 120:
    ##  SS_SS2_SS3 emmean    SE  df asymp.LCL asymp.UCL
    ##  2            5.09 0.213 Inf      4.62      5.57
    ##  3            5.75 0.228 Inf      5.24      6.26
    ## 
    ## isolation_SS2_SS3 = 480:
    ##  SS_SS2_SS3 emmean    SE  df asymp.LCL asymp.UCL
    ##  2            4.84 0.211 Inf      4.37      5.31
    ##  3            6.47 0.239 Inf      5.93      7.01
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: sidak method for 2 estimates 
    ## 
    ## $`pairwise differences of SS_SS2_SS3 | isolation_SS2_SS3`
    ## isolation_SS2_SS3 = 30:
    ##  contrast estimate    SE  df z.ratio p.value
    ##  2 - 3      -0.394 0.264 Inf -1.489  0.1365 
    ## 
    ## isolation_SS2_SS3 = 120:
    ##  contrast estimate    SE  df z.ratio p.value
    ##  2 - 3      -0.654 0.249 Inf -2.631  0.0085 
    ## 
    ## isolation_SS2_SS3 = 480:
    ##  contrast estimate    SE  df z.ratio p.value
    ##  2 - 3      -1.630 0.262 Inf -6.229  <.0001 
    ## 
    ## Results are averaged over the levels of: fish_SS2_SS3 
    ## Results are given on the log (not the response) scale.

It seems that community size grows with time. But it grows larger in
fishless ponds, and in higher isolation treatments.

Lets plot it:

``` r
isolation_SS <-rep(NA, length(All))
for(i in 1:length(All)){
  if(SS_SS2_SS3[i] == "2" & isolation_SS2_SS3[i] == "30"){isolation_SS[i] <- "2 030"}
  if(SS_SS2_SS3[i] == "2" & isolation_SS2_SS3[i] == "120"){isolation_SS[i] <- "2 120"}
  if(SS_SS2_SS3[i] == "2" & isolation_SS2_SS3[i] == "480"){isolation_SS[i] <- "2 480"}
  if(SS_SS2_SS3[i] == "3" & isolation_SS2_SS3[i] == "30"){isolation_SS[i] <- "3 030"}
  if(SS_SS2_SS3[i] == "3" & isolation_SS2_SS3[i] == "120"){isolation_SS[i] <- "3 120"}
  if(SS_SS2_SS3[i] == "3" & isolation_SS2_SS3[i] == "480"){isolation_SS[i] <- "3 480"}
}
isolation_SS <- as.factor(isolation_SS)
isolation_SS <- factor(isolation_SS, levels = c("2 030", "3 030","2 120", "3 120","2 480", "3 480" ))
new_All <- factor(All, levels = c("2 030 absent","3 030 absent","2 120 absent","3 120 absent","2 480 absent","3 480 absent",
                                  "2 030 present","3 030 present","2 120 present","3 120 present","2 480 present","3 480 present"))

boxplot(com_SS2_SS3_abundance~isolation_SS, outline = F, ylab = "Abundance", xlab = "", at = c(1,2,4,5,7,8), lwd = 1.5, ylim = c(0,1400), col = "transparent", xaxt="n")
mylevels <- levels(new_All)
levelProportions <- summary(new_All)/length(com_SS2_SS3_abundance)
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
axis(1,labels = c("1st", "2nd", "1st", "2nd", "1st", "2nd"), cex.axis = 0.8, at =c(1,2,4,5,7,8))
axis(1,labels = c("30 m","120 m", "480 m"), cex.axis = 1, at =c(1.5,4.5,7.5), line = 1.5, tick = F )
boxplot(com_SS2_SS3_abundance~isolation_SS, add = T, col = "transparent", outline = F,at = c(1,2,4,5,7,8), lwd = 1.5, xaxt="n")
box(lwd = 2.5)
```

![](Community-Size-Analyses_files/figure-gfm/Plotting%20effect%20of%20abundance-1.png)<!-- -->
