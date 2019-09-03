
# TO DO:

code.dir = "~/Dropbox/Personal computer/Independent studies/Nonparametric Phat (NPPhat)/Linked to OSF (NPPhat)/Code (on git)/Applied example"
setwd(code.dir)
source("temp_helper.R")  # to be replaced by calling analysis_RRR.R
source("helper_applied.R")


############################################ READ IN AND DO SANITY CHECKS ############################################

data.dir = "~/Dropbox/Personal computer/Independent studies/Nonparametric Phat (NPPhat)/Linked to OSF (NPPhat)/Applied example data"
results.dir = "~/Dropbox/Personal computer/Independent studies/Nonparametric Phat (NPPhat)/Linked to OSF (NPPhat)/Applied example results"
write.results = TRUE  # should we overwrite results?
setwd(data.dir)

# we scraped data from Croke's Figures 1-2
dm0 = read.csv("Croke_data.csv")
dt0 = read.csv("taylor_data.csv")


library(metafor)
library(MetaUtility)
library(boot)

dm = scrape_meta( type = "raw", 
                  est = dm0$est,
                  hi = dm0$hi)
dm$study = dm0$study

dt = scrape_meta( type = "raw", 
                  est = dt0$est,
                  hi = dt0$hi )
dt$study = dt0$study

# full sample
# reported: 0.13 [0.03, 0.24]
meta.m = rma.uni( yi = dm$yi,
                  vi = dm$vyi, 
                  method = "REML", 
                  knha = TRUE)
# matches if using DL

# TMSDG sample
# reported: 0.08 [-0.11, 0.27]
meta.t = rma.uni( yi = dt$yi,
                  vi = dt$vyi, 
                  method = "REML",
                  knha = TRUE)
# matches if using DL


############################################ FOREST PLOTS ############################################

forest(meta.m,
       slab = dm$study,
       xlab = "Estimated increase in body weight (kg)")
# 8 x 7" works well

forest(meta.t,
       slab = dt$study,
       xlab = "Estimated increase in body weight (kg)")

############################################ DENSITY PLOTS ############################################

library(ggplot2)

# both on one plot
yi.std.m = (dm$yi - c(meta.m$b)) / sqrt(c(meta.m$tau2) + dm$vyi)
yi.std.t = (dt$yi - c(meta.m$b)) / sqrt(c(meta.m$tau2) + dt$vyi)

# LOOKS GREAT! :)
p = ggplot( ) +
  
  geom_vline(xintercept = 0,
             color = "gray") +

  geom_density( data = data.frame(yi.std.m), 
                aes(x = yi.std.m,
                    lty = "Croke") ) +
  
  geom_density( data = data.frame(yi.std.t), 
                aes(x = yi.std.t,
                    lty = "Taylor-Robinson")) +
  
  labs(linetype="Meta-analysis") +
  
  xlab("Standardized point estimates") +
  ylab("Kernel density estimate") +
  
  scale_x_continuous( breaks = seq(-2, 3, 0.5) ) +
  
  theme_classic()

plot(p)
# 7 x 5



############################### ANALYZE BOTH METAS - PARAMETRIC AND NP ############################### 

# to run main analyses vs. supplementary analyses, just change the pub.bias parameter
#  then run the below script
# it will automatically name the files appropriately
pub.bias = FALSE
ql = as.list( c( 0, .1, .2, 0.5, -0.1, -0.2 ) )
boot.reps = 2000  # ~~~ increase later
if( exists("resE") ) rm(resE)


##### Croke #####
# bm
analyze_one_meta( dat = dm,
                  meta.name = "Croke", 
                  pub.bias = FALSE,
                  method = "parametric" )

analyze_one_meta( dat = dm,
                  meta.name = "Croke", 
                  pub.bias = FALSE,
                  method = "calibrated" )

##### TMSDG #####
analyze_one_meta( dat = dt,
                  meta.name = "Taylor-Robinson", 
                  pub.bias = FALSE,
                  method = "parametric" )

analyze_one_meta( dat = dt,
                  meta.name = "Taylor-Robinson", 
                  pub.bias = FALSE,
                  method = "calibrated" )

View(resE)

##### Write Results #####
if ( write.results == TRUE ) {
  setwd(results.dir)

  write.csv( t(resE), 
             paste( "results_table.csv", sep = "" ),
             row.names = FALSE )
  
  # for copy-pasting into LaTeX
  library(xtable)
  print( xtable(t(resE)), file = paste( "results_xtable.txt", sep = "" ) )
  
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           MANUAL SANITY CHECKS                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

############################################ PARAMETRIC ANALYSIS ############################################

library(MetaUtility)
# dm dataset seems nonnormal :) 
q = 0.1
Phat.m = MetaUtility::prop_stronger( q = q,
                                     M = meta.m$b,
                                     t2 = meta.m$tau2,
                                     se.M = meta.m$se,
                                     se.t2 = meta.m$se.tau2,
                                     tail = "above",
                                     
                                     # in order to get Shapiro p-value
                                     dat = dm,
                                     yi.name = "yi",
                                     vi.name = "vyi")

Phat.t = MetaUtility::prop_stronger( q = q,
                                     M = meta.t$b,
                                     t2 = meta.t$tau2,
                                     se.M = meta.t$se,
                                     se.t2 = meta.t$se.tau2,
                                     tail = "above",
                                     
                                     # in order to get Shapiro p-value
                                     dat = dt,
                                     yi.name = "yi",
                                     vi.name = "vyi")

# either way, it's 47-57% above q = 0.10
# :) 
# so can still agree with Croke about power, but this is an important example
# where "statistical significance" was part of the debate
# so that also helps with response to Ferguson


############################################ NONPARAMETRIC ANALYSIS ############################################

# check heterogeneity ratios
sqrt(meta.m$tau2) / meta.m$b
sqrt(meta.t$tau2) / meta.t$b
# definitely fine to use our metrics

q = 0.1
library(dplyr)
library(purrr)
PhatNP.m = prop_stronger_np(q = q,
                            yi = dm$yi,
                            vi = dm$vyi,
                            CI.level = 0.95, 
                            tail = "above",
                            R = 2000,
                            return.vectors = TRUE )
# compare to parametric
# point estimates of 48% vs. 54% (parametric is higher)
PhatNP.m$res; Phat.m

PhatNP.t = prop_stronger_np(q = q,
                            yi = dt$yi,
                            vi = dt$vyi,
                            CI.level = 0.95, 
                            tail = "above",
                            R = 2000,
                            return.vectors = TRUE )
# compare to parametric
# NP CI width is actually better
# point estimates of 34% vs. 48% (parametric is higher)
PhatNP.t$res; Phat.t


#### ~~ TEMP ONLY

boot.res.ens = boot( data = dt, 
                     parallel = "multicore",
                     R = boot.reps, 
                     statistic = function(original, indices) {
                       
                       b = original[indices,]
                       
                       ens.b = my_ens( yi = b$yi, 
                                       sei = sqrt(b$vyi) )
                       return( sum(ens.b > 0) / length(ens.b) )
                     }
)

bootCIs.ens = boot.ci(boot.res.ens, type="bca")
boot.lo.ens = bootCIs.ens$bca[4]
boot.hi.ens = bootCIs.ens$bca[5]
