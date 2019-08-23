
# TO DO:
#  - Make forest plots
#  - Make table as in Hedges response for different values of q (but only 1 table OR figure in text)
#   and show parametric for all of them

# temporarily taken from sim study

code.dir = "~/Dropbox/Personal computer/Independent studies/Nonparametric Phat (NPPhat)/Linked to OSF (NPPhat)/Code (on git)/Applied example"
setwd(code.dir)
source("temp_helper.R")

############################################ READ IN AND DO SANITY CHECKS ############################################

data.dir = "~/Dropbox/Personal computer/Independent studies/Nonparametric Phat (NPPhat)/Linked to OSF (NPPhat)/Applied example data"

setwd(data.dir)
dm = read.csv("miguel_data.csv")
dt = read.csv("taylor_data.csv")


library(metafor)
library(MetaUtility)

dm = scrape_meta( type = "raw", 
                  est = dm$est,
                  hi = dm$hi )

dt = scrape_meta( type = "raw", 
                  est = dt$est,
                  hi = dt$hi )

# full sample
# reported: 0.13 [0.03, 0.24]
meta.m = rma.uni( yi = dm$yi,
                  vi = dm$vyi, 
                  method = "REML" )
# matches if using DL

# TMSDG sample
# reported: 0.08 [-0.11, 0.27]
meta.t = rma.uni( yi = dt$yi,
                  vi = dt$vyi, 
                  method = "REML" )
# matches if using DL


############################################ FOREST PLOTS ############################################

forest(meta.m)
forest(meta.t)

############################################ DENSITY PLOTS ############################################

library(ggplot2)

# plot_yi_std = function(dat) {
#   
#   meta = rma.uni( yi = dat$yi,
#                     vi = dat$vyi, 
#                     method = "REML" )
#   
#   yi.std = (dat$yi - c(meta$b)) / sqrt(c(meta$tau2) + dat$vyi)
#   
#   p = ggplot( data = data.frame(yi.std), 
#               aes(x = yi.std) ) +
#     geom_density() +
#     geom_histogram(alpha = 0.4) +
#     theme_minimal()
#   
#   plot(p)
#   return(p)
# }
# 
# plot_yi_std(dm)
# 


# both on one plot

yi.std.m = (dm$yi - c(meta.m$b)) / sqrt(c(meta.m$tau2) + dm$vyi)
yi.std.t = (dt$yi - c(meta.m$b)) / sqrt(c(meta.m$tau2) + dt$vyi)

# LOOKS GREAT! :)
p = ggplot( ) +
  
  geom_vline(xintercept = 0,
             color = "gray") +
  # geom_vline(xintercept = 0.1,
  #            color = "red") +
  # geom_vline(xintercept = 0.2,
  #            color = "red") +
  
  geom_density( data = data.frame(yi.std.m), 
                aes(x = yi.std.m,
                    lty = "Miguel") ) +
  
  geom_density( data = data.frame(yi.std.t), 
                aes(x = yi.std.t,
                    lty = "Taylor")) +
  
  labs(linetype="Meta-analysis") +
  
  xlab("Standardized point estimates") +
  ylab("Kernel density estimate") +
  
  scale_x_continuous( breaks = seq(-2, 3, 0.5) ) +
  
  theme_classic()

plot(p)

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
# so can still agree with Miguel about power, but this is an important example
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


############################################ FROM "METAWARS" ############################################

# thresholds of scientific importance on z or standardized-bhat scale
ql = as.list( c( 0, .1, .2 ) )

# should we overwrite existing results?
write.results = FALSE


boot.reps = 1000

# for rounding
digits = 2

library(metafor)
library(weightr)
library(MetaUtility)

# I am re-analyzing them all using the same model
# Note that all of these have potentially correlated point estimates
# We could get cluster ID in order to fit a robust model for Prescott and
#  Anderson, but probably not Ferguson.
# See sensitivity analyses at the very end. 


############################### VIDEO GAME MAIN ANALYSES: PROPORTION ABOVE ############################### 

# to run main analyses vs. supplementary analyses, just change the pub.bias parameter
#  then run the below script
# it will automatically name the files appropriately
pub.bias = FALSE
ql = as.list( c( 0, .1, .2 ) )

rm(resE)

##### Anderson #####
# Anderson main
analyze_one_meta( dat = da,
                  meta.name = "Anderson main", 
                  pub.bias = pub.bias )

# Anderson controlled
da2 = da[ da$study_cat == "long",]  # longitudinal studies only
analyze_one_meta( dat = da2,
                  meta.name = "Anderson controlled", 
                  pub.bias = pub.bias )

##### Ferguson #####
# caveat: never managed to get number of effect sizes to match what was reported in paper,
#  despite contact with author
# but results are still quite similar to what was reported in paper

# Ferguson main
analyze_one_meta( dat = dfa,
                  meta.name = "Ferguson main", 
                  pub.bias = pub.bias )

# Ferguson controlled
analyze_one_meta( dat = dfc,
                  meta.name = "Ferguson controlled", 
                  pub.bias = pub.bias )

##### Prescott #####
# only one analysis because this meta only included longitudinal studies
# this has only 1 study that provided 2 effect sizes, so little issue of non-independence

# Prescott main and controlled
analyze_one_meta( dat = dp,
                  meta.name = "Prescott main and controlled", 
                  pub.bias = pub.bias )

##### Write Results #####

if ( write.results == TRUE ) {
  setwd(results.dir)
  
  if( pub.bias == FALSE ) suffix = "plain" else suffix = "weightr"
  
  write.csv( resE, 
             paste( "results_table_", suffix, ".csv", sep = "" ),
             row.names = FALSE )
  
  library(xtable)
  print( xtable(resE), file = paste( "results_xtable_", suffix, ".txt", sep = "" ) )
  
}


############################### VIDEO GAME SUPPLEMENTARY ANALYSES: PROPORTION BELOW ############################### 

# just define a new list of threshold values
pub.bias = FALSE
ql = as.list( c( -.1, -.2 ) )

rm(resE)

##### Anderson #####
# Anderson main
analyze_one_meta( dat = da,
                  meta.name = "Anderson main", 
                  pub.bias = pub.bias,
                  tail = "below" )

# Anderson controlled
da2 = da[ da$study_cat == "long",]  # longitudinal studies only
analyze_one_meta( dat = da2,
                  meta.name = "Anderson controlled", 
                  pub.bias = pub.bias,
                  tail = "below" )

##### Ferguson #####
# caveat: never managed to get number of effect sizes to match what was reported in paper,
#  despite contact with author
# but results are still quite similar to what was reported in paper

# Ferguson main
analyze_one_meta( dat = dfa,
                  meta.name = "Ferguson main", 
                  pub.bias = pub.bias,
                  tail = "below" )

# Ferguson controlled
analyze_one_meta( dat = dfc,
                  meta.name = "Ferguson controlled", 
                  pub.bias = pub.bias,
                  tail = "below" )

##### Prescott #####
# only one analysis because this meta only included longitudinal studies
# this has only 1 study that provided 2 effect sizes, so little issue of non-independence

# Prescott main and controlled
analyze_one_meta( dat = dp,
                  meta.name = "Prescott main and controlled", 
                  pub.bias = pub.bias,
                  tail = "below" )

##### Write Results #####

if ( write.results == TRUE ) {
  setwd(results.dir)
  
  if( pub.bias == FALSE ) suffix = "plain" else suffix = "weightr"
  
  write.csv( resE, 
             paste( "results_table_plain_propbelow.csv", sep = "" ),
             row.names = FALSE )
  
  library(xtable)
  print( xtable(resE), file = paste( "results_xtable_plain_propbelow.txt", sep = "" ) )
  
}


