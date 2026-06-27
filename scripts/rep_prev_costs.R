####
#Generating Prevalence Replication Table
####
packages <- c(
  "lfe",
  "xtable",
  "dplyr",
  "psych",
  "sjmisc",
  "locfit",
  "lme4",
  "purrr"
)

# Install missing packages
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))


###############
#Prevalence Annual Rates per 100.000
###############


###Data Wrangling and Loading

cy = read.csv("../source_data/cy.csv",as.is = T) # country-year frame with population and gdp variables.
load("../source_data/imp_results_20231110.RData") # load datasets with: (1) country rate variables; (2) country rates with imputed means; (3) sd's of imputed values in c.rates.
d2e = read.csv("../reproduction_tables/gdp_cntryr_cfloss.csv",as.is = T) #gdp loss estimates 
tang.cn = read.csv("../source_data/tangible.costs.cn.csv",as.is = T) # tangible cost estimates for both cv and ip violence

d = cr5[[1]]  # violence data with NAs ("raw")
d.imp = cr5[[2]]  # violence data with imputations (no NAs)

d$region = cy$region[match(d$ccode,cy$ccode)]
d.imp$region = cy$region[match(d$ccode,cy$ccode)]
d$pop19 = cy$unpop[cy$year==2019][match(d$ccode,cy$ccode[cy$year==2019])]
d.imp$pop19 = cy$unpop[cy$year==2019][match(d$ccode,cy$ccode[cy$year==2019])]

# Bring in npsv12 estimate, and terrorism.alone.meltt and .tac estimates.
npsv12 = read.csv("../source_data/npsv12_august_2023.csv",as.is=T)
d.imp$npsv12.est = npsv12$nspv12.max.hat[match(d.imp$ccode,npsv12$ccode)]
d$npsv12.est = npsv12$nspv12.max.hat[match(d$ccode,npsv12$ccode)]  
d$npsv12.est[d$ccode %in% npsv12$ccode[is.na(npsv12$npsv12.max)]] = NA
# setting to NA in unimputed data for the raw npsv12 estimates.  (Ie the raw data for 89 countries (DHS + OSCE + FRA + ICVS))

d$assault = d$gal.assault
d.imp$assault = d.imp$gal.assault

# SUBTRACT NPSV12 FROM ASSAULT
d$assault = with(d,assault - npsv12.est)
d.imp$assault = with(d.imp,assault - npsv12.est,na.rm=T)

# MAKE FINAL CV VARIABLE FORMS
# need to bring in c.rates because don't currently have terrorism.alone.tac in imp data.
c.rates = read.csv("../source_data/c.rates.csv",as.is = T)
d$terrorism.alone.tac = c.rates$terrorism.alone.tac[match(d$ccode,c.rates$ccode)]
d.imp$terrorism.alone.tac = c.rates$terrorism.alone.tac[match(d.imp$ccode,c.rates$ccode)]

d$colvio.meltt = with(d,ged.high + terrorism.alone.meltt)
d$colvio.tac = with(d,ged.high + terrorism.alone.tac)
d.imp$colvio.meltt = with(d.imp,ged.high + terrorism.alone.meltt)
d.imp$colvio.tac = with(d.imp,ged.high + terrorism.alone.tac)

# variables
viovars = c('colvio.meltt','colvio.tac','ged.high','terrorism.alone.meltt','terrorism.alone.tac','hom.max','who.su','assault','npsv12.est','vaw12.max','ipv12.max','sexvio12.max','ipvsexvio12.max','child.s.pun.max','fgm.sum','cm.female14','cm.female15.17','cm.male10.17')

othviovars = c('gal.assault.fem','child.viodisc.max','child.physpun.max','sexipvlife.un15','vaw.phys.life.un17',
               'sexviolife','ipvsexvio.life.un17','ipvlife.who18')

covars = c('country','ccode','region','pop19','pop.mean','median.age','gdpe.pc',
           'gws.tp','litres.pc','civ.gun.count','ggap','gii','ethfrac','polity2','wmpri','wmsec',
           'wnonag','fert','femagemp.ilo','femind.ilo','beatok','femlaborprieduc','femlaborseceduc','fuelexpsh','gini.disp')



##########################################
###Calculate Costs 

discdat = read.csv("../source_data/ip_cost_discount_factors.csv",as.is = T)#load discount rates
# need to bring in imputation estimates for ip vio at country level.
load("../source_data/imp_results_20231110.RData") # this loads cr5, which is a list of three data frames: (1) c.rates variables; (2) c.rates with imputed means; (3) sd's of imputed values in c.rates.
prev = cr5[[1]]  # "raw" prevalence estimates, with NAs, for 197 countries.
prev.imp = cr5[[2]] # imputed prevalence estimates.

prev$region = cy$region[match(prev$ccode,cy$ccode)]
prev.imp$region = cy$region[match(prev.imp$ccode,cy$ccode)]
prev$pop19 = cy$unpop[cy$year==2019][match(prev$ccode,cy$ccode[cy$year==2019])]
prev.imp$pop19 = cy$unpop[cy$year==2019][match(prev.imp$ccode,cy$ccode[cy$year==2019])]

### SET US VSL and produce other countries' 2019 VSLs from this.
usvsl = 10.9E6 # This is US Dept of Transportation VSL for 2019.
# compute global average (pop weighted) VSL for 2019, using gdpe
t = with(d2e[d2e$year==2019,],sum(gdpe.pc*pop/sum(pop,na.rm=T),na.rm=T))  # global per capita PWT income (expenditure) in 2019
usvsl*t/d2e$gdpe.pc[d2e$ccode==2 & d2e$year==2019]  # 2.81 million = global avg vsl

# need to compute vsle for all countries, not just conflict-affected.
vsls = d2e[d2e$year == 2019,c('country','ccode','gdpe.pc')] # a bunch of NAs here though.
vsls[is.na(vsls$gdpe.pc),]  # good. missing only some microstates and Venezuela and North Korea. For Venezuela use 10k?
vsls$gdpe.pc[vsls$ccode==101] = 10000
# North Korea:  CIA Factbook gives a Maddison PPP estiamte of $40 billion in 2015 (2015 ppp dollars).
# TradingEconomics.com say 18 billion USD in 2019.  Not PPP. Shows growth estimate of 18/16.28 = 1.106 from 2015, for agg GDP.
# unpop growth for this period is 25755.44/25258.01 = 1.019694 so per cap growth would be 1.106/1.019 = 1.085378.
# so if use 2015 maddison PPP estimate would have 1.085378*(40E9/25.25801E6) = 1719 for 2018. Which seems awfully small.
vsls$gdpe.pc[vsls$ccode==731] = 1719

# Choice of US VSL number for all that follows.
# using usvsl from above. (DOT number)
vsls$vsle = with(vsls,usvsl*gdpe.pc/gdpe.pc[ccode==2])  
prev.imp$vsle = vsls$vsle[match(prev.imp$ccode,vsls$ccode)]  

###Data Wrangling
# Bring in npsv12 estimate.
# npsv12.max estimates updated following import and use of full ICVS estimates, August 2023. See npsv_estimation.R for the approach.
npsv12 = read.csv("../source_data/npsv12_august_2023.csv",as.is=T)
prev.imp$npsv12.est = npsv12$nspv12.max.hat[match(prev.imp$ccode,npsv12$ccode)]  # misspelling in var name.
# SUBTRACT NPSV12 FROM gal.assault to make assault used in prevalence, and apply "discount"
# Using .35 based on evidence that "robbery" category quite frequently involves (what we call) physical assault.  Medical costs
# are usually higher than for assault, for example, in studies that have this.
prev.imp$assault = .35*(prev.imp$gal.assault - prev.imp$npsv12.est)

# display variables (from Table 3.1 in chapter 3)
dispvars = c('ged.high','terrorism.alone.meltt','hom.max','assault','ipv12.max','ipvsexvio12.max','npsv12.est',
             'child.s.pun.max','fgm.sum','cm.female14','cm.female15.17','who.su')
com3 = prev.imp[,c('country','ccode','region','pop19','vsle',dispvars)] # keeping com3, com4 names because used below.
# make a version with counts, not rates
#com3 is rates per 100.000
com4 = cbind(com3[,1:5],apply(com3[,6:17],2,function(x) x*com3$pop19/1E5)) # thus counts are in 1000s here.

#INTANGIBLE COSTS
#make data frame with share of VSL cost for each type of ip violence.
#convert to version matching com4 etc violence type names (which are what we use in graphs and tables and text)
dth.share = data.frame(vio = names(com4[6:17]))
dth.share = cbind(dth.share,discdat[match(dth.share$vio,discdat$ip.type),])
dth.share[1:2,4:5] = 1 # for colvio direct deaths in conflict and terrorism.

# Intangible cost estimates by country, with country specific vsl, for PWT rgdpe. dividing by 1000 to put in $millions. (because com4 counts are in 1000s)
ip.intang.costs.rgdpe = cbind(com4[,c(1:5)],sapply(names(com4)[6:17],function(x) .001*com4[,x]*com4$vsle*dth.share$dth.share.intang[(1:12)[x==names(com4)[6:17]]]))

# Intangible cost estimates by country, with global average VSL, for PWT rgdpna and PWT rgdpe
vsle.gl = with(com4,sum(vsle*pop19/sum(pop19,na.rm=T),na.rm=T))

ip.intang.costs.rgdpe.gl = cbind(com4[,c(1:5)],sapply(names(com4)[6:17],function(x) .001*com4[,x]*vsle.gl*dth.share$dth.share.intang[(1:12)[x==names(com4)[6:17]]]))

# to get intangible ip costs as share of gdpe, only need to divide all by aggregate gdpe in 2019.

# Tangible IP and CV costs.
# tang, and tang.cn already have discounted sums according to 
# a common discount rate drate.  Using drate = 0 except for robustness checks, to make comparable to treatment of 
# economic growth losses except in long fun (we only use social discount rate for post 2019, and
# these tangible cost estimates are based on 2010 averages.
 
tang = read.csv("../source_data/tangible.costs.cy.csv",as.is = T)
# reconfigure tang.cn so that it can be used in functions the same as ip.intang.costs.rgdpe
ip.tang.costs = ip.intang.costs.rgdpe[,1:5,]  # note this has the 197 measured countries. Note also we won't use vsle because tangible cost estimates never use this. But need for iptab() to work.
# rename variables to make ip.intang.costs.rgdpe
tang.cn$ged.high = tang.cn$tang.ged.costs
tang.cn$terrorism.alone.meltt = tang.cn$tang.terrorism.costs
tang.cn$hom.max = tang.cn$tang.hom.costs
tang.cn$who.su = tang.cn$tang.sui.costs

ip.tang.costs = cbind(ip.tang.costs,tang.cn[match(ip.tang.costs$ccode,tang.cn$ccode),names(ip.intang.costs.rgdpe)[6:17]])
# divide by 1E6 to put in millions, parallel to intangible costs.
ip.tang.costs[,6:17] = apply(ip.tang.costs[,6:17],2,function(x) x/1E6)

#############
##Prevalence
#############

prev_cn_aggr <- com3
prev_cn_aggr$vsle <- NULL
prev_cn_aggr$pop19 <- NULL
prev_cn_aggr$colvio <-prev_cn_aggr$ged.high + prev_cn_aggr$terrorism.alone.meltt

#rename variables
cols_prev <- c("armed.conflict", "terrorism", "hom", "assault", "intimpar.vio", "intimpar.sex.vio", 
               "nonpar.sex.vio", "child.pun", "fgm", "cm.girls.14", "cm.girls.15", "sui")
colnames(prev_cn_aggr)[4:15] <- cols_prev
write.csv(prev_cn_aggr, "../reproduction_tables/prevalence_cntr_averages.csv")


#############
## VSL Costs
#############

#Rename cols and drop population
mod_vsl <- function(df, suffix) {
  
  if ("pop19" %in% names(df)) {
    df$pop19 <- NULL
  }
  
  start_pos <- match("ged.high", names(df))
  
  if (is.na(start_pos)) {
    stop('"ged.high" not found in dataset')
  }
  
  end_pos <- start_pos + length(cols_prev) - 1
  
  if (end_pos > ncol(df)) {
    stop("Not enough columns after 'ged.high' to rename 12 columns")
  }
  
  names(df)[start_pos:end_pos] <- cols_prev
  colnames(df)[start_pos:end_pos] <- paste0(colnames(df)[start_pos:end_pos], suffix)
  
  return(df)
}


#Generate tables with VSL costs
ip.intang.costs.rgdpe <- mod_vsl(ip.intang.costs.rgdpe, ".intang.vsl.cn")
ip.intang.costs.rgdpe.gl <- mod_vsl(ip.intang.costs.rgdpe.gl, ".intang.vsl.gl")
ip.tang.costs <- mod_vsl(ip.tang.costs, ".tang.vsl.cn")

vsl_cost_cntr <- merge(ip.intang.costs.rgdpe, ip.intang.costs.rgdpe.gl, by = c("country", "ccode", "region", "vsle"))
vsl_cost_cntr <- merge(vsl_cost_cntr, ip.tang.costs, by = c("country", "ccode", "region", "vsle"))

write.csv(vsl_cost_cntr, "../reproduction_tables/vsl_costs_cntr.csv")