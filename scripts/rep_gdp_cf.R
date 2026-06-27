
####
#Replicating Generation of Counterfactual GDP Values and GDP Loss
####
packages <- c(
  "lfe",
  "xtable",
  "psych",
  "sjmisc",
  "locfit",
  "lme4",
  "plyr",
  "purrr",
  "dplyr"
)

# Install missing packages
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

d <- read.csv("../source_data/gdp_ests.csv")

########################################
#Estimate counterfactual GDPs & GDP Loss
####

# lag function
lg <- function(x,panelvar = d$ccode, first=FALSE) {
  y <- c(NA,x[1:(length(x)-1)])
  cm <- panelvar != c(0,panelvar[1:(length(panelvar)-1)])	
  y[cm] <- NA
  if (first==TRUE) { y[cm] <- x[cm] }
  y}

# (1) use median regional peaceful country growth rate in year t as cf for treatment years of country i
# (2) use 75th percentile gy of all years of country i as cf for treatment years. ("optimistic")
# (3) Mixed model that estimates linear effect of log(colvio.pc > 2.9), and is zero below this. Then set to zero to compute cf.
# (4) Mixed model which runs on peaceful years for each country, predicts for treatment years.

# Method 1:  Median growth non treatment countries in same region, for each year.  
t = ddply(d[d$treatment==0 & d$year > 1979,],.(region,year),summarize, 
          gye.reg = median(gye,na.rm=T),n = sum(!is.na(gye)))
# using median because better measure of central tendency for variables with erratic tails and possibly low N.

d$gye.reg = t$gye.reg[match(paste(d$region,d$year),paste(t$region,t$year))]

# now substitute region peaceful control estimate in treatment years.
d$cfgye.reg = d$gye.reg  # need to start with gye.reg and replace non cf year or get NAs in cfgye.reg when shouldn't.
d$cfgye.reg[d$treatment==0 & d$recovery==0 & d$year > 1988 & d$year < 2020] = d$gye[d$treatment==0 & d$recovery==0 & d$year > 1988 & d$year < 2020]
# note that using regional mean reduces annual volatility a lot, which tends to favor growth.

# METHOD 2: optimistic cf.  Take X pctile from country years after 1999, dropping recovery years because
# they often have crazy high growth rates. Leave cw = 1 years in, because some countries have very few 
# non-cw years, and some countries have stronger growth during cw than in the years outside it.

pctile = .55
t = ddply(d[d$recovery==0,],.(ccode),summarize,country = country[1],
          gye.pctile = quantile(gye[year > 1999],probs = pctile,na.rm=T),
          gye.trmedian = median(gye[year > 1999 & treatment==1],na.rm = T),
          n.treatment = sum(treatment[year > 1999],na.rm=T))
a = t[order(-t$gye.pctile),]
d$cfgye.pctile = d$gye

# MODEL 3: Mixed model cf estimates, estimating country specific intercepts and slopes for log(colvio.pc > 2.9)
# restrict to the 197 countries in our main analysis, to avoid excess influence of microstates on big country estimates.
cns = c(unique(d$ccode[d$territory==0]),347,996,915)

d$lcv = log(d$colvio.pc/2.9) #treshold for collective violence
d$lcv[d$colvio.pc < 2.9] = 0
d$lcv[is.na(d$colvio.pc)] = NA

t = lmer(gye ~ (1 + lcv|ccode) + factor(region)*factor(year),data = d,subset = ccode %in% cns & year > 1998 & year < 2020)

# what we want to do is compute for each c-year, y^{cf}_it = y_it - beta_lcv_i*lcv_it
#dhat = d[d$ccode %in% cns & d$year > 1989 & d$year < 2020,] 
dhat = d[d$ccode %in% cns & d$year > 1998 & d$year < 2020 & !d$ccode %in% c(223,731),]  # zero gye data for Liechtenstein and NKorea.
dhat$lcv = 0; dhat$recovery = 0
z = predict(t,newdata = dhat)  #This is y_hat_it under counterfactual of cv below 2.9 = 0
z.resid = resid(t) # add to prediction with newdata to get counterfactual. These are the other determinants of gye_it 
# than cv in that year

dhat$gye.hat = z # predicted value under cf of lcv = 0 and recovery = 0.
d$gye.hat = NA
d$gye.hat = dhat$gye.hat[match(rownames(d),rownames(dhat))]
# this is the predicted value.  add resid to get cf value.
d$cfgye.mm1 = d$gye.hat + z.resid[match(rownames(d),names(z.resid))]
# this automatically gets right cf for us in zero conflict years.
d$resid = z.resid[match(rownames(d),names(z.resid))]
q = predict(t)
d$gye.hat.notcf = q[match(rownames(d),names(q))]

# METHOD 4, use mixed model to predict non-treatment year gye, then predict from this for treatment years.
# avg gy of regional states other than i in year t. Used as a predictor for country it gy in treatment == 0 years.
d$gyereg.noti = NA
for (i in 1:nrow(d)) {j = d$region==d$region[i] & d$ccode != d$ccode[i] & d$year==d$year[i]
d$gyereg.noti[i] = mean(d$gye[j],na.rm=T)}

## MIXED MODEL APPROACH. Peaceful years only, random effects for country intercept and slope of effect of regional growth (other than i), plus fixed effects for years.
t = lmer(gye ~ (1 + log(gyereg.noti)|ccode) + I(log(gyereg.noti)) + factor(year),data = d,subset = year > 1998 & treatment==0)

dhat = d[d$year > 1998 & d$year < 2020 & d$ccode %in% cns & !d$ccode %in% c(223, 626, 700,731),] 
z = predict(t,newdata = dhat)
d$cfgye.mm2 = z[match(rownames(d),names(z))]
d$cfgye.mm2[d$treatment == 0 & d$year > 1998 & d$year < 2020] = d$gye[d$treatment==0 & d$year > 1998 & d$year < 2020]

# these are all fairly different approaches to generating the counterfactual.
# reg and mm2 are most similar, because mm2 uses reg growth rates to predict cf treatment growth.

# Next, total growth rate for treatment years by method, for 2000-2019
z = ddply(d[d$year > 1999,],.(ccode),country= country[1], treatment = sum(treatment,na.rm=T),summarize,
          gye = prod(gye,na.rm=T),cfgye.reg = prod(cfgye.reg,na.rm=T),
          cfgye.pctile = prod(cfgye.pctile,na.rm = T),
          cfgye.mm1 = prod(cfgye.mm1,na.rm=T),cfgye.mm2 = prod(cfgye.mm2,na.rm=T))
z[z$treatment > 0,]
z[,4:8] = apply(z[,4:8],2,function(x) ifelse(x==1,NA,x))

# for good measure, equivalent OLS ("fixed effects") prediction model
# Run model on treatment = 0 years, and then generate predictions on the full data (to avoid problems with prediction from 
# rank-deficient matrix)

q = lm(gye ~ log(gyereg.noti)*factor(ccode) + factor(year),data = d,subset = year > 1989 & treatment == 0,na.action = na.exclude)
cns2 = unique(q$model$'factor(ccode)')
dhat = d[d$ccode %in% cns2 & d$year > 1989 & d$year < 2020,c('country','ccode','year','gyereg.noti')]
z2 = predict(q,newdata = dhat)
d$cfgye.ols = z2[match(rownames(d),names(z2))]   ###
d$cfgye.ols[d$treatment == 0 & d$year > 1988 & d$year < 2020] = d$gye[d$treatment==0 & d$year > 1988  & d$year < 2020]

# total growth over period, again
z = ddply(d[d$year > 1999 & d$year < 2020,],.(ccode),country= country[1], treatment = sum(treatment,na.rm=T),summarize,
          gye = prod(gye,na.rm=T),cfgye.reg = prod(cfgye.reg,na.rm=T),
          cfgye.pctile = prod(cfgye.pctile,na.rm = T),
          cfgye.mm1 = prod(cfgye.mm1,na.rm=T),cfgye.mm2 = prod(cfgye.mm2,na.rm=T),
          cfgye.ols = prod(cfgye.ols,na.rm=T))
z[z$treatment > 0,]
z[,4:9] = apply(z[,4:9],2,function(x) ifelse(x==1,NA,x))
summary(z[z$treatment > 0,])
# gets some crazy estimates, like Somalia


###  Need to deal with NAs for some methods.  But only for treated  cyears.
t = with(d[d$year > 1999 & d$year < 2020,],tapply(treatment,ccode,function(x) as.numeric(sum(x,na.rm=T))))
tr.cns = as.numeric(names(t[t > 0]))

summary(d[d$year > 1998 & d$year < 2020 & d$ccode %in% tr.cns,c('cfgye.reg','cfgye.pctile','cfgye.mm1','cfgye.mm2')])
# only mm2.
d[d$year > 1998 & d$year < 2020 & d$ccode %in% tr.cns & is.na(d$cfgye.mm2),c('country','ccode','year','gdpe','gye','treatment','cfgye.reg','cfgye.pctile','cfgye.mm1','cfgye.mm2')]  # all Sudan, all Afghanistan, and one year of Palestine.
# reason for Sudan and Afghanistan for mm2 is that they are treated in all years, so that we can't estimate relation between
# regional growth rate and country growth rate in peace years.
# therefore, natural and best to just use regional peaceful growth for these cases, to avoid dropping (and so
# not counting) these cases.
d$cfgye.mm2[d$ccode == 626] = d$cfgye.reg[d$ccode==626]
d$cfgye.mm2[d$ccode == 700] = d$cfgye.reg[d$ccode==700]


#### COST ESTIMATES

# produce estimates for 2000-2019 using these four methods in one table, then a table with rule that sets negative costs to zero.
# write function that (1) produces cfgdpe path and cost estimates for a given country and cfgye series, from a start year;
# (2) alternative estimate that goes year by year and applies cfgye for max of cf and actual gdpe in that year. 
# (2) is an optimistic approach that assumes could never have cf income below actual (which is probably wrong from some cases,
# like Afghanistan.)

# convenient to create lagged gdpe
d$gdpel = lg(d$gdpe)
# fix missing lag for South Sudan 2010 (see mk_gdp.R for the imputation, which produced a gdpe.imp est for S Sudan 2010 there.)
d$gdpel[d$ccode==626 & d$year==2011] = 5859.5
d[d$ccode==626,c('country','year','gdpe','gdpel','gye')]

# Function to compute for ccode country cf gdpe in year based on cfgye estimates in case (1) allow actual gye > cfgye,
# and (2) assume actual gye <= cfgye.
# this function only works for countries with some gdpe data, and may for treatment = 0 countries in some cases. Below, 
# only applied to treatment > 0 countries but we don't need to compute cf's for non-treatment countries.

cf = function(ccode,cfgye = "cfgye.reg",styear = 2000) {
  a = d[d$ccode == ccode & d$year >= styear & d$year < 2020,c('country','ccode','year','treatment','recovery','gdpe','gdpel','gye',cfgye)]
  a$cfgye = a[,cfgye]
  #a$cfgye[is.na(a$cfgye)] = 1 # good way to deal with missing obs at start of some series for countries that enter after 2000.
  if (any(!is.na(a$gye))) {styr = min(a$year[!is.na(a$gye)],na.rm=T)}  # some countries missing these, messes up series.
  a$cfgdpe[a$year >= styr] = a$gdpel[a$year==styr]*cumprod(a$cfgye[a$year >= styr])
  # now do method (2).
  a$cfgdpe2 = a$gdpe; a$cfgdpe2[1] = with(a,max(cfgdpe2[1],gdpel[1]*cfgye[1],na.rm=T))
  for (i in 2:nrow(a)) {if ((a$treatment[i]==1|a$recovery[i]==1)  & !is.na(a$cfgye[i])) {  # Note:  dropping & !is.na(a$cfgdpe2[i]) from the condition, so that we impute from cfgye series when there is missing gdpe data, as for Somalia, Eritrea, South Sudan.
    a$cfgdpe2[i] = max(a$gdpe[i],a$cfgdpe2[i-1]*a$cfgye[i],na.rm=T)}
    if ((a$treatment[i]==0 & a$recovery[i]==0) & !is.na(a$cfgdpe[i]) & !is.na(a$cfgye[i])) {
      a$cfgdpe2[i] = a$cfgdpe2[i-1]*a$gye[i]}  # using observed growth rate for non-treatment years}
  }
  a}

# make for all treated countries. to avoid problems with countries have no data at all, do for countries at least one treatment year,
# and set others to cfgdpe.reg = gdpe.  Actually, don't even need these to calculate costs.
t = with(d[d$year > 1999 & d$year < 2020,],tapply(treatment,ccode,function(x) as.numeric(sum(x,na.rm=T))))
tr.cns = as.numeric(names(t[t > 0]))

cfgye.meth = c('cfgye.reg','cfgye.pctile','cfgye.mm1','cfgye.mm2')
meth1 = c('cfgdpe.reg','cfgdpe.pctile','cfgdpe.mm1','cfgdpe.mm2')
meth2 = c('cfgdpe2.reg','cfgdpe2.pctile','cfgdpe2.mm1','cfgdpe2.mm2')

d.cfgdps = d[,c('country','ccode','year')]  # will put variables here and then move to d.

for (k in cfgye.meth) {num = (1:4)[k==cfgye.meth]
for (i in tr.cns) {j = d.cfgdps$ccode == i & d.cfgdps$year > 1999 & d.cfgdps$year < 2020
z = cf(i,k)
d.cfgdps[j,meth1[num]] = z$cfgdpe
d.cfgdps[j,meth2[num]] = z$cfgdpe2
}}

d = cbind(d,d.cfgdps[,c(meth1,meth2)]) 
rm(d.cfgdps)

# initial cost estimates.  using only countries with treatment years.
d2 = d[d$ccode %in% tr.cns & d$year > 1999 & d$year < 2020,c('country','ccode','year','territory','region','unpop',
                                                             'colvio','colvio.pc','treatment','recovery','gdpe','gye',meth1,meth2)]
t = apply(d2[,c(meth1,meth2)],2,function(x) round(x - d2$gdpe))
colnames(t) = paste0("diff.",colnames(t))
d2 = cbind(d2,t)

summary(d2)
summary(d2[d2$treatment==1,])

# for later use:
ip.cns = unique(d$ccode[d$territory==0|(d$ccode %in% c(915,996))])  # UN member states plus Palestine and Hong Kong (ip violence countries)

# The concept here: We estimate the economic costs of conflict in a period from year x to year y by starting from gdp/cap
# in year x, and computing growth from that point on using cfgy's for treatment years in the period in question.
# compute difference in gdp and cfgdp for each country and divide by years in period (y - x + 1) gives an estimate 
# of the annual average costs of conflict in that period.  Set to zero years in which cfgdp < gdp, on ground that these
# are prima facie underestimates of costs, for most costs.

# First prepare data sets with gdp and gdpe base PWT numbers. d2e = for gdpe.
# GDPe
d2e = d[d$year > 1999 & d$year < 2020,c('country','ccode','year','region','colvio.pc','treatment','recovery','unpop','gdpe','gye',meth1,meth2)]
summary(d2e)
summary(d2e[d2e$treatment==1,])

# compute total estimated costs (not per capita) for each country year (in millions)
t = with(d2e,round((d2e[,c(meth1,meth2)] - d2e$gdpe)*d2e$unpop/1E3)) 
summary(t[,])
# add a column with aggregate GDP in millions.
t = cbind(t,d2e$gdpe*d2e$unpop/1E3)
names(t)[9] = "GDPe"

d2e[t$cfgdpe.pctile < -290000 & !is.na(t$cfgdpe.pctile),]
d2e[t$cfgdpe.mm2 < -290000 & !is.na(t$cfgdpe.mm2),]
# some Russia years, and Iraq 2012, 2013.

# add loss estimates to d2e
d2e[,c('loss.reg','loss.pctile','loss.mm1','loss.mm2','loss2.reg','loss2.pctile','loss2.mm1','loss2.mm2','GDPe')] = t 

# Approach: compute pdv costs for each country separately, along with VSL-mortality equivalent using VSL income elasticity.
# typical global estimate is 1.0.
# need to do this at per capita level or else need to deal with what to assume about population growth.

# treated countries
d2e$treated.cn = as.numeric(d2e$ccode %in% tr.cns)
colnames(d2e)[colnames(d2e) == "gdpe"] <- "gdpe.pc"
colnames(d2e)[colnames(d2e) == "unpop"] <- "pop"
colnames(d2e) <- gsub("2\\.", ".scn2.", colnames(d2e))
colnames(d2e)
write.csv(d2e,"../reproduction_tables/gdp_cntryr_cfloss.csv")



##########################################
#Calculate GDP Country-Year Loss Estimates
####


#Load data
d2e = read.csv("../reproduction_tables/gdp_cntryr_cfloss.csv", row.names = 1)
cy = read.csv("../source_data/cy.csv",as.is = T) # country-year frame with population and gdp variables.

# function to create country level econ loss estimates for given GDP method . Both amounts and shares of country gdp
# note: this produces estimates only for treated countries (at least 1 conflict year)
c.loss = function(var = "cfgdpe.reg", delta = .02, eta = 1.4) {
  a = d2e[d2e$treated.cn==1,]; a$cfy = a[,var] # cf gdp/capita.
  a$loss = with(a,cfy - gdpe.pc)
  a$delta = delta
  a$eta = eta
  
  t = ddply(a,.(ccode),summarize,country = country[1],region = region[1],pop19 = pop[year==2019],years = length(ccode),
            loss.agg = mean(loss*pop/1E3,na.rm=T),  # putting aggregate loss in millions.
            loss.pc = mean(loss,na.rm = T),
            loss.pc.sh = mean(loss/gdpe.pc,na.rm=T),  # loss as share of gdpe.pc
            gdpgap = cfy[year==2019] - gdpe.pc[year==2019],
            styr = min(year), # first year in data
            cfyst = cfy[year==styr],
            fut.gdpgrw = (cfy[year==2019]/cfyst)^(1/(2020-styr)) - 1, # this allows negative cf growth.
            fut.gdpgrw = max(1,(cfy[year==2019]/cfyst)^(1/(2020-styr))) - 1, # this forbids negative cf growth.
            sdr = (1 + fut.gdpgrw)/(1 + delta[1] + eta[1]*fut.gdpgrw),  # social discount rate for future gdp per capita
            fut.gdp.loss.pc = gdpgap/(1 - sdr)/years, # dividing by # years to get avg annual loss
            fut.gdp.loss.agg = fut.gdp.loss.pc*pop19/1E3, # putting in millions of dollars
            fut.gdp.loss.pc.sh = fut.gdp.loss.pc/gdpe.pc[year==2019]
  )
  
  t
}



# TABLES
# code to make presentable cost tables showing regional and world estimates, using output of c.loss
etab = function(pdv.cn.data = t) {
  a = ddply(pdv.cn.data,.(region),summarize,loss = sum(loss.agg,na.rm=T),fut.loss = sum(fut.gdp.loss.agg,na.rm=T)) # in millions
  world = apply(a[,2:3],2,sum)
  a = round(rbind(a[,2:3],world))
  rownames(a) = c('Asia','EEur','LA/Ca','NA/ME','SSA','World')
  a}

etab(c.loss())
etab(c.loss("cfgdpe.reg", 0.2, 1.4))

meth1 = c('cfgdpe.reg','cfgdpe.pctile','cfgdpe.mm1','cfgdpe.mm2')
meth2 = c('cfgdpe.scn2.reg','cfgdpe.scn2.pctile','cfgdpe.scn2.mm1','cfgdpe.scn2.mm2')

etabs = lapply(c(meth1, meth2),function(x) etab(c.loss(x)))

# now for share of per capita gdpe. For this, we need to include all the zero (non treated countries).
cns = c(unique(cy$ccode[cy$territory==0]),915,996)

etab.sh = function(pdv.cn.data) {a = d2e[d2e$ccode%in% cns & d2e$year==2019,] 
  a$pop19 = d2e$pop[d2e$year==2019][match(a$ccode,d2e$ccode[d2e$year==2019])]
  a = cbind(a,pdv.cn.data[match(a$ccode,pdv.cn.data$ccode),c('loss.pc.sh','fut.gdp.loss.pc.sh')])
  a$loss.pc.sh[is.na(a$loss.pc.sh)] = 0; a$fut.gdp.loss.pc.sh[is.na(a$fut.gdp.loss.pc.sh)] = 0
  b = ddply(a,.(region),summarize,loss = sum(loss.pc.sh*pop19/(sum(pop19,na.rm=T)),na.rm=T), # in millions
          fut.loss = sum(fut.gdp.loss.pc.sh*pop19/sum(pop19,na.rm=T),na.rm=T))
  world.loss.pc.sh = with(a,sum(loss.pc.sh*pop19,na.rm=T)/sum(pop19,na.rm=T))
  world.fut.gdp.loss.pc.sh = with(a,sum(fut.gdp.loss.pc.sh*pop19,na.rm=T)/sum(pop19,na.rm=T))
  b = rbind(b,c(NA,world.loss.pc.sh,world.fut.gdp.loss.pc.sh))
  b = b[,c(2,3)]
  rownames(b) = c('Asia','EEur','LA/Ca','NA/ME','SSA','West','World')
  b}

etab.sh(c.loss("cfgdpe.reg"))
etab.sh(c.loss("cfgdpe.scn2.pctile"))

etabs.sh = lapply(c(meth1,meth2),function(x) etab.sh(c.loss(x)))


###########################################
#Calculate GDP Total Country Loss Estimates
#############

make_cfgdpe <- function(df_name, suffix) {
  df <- c.loss(df_name)
  colnames(df)[6:ncol(df)] <- paste0(colnames(df)[6:ncol(df)], suffix)
  df
}

cfgdpe.reg  <- make_cfgdpe("cfgdpe.reg",       ".reg")
cfgdpe.pctl <- make_cfgdpe("cfgdpe.pctile",    ".pctl")
cfgdpe.mm1  <- make_cfgdpe("cfgdpe.mm1",       ".mm1")
cfgdpe.mm2  <- make_cfgdpe("cfgdpe.mm2",       ".mm2")

cfgdpe2.reg  <- make_cfgdpe("cfgdpe.scn2.reg",    ".scn2.reg")
cfgdpe2.pctl <- make_cfgdpe("cfgdpe.scn2.pctile", ".scn2.pctl")
cfgdpe2.mm1  <- make_cfgdpe("cfgdpe.scn2.mm1",    ".scn2.mm1")
cfgdpe2.mm2  <- make_cfgdpe("cfgdpe.scn2.mm2",    ".scn2.mm2")


df_list <- list(
  cfgdpe.reg,
  cfgdpe.pctl,
  cfgdpe.mm1,
  cfgdpe.mm2,
  cfgdpe2.reg,
  cfgdpe2.pctl,
  cfgdpe2.mm1,
  cfgdpe2.mm2
)
by_vars <- c("country", "ccode", "region", "pop19", "years")
gdp_loss_tot <- reduce(df_list, full_join, by = by_vars)

write.csv(gdp_loss_tot, "../reproduction_tables/gdp_cntr_totloss.csv")

