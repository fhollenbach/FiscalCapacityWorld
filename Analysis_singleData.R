#####################################################
# Fiscal Capacity and Inequality #
#####################################################
rm(list = ls())
options(stringsAsFactors = FALSE)

### Libraries
library(ggplot2)
library(stargazer)
library(haven)
library(lfe)
library(dplyr)
library(lmtest)
library(sandwich)
library(rstanarm)

## Clutered-SE Function 
robustse.f <- function(model, cluster, df_correction) {
  ## Huber-White heteroskedasticity-robust standard error calculation and table generation code for lm and glm models in R.
  ##Written by Joshua Gubler ~  http://joshuagubler.com.  Note that the second half of this function is just a wrapper for the excellent "multiwaycov" package here: https://cran.r-project.org/web/packages/multiwayvcov/multiwayvcov.pdf .  Love the authors of that package...
  ##Last updated: 3 November 2016  
  
  #model = the model you estimated, now calculated with robust standard errors
  #cluster = the name of the entity on which you will cluster (need to put full path: e.f. df$cluster).  If you don't put in a cluster, you will simply get huber-white robust errors.
  #df_correction: If you do not want the number of levels in your cluster variable to count against your degrees of freedom (like the xt- options in Stata), then type "F".  Otherwise, type "T" and these levels will be counted against your degrees of freedom
  
  require(sandwich)
  require(lmtest)
  require(multiwayvcov)
  if(missing(cluster)) {
    name <- deparse(substitute(model))
    modelname <- paste(name,"rob",sep=".")
    model$se <- coeftest(model, vcov=vcovHC(model,"HC1"))[,2]
    model$vcovHC <- vcovHC(model,"HC1")
    assign(modelname,model,envir = .GlobalEnv)
    coeftest(model, vcov=vcovHC(model,"HC1"))
  } else {
    name <- deparse(substitute(model))
    modelname <- paste(name,"rob",sep=".")
    vcovCL <- cluster.vcov(model, cluster, df_correction = df_correction)
    model$vcovCL <- vcovCL
    modelname <- paste(name,"clustrob",sep=".")
    model$se <- coeftest(model, vcovCL)[,2]
    assign(modelname,model,envir = .GlobalEnv)
    coeftest(model, vcovCL)
  }
}







#################################################
#### Load Data
#################################################
setwd("~/Dropbox/FiscalCapacityCrossNational/Data/")

load("data_mean_inequ.rda")
data <- data_mean_inequ[data_mean_inequ$Year >1989 & is.na(data_mean_inequ$Year)==F,]
data <- subset(data, is.na(Year) == F)


# control variables
data$lpop <- log(data$PopulationTotal)
data$lgdp <- log(data$GDPConst2005US)



source("~/Dropbox/FiscalCapacityCrossNational/Rcode/demean_fun.R")
x <- demean_function(data, "COWCODE", c("lpop","lgdp", "gini_disp", "tradePctGDP", "AgeDependency"))





control_list_1 <- c("lag(tradePctGDP)", "lag(lpop,1)", "lag(lgdp,1)")


control_list_2<- c("lag(e_miinterc,1)","lag(e_miinteco,1)",
                   "lag(e_Civil_War,1)", "lag(e_migdpgro,1)", "lag(e_miurbani,1)")

control_list_3 <- c("lag(AgeDependency,1)","lag(PoliticalTerrorState,1)")










x <- subset(data, Year < 2006 & Year > 2000)
cross <- ddply(x, .(COWCODE), summarize, gdp = mean(GDPConst2005US, na.rm = T), ageDep = mean(AgeDependency, na.rm = T), bureauCap = mean(BQuality, na.rm = T), tradeTaxTotalTax = mean(TradeTaxPctTotalTax, na.rm = T), TradeTaxIncomeTax = mean(TradeTaxPctIncomeTax, na.rm = T), amnesty = median(PoliticalTerrorAmnesty, na.rm = T), gini = mean(gini_disp, na.rm = T), e_miinteco = max(e_miinteco, na.rm = T), e_miinterc = max(e_miinterc, na.rm = T), polTerrorState = median(PoliticalTerrorState), aid = mean(NetODAPctGNI, na.rm = T), population = mean(PopulationTotal, na.rm = T), trade = mean(tradePctGDP, na.rm = T) )


cross$lpop <- log(cross$population)
cross$lgdp <- log(cross$gdp)
m1 = lm(tradeTaxTotalTax ~ gini + lpop + lgdp + trade, data =cross)
summary(m1)
m2 = lm(TradeTaxIncomeTax ~ gini + lpop + lgdp + trade, data =cross)
summary(m2)
m3 = lm(bureauCap ~ gini + lpop + lgdp + trade, data =cross)
summary(m2)




mod1 = plm(TradeTaxPctIncomeTax ~ lag(TradeTaxPctIncomeTax,1) +lag(gini_disp,1) + lag(lgdp,1) + lag(NetODAPctGNI,1) + lag(tradePctGDP,1) + lag(lpop,1), data=data , index=c("COWCODE", "Year"), model = "within", effect = "twoway")
summary(mod1)
mod2 = plm(TradeTaxPctTotalTax ~ lag(TradeTaxPctTotalTax,1) +lag(gini_disp,1) + lag(lgdp,1) + lag(NetODAPctGNI,1) + lag(tradePctGDP,1) + lag(lpop,1), data=data , index=c("COWCODE", "Year"), model = "within", effect = "twoway")
summary(mod2)



m1 = lm(TradeTaxPctIncomeTax ~  gini_disp_mean + lgdp_mean + lpop_demeaned + tradePctGDP_demeaned + AgeDependency_demeaned + lpop_demeaned + lgdp_demeaned + gini_disp_demeaned + tradePctGDP_demeaned + AgeDependency_demeaned + as.factor(COWCODE), data = x)
summary(m1)

#################################################
#### Transform Data
#################################################

# Cheibub democracy
data$chga_demo_dum <- ifelse(data$CheibubDemocracy=="1. Democracy",1,NA)
data$chga_demo_dum <- ifelse(data$CheibubDemocracy=="0. Dictatorship",0,data$chga_demo_dum)

# Specific year dummies
data$decades <- ifelse(data$Year>1949 &data$Year<1959,"1950s",NA)
data$decades <- ifelse(data$Year>1959 &data$Year<1969,"1960s",data$decades)
data$decades <- ifelse(data$Year>1969 &data$Year<1979,"1970s",data$decades)
data$decades <- ifelse(data$Year>1979 &data$Year<1989,"1980s",data$decades)
data$decades <- ifelse(data$Year>1989 &data$Year<1999,"1990s",data$decades)
data$decades <- ifelse(data$Year>1999 &data$Year<2010,"2000s",data$decades)
data$decades <- ifelse(data$Year>2010 ,"2010s",data$decades)





# control variables
data$lpop <- log(data$PopulationTotal)
data$lgdp <- log(data$GDPConst2005US)
# Define variables for analysis
data$interaction <- data$gini_disp*data$democracy
treat_list <- c("gini_disp", "democracy", "interaction")

control_list_1 <- c("lag(tradePctGDP)", "lag(lpop,1)", "lag(lgdp,1)")


control_list_2<- c("lag(e_miinterc,1)","lag(e_miinteco,1)",
                        "lag(e_Civil_War,1)", "lag(e_migdpgro,1)", "lag(e_miurbani,1)")

control_list_3 <- c("lag(AgeDependency,1)","lag(PoliticalTerrorState,1)")


dv_list <- c("TradeTaxPctIncomeTax","TradeTaxPctTotalTax")



#################################################
#### Analysis
#################################################

#data = data[data$democracy == 0,]
# Table 1
fmla1 <- as.formula(paste(paste(dv_list[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| Year | 0 | COWCODE")))

fmla2 <- as.formula(paste(paste(dv_list[2],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| Year | 0 | COWCODE")))

fmla3 <- as.formula(paste(paste(dv_list[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| Year | 0 | COWCODE")))

fmla4 <- as.formula(paste(paste(dv_list[2],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| Year | 0 | COWCODE")))


fmla5 <- as.formula(paste(paste(dv_list[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| COWCODE + Year | 0 | COWCODE")))

fmla6 <- as.formula(paste(paste(dv_list[2],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| COWCODE + Year | 0 | COWCODE")))

fmla7 <- as.formula(paste(paste(dv_list[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| COWCODE + Year | 0 | COWCODE")))

fmla8 <- as.formula(paste(paste(dv_list[2],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| COWCODE + Year | 0 | COWCODE")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)
# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "fiscal capacity",
          #covariate.labels = c("Universities pc",
           #                    "GDP pc",
          #                     "Armed Conflict Internal",
          #                     "Armed Conflict International",
          #                     "Log(Population)",
           #                    "Avg Years of Schooling",
           #                    "Inequality",
           #                    "Inequality Squared"),
          dep.var.caption  = "Fiscal Capacity",
          dep.var.labels.include = FALSE,
          column.labels   = c("TradeTaxPctIncomeTax","TradeTaxPctTotalTax","TradeTaxPctIncomeTax","TradeTaxPctTotalTax","TradeTaxPctIncomeTax","TradeTaxPctTotalTax","TradeTaxPctIncomeTax","TradeTaxPctTotalTax"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "NO", "No", "No", "No","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes")),
          digits=2,
          out="~/Desktop/table1_main.tex"
)


# Table 1b
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)

# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse

# Table
stargazer(m1,m2,m3,m4,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table1b_main.tex"
)

# Table 1c
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_4,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared",
                               "Rev % of GDP"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table1c_main.tex"
)


# Table 1c
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),"+",
                          paste(control_list_5,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared",
                               "Oil Income pc"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table1d_main.tex"
)



# Table 2
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))


fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Public Universities pc",
                               "Private Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table2_main.tex"
)



# Table 2
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))


fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Public Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table2b_main.tex"
)

# Table 2
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))


fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Private Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table2c_main.tex"
)



# Table 3
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))


fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table3_main.tex"
)







# Table 4
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))


fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Public Universities",
                               "Private Universities",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table4_main.tex"
)




# Table 5
fmla1 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[4],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[19],"~"),
                          paste("lag(",treat_list[5],",1)","+"),
                          paste("lag(",treat_list[6],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "Universities",
                               "Public Universities pc",
                               "Private Universities pc",
                               "Public Universities",
                               "Private Universities",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Coups","Coups","Coups","Coups","Coups","Coups","Coups","Coups"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table5_coups.tex"
)


##############################################
##############################################
##### Heterogenous Effect by Time?
##############################################
##############################################

# Table 6
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+year+","lag(",treat_list[1],",1)*year","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "Year",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared",
                               "Year*Universities pc"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table6_main_time.tex"
)


# Table 6
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)+","lag(",treat_list[1],",1)*decades","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "1910s",
                               "1920s",
                               "1930s",
                               "1940s",
                               "1950s",
                               "1960s",
                               "1970s",
                               "1980s",
                               "1990s",
                               "2000s",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared",
                               "1910s*Universities pc",
                               "1920s*Universities pc",
                               "1930s*Universities pc",
                               "1940s*Universities pc",
                               "1950s*Universities pc",
                               "1960s*Universities pc",
                               "1970s*Universities pc",
                               "1980s*Universities pc",
                               "1990s*Universities pc",
                               "2000s*Universities pc"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table7_main_time2.tex"
)


# Table 8
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table8_main.tex"
)



################################
### State Capacity
################################

# Table 9
fmla1 <- as.formula(paste(paste(dv_list_state[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_state[2],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_state[3],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_state[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_state[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_state[2],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_state[3],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_state[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Public Sector Corr","Executive Corr","Public Sector Theft","Impartial","Public Sector Corr","Executive Corr","Public Sector Theft","Impartial"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table9_state.tex"
)


# Table 9
fmla1 <- as.formula(paste(paste(dv_list_state[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_state[2],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_state[3],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_state[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_state[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_state[2],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_state[3],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_state[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Public Universities pc",
                               "Private Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Public Sector Corr","Executive Corr","Public Sector Theft","Impartial","Public Sector Corr","Executive Corr","Public Sector Theft","Impartial"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table9b_state.tex"
)





# Table 10
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[1],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared",
                               "Impartial Admin"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table10_main_state.tex"
)



# Table 10
fmla1 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla2 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla3 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla4 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla5 <- as.formula(paste(paste(dv_list_demo[1],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla6 <- as.formula(paste(paste(dv_list_demo[4],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla7 <- as.formula(paste(paste(dv_list_demo[13],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))

fmla8 <- as.formula(paste(paste(dv_list_demo[16],"~"),
                          paste("lag(",treat_list[2],",1)","+"),
                          paste("lag(",treat_list[3],",1)","+"),
                          paste(control_list_1,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_2,collapse="+"),
                          paste("+",sep=""),
                          paste(control_list_3[4],collapse="+"),
                          paste("| CountryCode + year | 0 | CountryCode")))



# Estimation
m1 <- felm(fmla1, data=data)
m2 <- felm(fmla2, data=data)
m3 <- felm(fmla3, data=data)
m4 <- felm(fmla4, data=data)
m5 <- felm(fmla5, data=data)
m6 <- felm(fmla6, data=data)
m7 <- felm(fmla7, data=data)
m8 <- felm(fmla8, data=data)


# Standard errors
se1 <- m1$cse
se2 <- m2$cse
se3 <- m3$cse
se4 <- m4$cse
se5 <- m5$cse
se6 <- m6$cse
se7 <- m7$cse
se8 <- m8$cse

# Table
stargazer(m1,m2,m3,m4,m5,m6,m7,m8,type="latex",style="qje", 
          title            = "Democratization",
          covariate.labels = c("Public Universities pc",
                               "Private Universities pc",
                               "GDP pc",
                               "Armed Conflict Internal",
                               "Armed Conflict International",
                               "Log(Population)",
                               "Avg Years of Schooling",
                               "Inequality",
                               "Inequality Squared",
                               "Impartial Admin"),
          dep.var.caption  = "Democratization",
          dep.var.labels.include = FALSE,
          column.labels   = c("Boix","Polity2","Vdem","Cheibub","Boix","Polity2","Vdem","Cheibub"),
          se=list(se1,se2,se3,se4,se5,se6,se7,se8),
          add.lines = list(c("Country FE", "Yes", "Yes","Yes", "Yes","Yes", "Yes","Yes", "Yes"),
                           c("Year FE","Yes","Yes","Yes","Yes","Yes","Yes","Yes", "Yes")),
          digits=2,
          out="~/Dropbox/UniversitiesCapacity/Tables/table10b_main_state.tex"
)



## run public/private separately
# run biavriate model
# public/private for state capacity
# add tax/dgp and oil as controls

