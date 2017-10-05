rm(list=ls(all=TRUE))

library(cshapes)
library(countrycode)
library(WDI)
library(foreign)
library(sbgcop)
library(tidyverse)
library(haven)
source("~/Dropbox/Duke/Rcode/lagger.R")

setwd("~/Dropbox/FiscalCapacityCrossNational/Data/")

### extract world bank data and rename variables
wb.tax = c("NY.GDP.MKTP.CN","NE.TRD.GNFS.ZS","NV.AGR.TOTL.ZS","NY.GDP.MKTP.KD","NY.GDP.PCAP.KD","NV.IND.TOTL.ZS","NV.IND.MANF.ZS","MS.MIL.XPND.GD.ZS","SP.POP.DPND","AG.LND.AGRI.ZS","SL.AGR.EMPL.ZS","SL.IND.EMPL.ZS","NE.GDI.FTOT.ZS","DT.ODA.ODAT.GN.ZS","NY.GDP.PETR.RT.ZS","SP.POP.TOTL")
wb.code = c("GDPLCU","tradePctGDP","AgroValueAddedGDP","GDPConst2005US","GDPpCConst2005US","IndustValueAddedGDP","ManufacValueAddedGDP","MilExpGDP","AgeDependency","AgrpPctLand","EmploymentPctAgro","EmploymentPctIndust","GrossFixedCapFormationPctGDP","NetODAPctGNI","NetOilRentPctGDP","PopulationTotal")
wb = WDI(country="all",indicator=wb.tax,start=1980,end=2014,extra=TRUE)     
summary(wb)
names(wb)[3] = "Year"
names(wb)[4:19] = wb.code[wb.tax==names(wb)[4:19]]


### drop variables that are not needed
notwant = which(names(wb)%in%c("GDPLCU" ,"region","capital","income","lending"))
wb = wb[,-notwant]
### drop regions 
notwant = which(wb$country %in% c("Arab World", "Europe & Central Asia (excluding high income)", "Latin America & the Caribbean (IDA & IBRD countries)", "Middle East & North Africa (IDA & IBRD countries)", 
                                  "East Asia & Pacific (IDA & IBRD countries)", "South Asia (IDA & IBRD)", "Sub-Saharan Africa (IDA & IBRD countries)", "Europe & Central Asia (IDA & IBRD countries)", 
                                  "East Asia & Pacific (excluding high income)", "World", "East Asia & Pacific (developing only)", "Europe & Central Asia (developing only)", "South Asia", "Central Europe and the Baltics",
                                  "European Union", "Fragile and conflict affected situations", "OECD members",  "Small states", "Pacific island small states", "Caribbean small states", "Other small states",
                                  "Euro area", "High income", "Heavily indebted poor countries (HIPC)", "Latin America & Caribbean (developing only)","Least developed countries: UN classification",
                                  "Low income", "Lower middle income", "Low & middle income", "Middle income", "Middle East & North Africa (developing only)", "High income: nonOECD", "High income: OECD", "Upper middle income", "North America", "Not classified",  "East Asia & Pacific (all income levels)", "Europe & Central Asia (all income levels)",  "Sub-Saharan Africa (developing only)", "Sub-Saharan Africa (all income levels)"
                                  , "Latin America & Caribbean (all income levels)", "Middle East & North Africa (all income levels)", "IBRD only", "IDA total", "IDA blend", "IDA & IBRD total","Latin America & Caribbean (excluding high income)",
                                  "IDA only", "Middle East & North Africa (excluding high income)", "Sub-Saharan Africa" , "Latin America & Caribbean", "Sub-Saharan Africa (excluding high income)", "Middle East & North Africa", 
                                  "Pre-demographic dividend", "Early-demographic dividend", "Late-demographic dividend", "East Asia & Pacific", "Europe & Central Asia", "Post-demographic dividend"
                                  ))

wb = wb[-notwant,]

### now add cowcode to wb data
wb$COWCODE = countrycode(wb$iso2c, "iso2c", "cown", warn=T)
### which ones are problematic
 ### Some values were not matched unambiguously: AS, AW, BM, CW, FO, GI, GL, GU, HK, IM, JG, KY, MF, MO, MP, NC, PF, PR, PS, RS, SX, TC, VI, XK
wb$COWCODE[wb$iso2c == "AS"] = NA ## american samoa
wb$COWCODE[wb$iso2c == "AW"] = NA## aruba
wb$COWCODE[wb$iso2c == "BM"] = NA ## Bermuda
wb$COWCODE[wb$iso2c == "CW"] = NA## curacao
wb$COWCODE[wb$iso2c == "FO"] = NA## Faroe Islands
wb$COWCODE[wb$iso2c == "GI"] = NA## GIBRALTAR
wb$COWCODE[wb$iso2c == "GL"] = NA## Greenland
wb$COWCODE[wb$iso2c == "GU"] = NA## GUAM
wb$COWCODE[wb$iso2c == "HK"] = NA## Hong kong
wb$COWCODE[wb$iso2c == "IM"] = NA## Isle of man
wb$COWCODE[wb$iso2c == "JG"] = NA## Channel Islands
wb$COWCODE[wb$iso2c == "KY"] = NA## Cayman Islands
wb$COWCODE[wb$iso2c == "MF"] = NA## St. Martin (French part)
wb$COWCODE[wb$iso2c == "MO"] = NA## Macao SAR, China
wb$COWCODE[wb$iso2c == "MP"] = NA## Northern Mariana Islands
wb$COWCODE[wb$iso2c == "NC"] = NA## New Caledonia 
wb$COWCODE[wb$iso2c == "PF"] = NA## French Polynesia
wb$COWCODE[wb$iso2c == "PR"] = NA## Puerto Rico 
wb$COWCODE[wb$iso2c == "PS"] = NA## West Bank and Gaza
wb$COWCODE[wb$iso2c == "RS"] = 345 ## Serbia
wb$COWCODE[wb$iso2c == "SX"] = NA## Sint Maarten (Dutch part) 
wb$COWCODE[wb$iso2c == "TC"] = NA##  Turks and Caicos Islands
wb$COWCODE[wb$iso2c == "VI"] = NA## Virgin Islands (U.S.)
wb$COWCODE[wb$iso2c == "XK" & wb$Year > 2007 ] = 347## Kosovo
### drop observations without COWCODE
wb = wb[is.na(wb$COWCODE) == FALSE, ]



##### adding tax data from http://www.ictd.ac/en/ictd-government-revenue-dataset
ictd = read_dta("~/Dropbox/DATA/ICTD/StataFiles/merged.dta")
summary(ictd)
ictd = ictd[, c("country","iso", "year", "totrev",  "tottax", "totnontax", "direct", "income", "indiv", "corp", "indirect", "gst", "trade", "caution2resourcerevenuestax", "caution1accuracyqualityorco", "caution3unexcludedresourcere", "caution4inconsistencieswiths")]




### clearning ictd

#names(ictd)[names(ictd) == "iso"] = "wbcode"
names(ictd)[names(ictd) == "year"] = "Year"

ictd$COWCODE = countrycode(ictd$iso, "wb", "cown", warn=T)
###  ABW, AIA, COD, HKG, MAC, ROU, SRB, WBG
ictd$COWCODE[ictd$iso == "ABW"] = NA #Aruba
ictd$COWCODE[ictd$iso == "AIA"] = NA #Anguila
ictd$COWCODE[ictd$iso == "COD"] = NA #?

ictd$COWCODE[ictd$iso == "HKG"] = NA #Hong Kong
ictd$COWCODE[ictd$iso == "MAC"] = NA #Macao
ictd$COWCODE[ictd$iso == "ROU"] = NA #?
ictd$COWCODE[ictd$iso == "SRB"] = 345 #Serbia
ictd$COWCODE[ictd$iso == "WBG"] = NA #West bank, gaza strip
### drop observations without COWCODE
ictd = ictd[is.na(ictd$COWCODE) == FALSE, ]


names(ictd)
### additional tax variables to create
ictd$TradeTaxPctTotalTax = ictd$trade/ ictd$tottax
ictd$IncomeTaxPctTotalTax = ictd$income/ ictd$tottax
ictd$IncomeTaxPctTotalTax = ifelse(is.infinite(ictd$IncomeTaxPctTotalTax) == T, NA, ictd$IncomeTaxPctTotalTax)
ictd$TradeTaxPctIncomeTax = ictd$trade/ ictd$income
ictd$TradeTaxPctIncomeTax = ifelse(is.infinite(ictd$TradeTaxPctIncomeTax) == T, NA, ictd$TradeTaxPctIncomeTax)

### join ictd and wb

data = join(wb,ictd, by = c("COWCODE","Year"), type = "left")
data[duplicated(data[,c("COWCODE","Year")]),c("COWCODE","Year")]
data <- data[order(data$COWCODE, data$Year),]

View(data[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade")])

#### regime type data set from geddes, wright, frantz, 
gwf = read_dta("~/Dropbox/DataForCrossNational/GWFtscs.dta")
summary(gwf)
names(gwf)[c(1,2)] = c("COWCODE", "Year")
data = join(data,gwf, by = c("COWCODE","Year"), type = "full")
data[duplicated(data[,c("COWCODE","Year")]),c("COWCODE","Year")]



View(data[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade")])


### boix, miller, rosato data
dem.boix = read.dta("~/Dropbox/DATA/Boix_etal/democracy-v2.0.dta")
dem.boix = dem.boix[,c(1,2,4,5,6,9)]
names( dem.boix)[2] = "COWCODE"
names( dem.boix)[3] = "Year"
dem.boix$COWCODE[dem.boix$country == "GERMANY, WEST" &  dem.boix$Year<1990] = 255
dem.boix$COWCODE[dem.boix$country == "USSR" &  dem.boix$Year<1993] = 365

data = join(data,dem.boix,by=c("COWCODE","Year"), type = "left")
data[duplicated(data[,c("COWCODE","Year")]),c("COWCODE","Year")]


data <- data[order(data$COWCODE, data$Year),]
View(data[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade", "democracy")])



### ICRG bureaucratic quality data 

bq = read.csv("~/Dropbox/DATA/ICRG/bureaucratic_quality.csv")
bq$PubDate = as.Date(bq$PubDate, format = "%d/%m/%Y" )
bq$Year = format(bq$PubDate, "%Y")
bq = bq[,-1]
### rearrange dataset
bq = gather(bq,"Year")
glimpse(bq)
head(bq)
names(bq)  = c("Year","Country","BureaucraticQuality")

bq = ddply(bq, .variables = c("Country","Year"), summarise, BQuality = mean(BureaucraticQuality, na.rm = TRUE))
bq$BQuality[is.nan(bq$BQuality)==TRUE] = NA

bq$COWCODE = countrycode(bq$Country, "country.name", "cown", warn=T)
### warning: Some values were not matched: Hong.Kong, New.Caledonia, Serbia, Serbia...Montenegro, some strings were matched more than once: West.Germany,260,255
bq$COWCODE[bq$Country == "Serbia"] = 345 #Serbia
bq$COWCODE[bq$Country == "Germany" & bq$Year < 1991] = NA
bq$COWCODE[bq$Country == "West.Germany" & bq$Year > 1990] = NA
bq = bq[is.na(bq$COWCODE) == FALSE, ]
bq$COWCODE[bq$Country == "Russia" & bq$Year<1992] = NA
bq$COWCODE[bq$Country == "USSR" & bq$Year>1991] = NA
bq$COWCODE[bq$Country == "Korea..DPR"] = 731
data = join(data, bq, by=c("COWCODE","Year"), type = "left")
data[duplicated(data[, c("COWCODE", "Year")]), c("Country", "COWCODE", "Year")]

data <- data[order(data$COWCODE, data$Year),]
View(data[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade", "democracy", "BQuality")])



##### QoG data
qog = read.dta("~/Dropbox/DATA/QoG/qog_std_ts_jan15.dta")

names(qog)
WANT_qog = which(names(qog) %in%c("ccodecow","year", "cam_contest", "chga_demo", "chga_hinst", "ciri_physint", "ciri_polpris", "fh_ipolity2","gd_ptsa", "gd_ptss", "p_polity2", "uds_mean", "uds_median", "wdi_armedfper", "wdi_expmilgdp"))
qog_merge = qog[, WANT_qog]
names(qog_merge)
names(qog_merge) = c("Year", "COWCODE", "Contestation", "CheibubDemocracy", "CheibubRegimeInstitution", "CiriPhysicalIntegrity", "CiriPoliticalImprisonment", "FHPolity2Imputed", "PoliticalTerrorAmnesty", "PoliticalTerrorState", "Polity2","UnifiedDemocracyMean", "UnifiedDemocracyMedian", "ArmedPersonalPctLabor", "ExpenditureMilitary")

data = join(data,qog_merge, by=c("Year","COWCODE"), type = "left")
rm(qog_merge)
rm(qog)
data[duplicated(data[, c("COWCODE", "Year")]), c("Country", "COWCODE", "Year")]

View(data[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade", "democracy", "BQuality", "CheibubDemocracy")])

#### now create war variable based on ucdp data
vdem = read_dta("~/Dropbox/DATA/VDEM/Country_Year_V-Dem_other_STATA_v7.1/V-Dem-DS-CY+Others-v7.1.dta")
#vdem[] <- lapply(vdem, unclass)
vdem = vdem[, c("year", "COWcode", "e_Regime", "e_polity_s", "e_peginiwi", "e_region_dem_diffuse", "e_area", "e_migdppc", "e_migdpgro", "e_cow_exports", "e_cow_imports", "e_population", "e_miurbani", "e_miinteco", "e_miinterc", "e_Civil_War")]
names(vdem)[names(vdem) == "COWcode"] = "COWCODE"
names(vdem)[names(vdem) == "year"] = "Year"

data = join(data, vdem, by =c("COWCODE", "Year"), type = "left")
summary(data)



#View(data[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade", "democracy", "BQuality", "CheibubDemocracy", "e_Regime")])
data[duplicated(data[, c("COWCODE", "Year")]), c("Country", "COWCODE", "Year")]

write_dta(data, "FiscalCapacityData.dta")
save(data, file = "FiscalCapacityData.rda")



# 
# 
# inequ <- read_dta("~/Dropbox/DATA/wiid/WIID3.4_19JAN2017.dta")
# head(inequ)
# summary(inequ)
# inequ <- inequ[, c("Countrycode3", "Year", "Gini" )]
# inequ$COWCODE = countrycode(inequ$Countrycode3, "wb", "cown", warn=T)
# ###DHY, HKG, KSV, PRI, REU, SCG, SRB, SUN, WBG
# inequ$COWCODE[inequ$Countrycode3 == "SRB"] = 345 #Serbia
# data = join(data, inequ, by =c("COWCODE", "Year"), type = "left")
# 
# write_dta(data, "FiscalCapacityData.dta")
# save(data, file = "FiscalCapacityData.rda")

### SWIID multiple datasets

load("~/Dropbox/DATA/swiid6_0/swiid6_0.rda")
# 
# ### write function to do cleaning and merging of swiid, can then be used for all 100 datasets
# prep_SWIID = function(x){  
#   x$COWCODE=countrycode(x$country,"country.name","cown",warn=F)
#   ### problematic Anguilla, Cayman Islands, Hong Kong, Puerto Rico, Serbia, Serbia and Montenegro, Turks and Caicos
#   x$COWCODE[x$country == "USSR"] = NA 
#   x$COWCODE[x$country == "Czechoslovakia" & x$year>1986] = NA
#   x$COWCODE[x$country == "Montenegro"] = NA
#   x$COWCODE[x$country == "Serbia"] = 345
#   
#   ineq = x[,c("year", "gini_disp", "gini_mkt", "COWCODE")]
#   names(ineq) = c("Year","GiniNet", "GiniMarket", "COWCODE")
#   
#   x = join(data,ineq, by=c("Year","COWCODE"), type = "left")
#   countries = unique(x$COWCODE)
#   x$all.missing = 0
#   for(i in countries){
#     x$all.missing[x$COWCODE==i] = ifelse(all(is.na(x$TradeTaxPctTotalTax[x$COWCODE==i]))==TRUE,1,x$all.missing[x$COWCODE==i])    
#     x$all.missing[x$COWCODE==i] = ifelse(all(is.na(x$GiniMarket[x$COWCODE==i]))==TRUE,1,x$all.missing[x$COWCODE==i])    
#   }
#   x = subset(x, all.missing==0)
#   ### delete all the country names
#   notwant = which(names(x) %in%c("Country", "wbcode", "country", "iso2c", "iso3c", "country1", "Country.1", "cty_name", "iso3numeric"))
#   x = x[,-notwant]
#   return(x)
# }
# 
# 
# multi_data = lapply(swiid, FUN=prep_SWIID)
# for(i in 1:100){
#   print(dim(multi_data[[i]]))
# }
# save(multi_data, file = "multi_data.rda")
# 
# #### create a dataset with the gini mean


names(swiid_summary)[names(swiid_summary) == "year"] =  "Year"
swiid_summary$COWCODE=countrycode(swiid_summary$country,"country.name","cown",warn=T)

###Anguilla, Hong Kong, Micronesia, Palestine, Puerto Rico, Serbia, Turks and Caicos Islands 
swiid_summary$COWCODE[swiid_summary$country == "Hong Kong"] = NA 
swiid_summary$COWCODE[swiid_summary$country == "Micronesia"] = NA
swiid_summary$COWCODE[swiid_summary$country == "Palestine"] = NA
swiid_summary$COWCODE[swiid_summary$country == "Puerto Rico"] = NA
swiid_summary$COWCODE[swiid_summary$country == "Turks and Caicos Islands"] = NA
swiid_summary$COWCODE[swiid_summary$country == "Serbia"] = 345

### duplicate values for Russia for 1988, 1989, 1990 because of su and russia
swiid_summary$COWCODE[swiid_summary$country == "Russia"  & swiid_summary$Year < 1991] = NA

data_mean_inequ = join(data,swiid_summary, by=c("Year","COWCODE"), type = "left")

data_mean_inequ <- data_mean_inequ[order(data_mean_inequ$COWCODE, data_mean_inequ$Year),]
View(data_mean_inequ[,c("COWCODE", "Year", "country", "TradeTaxPctIncomeTax", "trade", "democracy", "BQuality", "gini_disp", "Gini")])





# countries = unique(data$COWCODE)
# data_mean_inequ$all.missing = 0
# for(i in countries){
#   data_mean_inequ$all.missing[data_mean_inequ$COWCODE==i] = ifelse(all(is.na(data_mean_inequ$TradeTaxPctTotalTax[data_mean_inequ$COWCODE==i]))==TRUE, 1, data_mean_inequ$all.missing[data_mean_inequ$COWCODE==i])    
#  data_mean_inequ$Country[data_mean_inequ$COWCODE==i] = ifelse(all(is.na(data_mean_inequ$gini_disp[data_mean_inequ$COWCODE==i]))==TRUE, 1, data_mean_inequ$all.missing[data_mean_inequ$COWCODE==i])    
# }
# data_mean_inequ = subset(data_mean_inequ, all.missing==0)
# ### delete all the country names
# notwant = which(names(data_mean_inequ) %in%c("Country", "wbcode", "country", "iso2c", "iso3c", "country1", "Country.1", "cty_name", "iso3numeric"))
# data_mean_inequ = data_mean_inequ[,-notwant]
save(data_mean_inequ, file = "data_mean_inequ.rda")
write_dta(data_mean_inequ, "data_inequ.dta")

