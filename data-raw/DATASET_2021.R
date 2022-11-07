## code to prepare `DATASET` dataset goes here

#EXPLORE WATERBASE WATERQUALITY DATABASE. COMPUTE MCR AND MAF. PRODUCE PLOTS AND SUMMARY TABLES
#ismaelm.rodeapalomares@bayer.com
#03/26/2021

#cleaning the environemt
rm(list=ls())

#Define and set paths
## this need to be changed by the user if reproducing the base data in the package
Dir_Data = "~/s3raw/WB-Zhenglei/"


#Required packages

library(tidyverse)


#####################################################################################################
# LOAD DATASETS
#####################################################################################################

# Load CAS for heavy metals and priority pollutants
CAS.Metals = read.csv(paste0(Dir_Data,"CAS.Metals.csv"))
CAS.Metals$CAS = as.character(CAS.Metals$CAS)
CAS.Priority = read.csv(paste0(Dir_Data,"EU_WFD_45PriorityPollutants.csv"))
CAS.Priority$CAS = gsub("CAS_|-|_", "", CAS.Priority$CAS)
CAS.Priority$CAS = as.character(CAS.Priority$CAS)

# Load Endpoint data
#L.Posthuma SSDs
data(Data.SSD)
Data.SSD <- read.csv(paste0(Dir_Data,"Copy of etc4373-sup-0002-supmat.csv"))

#Load data. Basic exploration of stations information

Stations = read.csv(paste0(Dir_Data,"Waterbase_v2020_1_S_WISE6_SpatialObject_DerivedData.csv")) # v2021 ==> too many missing information use the previous one instead.
Stations = read.csv(paste0(Dir_Data,"Waterbase_v2019_1_S_WISE6_SpatialObject_DerivedData_repaired.csv")) # v2021


names(Stations)
# names(Stations)[1] = "monitoringSiteIdentifier"

# "" to NA


#How many Stations?
Stations$monitoringSiteIdentifier = as.factor(Stations$monitoringSiteIdentifier)
length(levels(Stations$monitoringSiteIdentifier))
#60776
#How many have complete lat/long info?
sum(complete.cases(Stations[names(Stations)%in%c("lon", "lat")]))
#48563


###################################################################################################
# Chemical monitoring data curation
###################################################################################################

library(data.table)
Data.AggBySiteID = fread(paste0(Dir_Data,"Waterbase_v2020_1_T_WISE6_AggregatedData.csv"))
## Data.AggBySiteID = read.csv(paste0(Dir_Data,"Waterbase_v2020_1_T_WISE6_AggregatedData.csv")) #2021 #3962763 rows

dim(Data.AggBySiteID) ## [1] 3510775      31
# Data.AggBySiteID = Data.AggByWB.1
# rm(Data.AggByWB.1)
#2021 #3962763 rows


# Data.AggBySiteID = Data.AggByWB.1
# rm(Data.AggByWB.1)

names(Data.AggBySiteID)[1] = "monitoringSiteIdentifier"

# Data.AggBySiteID[Data.AggBySiteID==""] = NA

Data.AggBySiteID$monitoringSiteIdentifier = as.factor(Data.AggBySiteID$monitoringSiteIdentifier)

# How many of the sites are in the "Sites" dataset?
sum(levels(Data.AggBySiteID$monitoringSiteIdentifier)%in%levels(Stations$monitoringSiteIdentifier))
# 15990 Sites
# How many have lat long info?
sum(levels(Data.AggBySiteID$monitoringSiteIdentifier)%in%Stations$monitoringSiteIdentifier[complete.cases(Stations[names(Stations)%in%c("lon", "lat")])])
#13410 Sites have lat long info

head(Data.AggBySiteID)
names(Data.AggBySiteID)
# [1] "monitoringSiteIdentifier"             "monitoringSiteIdentifierScheme"       "parameterWaterBodyCategory"
# [4] "observedPropertyDeterminandCode"      "observedPropertyDeterminandLabel"     "procedureAnalysedMatrix"
# [7] "resultUom"                            "phenomenonTimeReferenceYear"          "parameterSamplingPeriod"
# [10] "procedureLOQValue"                    "resultNumberOfSamples"                "resultQualityNumberOfSamplesBelowLOQ"
# [13] "resultQualityMinimumBelowLOQ"         "resultMinimumValue"                   "resultQualityMeanBelowLOQ"
# [16] "resultMeanValue"                      "resultQualityMaximumBelowLOQ"         "resultMaximumValue"
# [19] "resultQualityMedianBelowLOQ"          "resultMedianValue"                    "resultStandardDeviationValue"
# [22] "procedureAnalyticalMethod"            "parameterSampleDepth"                 "resultObservationStatus"
# [25] "remarks"                              "metadata_versionId"                   "metadata_beginLifeSpanVersion"
# [28] "metadata_statusCode"                  "metadata_observationStatus"           "metadata_statements"
# [31] "UID"
summary(Data.AggBySiteID$phenomenonTimeReferenceYear)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1931    2007    2011    2009    2015    2019
barplot(table(Data.AggBySiteID$phenomenonTimeReferenceYear))

# How many determinads included in the L.Posthuma SSDs?
# Prepare CAS numbers
Data.AggBySiteID$CAS = gsub("CAS_|-|_", "", Data.AggBySiteID$observedPropertyDeterminandCode)

#How many CAS numbers in SSDs?
length(unique(Data.AggBySiteID$CAS)) #552 determinads


# How many entries from the exisiting data we will be able to use?
CAS.SSD = unique(Data.AggBySiteID$CAS)[unique(Data.AggBySiteID$CAS)%in%as.character(unique(Data.SSD$CAS.))]
# 373 CAS
length(Data.AggBySiteID$CAS[Data.AggBySiteID$CAS%in%CAS.SSD])
#2254186 entires out of 3510775
length(Data.AggBySiteID$CAS[Data.AggBySiteID$CAS%in%CAS.SSD])/nrow(Data.AggBySiteID)
# 64 % of the data

#Wich are the "CAS" not matched?
CAS.SSD.no = unique(Data.AggBySiteID$CAS)[!unique(Data.AggBySiteID$CAS)%in%as.character(unique(Data.SSD$CAS.))]
table(grepl("EE",CAS.SSD.no))
# FALSE  TRUE
# 65    58
# 65 determinads have a wierd code instead of the CAS number

# Further reduce by the type of sample analyze
table(Data.AggBySiteID$procedureAnalysedMatrix)
# W   W-DIS   W-SPM
# 3469602  492896     265
# http://dd.eionet.europa.eu/vocabulary/wise/Matrix/view?page=2#vocabularyConceptResults
# W: Water - Total
# W-DIS: Water - Dissolved (filtered)
# W-SPM:  	Water - Suspended particulate matter
# Keep only "W" to avoid replications
Data.AggBySiteID = Data.AggBySiteID[Data.AggBySiteID$procedureAnalysedMatrix=="W",] #v2020

# Keep only correct entries without QA flags
# resultObservationStatus,"Status of the observed value in terms of its availability, relevancy, correctness or specifics of its source category.","A: Record is confirmed as correct;
## http://dd.eionet.europa.eu/vocabulary/wise/Matrix/view?page=2#vocabularyConceptResults
# A: Record is confirmed as correct;
# L: Missing observed value, the data were not collected;
# M: Missing observed value, the data can not exist;
# N: Missing observed value, observed value is not relevant or not significant;
# O: Missing observed value, no further information is available or record reported in the past should be deleted;
# W: Missing observed value, data are included in another source category;
# X: Reported value includes data from another source category (categories);
# Y: The source category does not exactly match the standard definition
table(Data.AggBySiteID$resultObservationStatus)
# records validated as correct very few, but no records verified as incorrect
# A
# 383780
Data.AggBySiteID =Data.AggBySiteID[Data.AggBySiteID$resultObservationStatus!="O",]

# metadata_observationStatus,Status of the record regarding its reliability.
# A: Normal record;
# U: Record with lower reliability;
# V: Unvalidated record"
# metadata_statements,Aditional information or statements regarding the reliability of the feature or record.
table(Data.AggBySiteID$metadata_observationStatus)
# A       U
# 3343206  126396
#126396  unreliable records
Data.AggBySiteID =Data.AggBySiteID[Data.AggBySiteID$metadata_observationStatus!="U",]

Data.AggBySiteID = droplevels(Data.AggBySiteID)

#metadata_statements
Data.AggBySiteID$metadata_statements = as.factor(Data.AggBySiteID$metadata_statements)
metadata_statements = levels(Data.AggBySiteID$metadata_statements)
# [1] ""
# [2] "NOTE_LEGACY: Measurement has been confirmed by country to be taken from a highly polluted area "
# [3] "NOTE_LEGACY: Outlier has been confirmed by country as correct value"
# [4] "NOTE_LEGACY: The parameterWaterBodyCategory changed from LW to RW based on the the monitoring site data"
# [5] "NOTE_LEGACY: The parameterWaterBodyCategory changed from RW to LW based on the the monitoring site data"
# [6] "NOTE_WBCAT_CHANGE: The parameterWaterBodyCategory has been changed to match the one officialy reported for the given spatial identifier (old: LW, new: RW)"
# [7] "NOTE_WBCAT_CHANGE: The parameterWaterBodyCategory has been changed to match the one officialy reported for the given spatial identifier (old: LW, new: TW)"

levels(Data.AggBySiteID$metadata_statements)[2]
table(Data.AggBySiteID$metadata_statements!="NOTE_LEGACY: Measurement has been confirmed by country to be taken from a highly polluted area ")
# FALSE    TRUE
# 85 3343121
#Remove 85 because MAF is not to cover hot spot contamination
Data.AggBySiteID =Data.AggBySiteID[Data.AggBySiteID$metadata_statements!="NOTE_LEGACY: Measurement has been confirmed by country to be taken from a highly polluted area",]

# Waterbody category
# http://dd.eionet.europa.eu/vocabulary/wise/WFDWaterBodyCategory/
# Coastal water body		CW
# Groundwater body		GW
# Lake water body		LW
# River water body		RW
# Transitional water body		TW
# Territorial waters		TeW
# Marine waters		MW

#only LW, RW
table(Data.AggBySiteID$parameterWaterBodyCategory)
# CW      GW      LW      RW      TW
# 322  343589  387473 2602523    9214
Data.AggBySiteID =Data.AggBySiteID[Data.AggBySiteID$parameterWaterBodyCategory%in%c("LW", "RW"),]


length(unique(Data.AggBySiteID$monitoringSiteIdentifier))
#[1] 17984
length(unique(Data.AggBySiteID$CAS))
#[1] 496



st_tmp <- Stations[Stations$monitoringSiteIdentifier %in% unique(Data.AggBySiteID$monitoringSiteIdentifier),]

length(unique(st_tmp$countryCode))
length(unique(st_tmp$waterBodyIdentifier))




################################ End of Data Selection ##################################










# CAREFULL, THIS SUBSET IS INTERIM, ONLY TO MAKE IT EASY TO MANIPULATE -------



# Subset for the AIs that were matched (Data.AggSiteID.1)
Data.AggBySiteID.1 = Data.AggBySiteID[Data.AggBySiteID$CAS%in%CAS.SSD,];
Data.AggBySiteID.1 = droplevels(Data.AggBySiteID.1) #1816359
Data.SSD.1 = Data.SSD[Data.SSD$CAS.%in%CAS.SSD,]; Data.SSD.1 = droplevels(Data.SSD.1) #373

# Units
table(Data.AggBySiteID.1$resultUom)
# mg/L mg{NH4}/L   mg{P}/L   mg{S}/L      ug/L
# 97814    102903    186232        10   1429400

#Remove all but mg/L and ug/L and transform all to mg/L
Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$resultUom%in%c("mg/L", "ug/L"),] #1469931
#transform to ug/L
Data.AggBySiteID.1$resultMinimumValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMinimumValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
Data.AggBySiteID.1$resultMeanValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMeanValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
Data.AggBySiteID.1$resultMedianValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMedianValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
Data.AggBySiteID.1$resultMaximumValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMaximumValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
# NOTE: Unit label is changed later after removing unreliable LOQs

# Subset for complete data for Average and Maximun values
Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultMeanValue),]
Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultMaximumValue),]

###################################################################################
# SSD data quality and curation (origing, and number of taxa per SSD)----------
###################################################################################

table(Data.SSD.1$OrigenQuality.Acute.EC50)
#Acute
# AllData
# 368
# ReadAcross
# 2
# Yellow highlight-Potential to extrapolate acute SSD median concentration from poorly represented chronic NOEC data (AcuteMu = ChronicMu+1, AcuteSigma = ChronicSigma or 0.7 as a default)
# 3
#Chronic
table(Data.SSD.1$OrigenQuality.Chronic.NOEC)
# AllData
# 351
# Green highlight-Potential to extrapolate Chronic NOEC SSD median concentration from poorly represented Acute EC50 data (ChronicMu = AcuteMu-1, ChronicSigma = AcuteSigma or 0.7 as a default)
# 22

#Number of taxa per SSD
summary(Data.SSD.1$X.Species.Acute.EC50)
#Acute
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    1.00    8.00   16.00   23.76   26.00  210.00       3 hist(log10(Data.SSD.1$X.Species.Acute.EC50))
hist(Data.SSD.1$X.TaxClass.Acute.EC50)
# How many more than 3 organism
table(Data.SSD.1$X.Species.Acute.EC50>=3)
# FALSE  TRUE
# 22   348
#CHR
summary(Data.SSD.1$X.Species.Chronic.NOEC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    1.00   10.00   15.00   19.44   24.00  182.00      22
hist(log10(Data.SSD.1$X.Species.Chronic.NOEC))
table(Data.SSD.1$X.Species.Chronic.NOEC>=3)
# FALSE  TRUE
# 5   346

# Remove chem that have less than 3 sp per SSD
Data.SSD.1 = Data.SSD.1[Data.SSD.1$X.Species.Acute.EC50>=3,]
Data.SSD.1 = Data.SSD.1[Data.SSD.1$X.Species.Chronic.NOEC>=3,]
Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X.Species.Chronic.NOEC),]
Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X.Species.Acute.EC50),] #346 chemicals


# Data quality
table(Data.SSD.1$OrigenQuality.Chronic.NOEC)
table(Data.SSD.1$OrigenQuality.Acute.EC50)

# Highest data quality in all cases
summary(Data.SSD.1$X.Species.Acute.EC50)
summary(Data.SSD.1$X.Species.Chronic.NOEC)
summary(Data.SSD.1$X.TaxClass.Acute.EC50)

# Taxa clasess
# Different taxa classes, need to evaluate what exaclty is this.
table(Data.SSD.1$X.TaxClass.Acute.EC50>=3)
# FALSE  TRUE
# 39   307
table(Data.SSD.1$X.TaxClass.Chronic.NOEC>=3)
# FALSE  TRUE
# 61   285
# Final interim set of chemicals for analysis
CAS.SSD.1 = CAS.SSD[CAS.SSD%in%Data.SSD.1$CAS.]
Data.SSD.1 = Data.SSD.1[Data.SSD.1$CAS.%in%CAS.SSD.1,]; Data.SSD.1 = droplevels(Data.SSD.1)

# SSD data with Slope data
Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X10LogSSDSlope.ug.L..SigmaChronic.NOEC),]
Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X10LogSSDSlope.ug.L..SigmaAcute.EC50),]

#########################################################################
# Subset Water-base rivers based on available SSDs

Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$CAS%in%CAS.SSD.1,];
Data.AggBySiteID.1 = droplevels(Data.AggBySiteID.1) #1679005 entries

summary(Data.SSD.1$X.Species.Chronic.NOEC)
summary(Data.SSD.1$X.Species.Acute.EC50)

##########################################################################
# Subset the dataset based on the number of measurements per CAS per year

#Number of samples
summary(Data.AggBySiteID.1$resultNumberOfSamples)
#Stations reporting number of samples
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   1.000   4.000  10.000   8.705  12.000 366.000      47

# Remove stations do not reporting number of samples
Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultNumberOfSamples),]

# Stations reporting more than 5 samples per year to report annual aggregates
Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$resultNumberOfSamples>=5,] #962409 v2021

###########################################################################
# LOQs: Remove unreliable entries based on missing/unreliable LOQs reporting

# Remove data for which LOQ is not reported
Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$procedureLOQValue),]

# Remove data with LOQ <= 0
Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$procedureLOQValue> 0,]

# N Samples below LOQ
summary(Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ)
# Remove data that do not report "resultQualityNumberOfSamplesBelowLOQ"
Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ),] #640197

#N samples below LOQ
summary(Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    6.00   11.00    9.24   12.00  366.00  283332

# Number of samples above LOQ
table(Data.AggBySiteID.1$resultNumberOfSamples > Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ)
# FALSE   TRUE
# 471043 169154
100*(471043/(471043+169154)) #73.5% of samples below LOQ; 30% actual measurements
# This is very important to see what to do wit BelowLOQ values!

#if N samples > 1, and N.detects > 2, then max > mean. Therefore, max > mean can be used to separate the samples with detects from those with non-detects
table(Data.AggBySiteID.1$resultMaximumValue>Data.AggBySiteID.1$resultMeanValue)
# FALSE   TRUE
# 461365 178832
# 100*(429987/(429987+160803))
# 72.78 below LOQ, 37% avobe LOQ
# The number should be equal to the value in the previous logic

# Convert reported LOQs to ug/L
Data.AggBySiteID.1$procedureLOQValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$procedureLOQValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
# Change unit label
Data.AggBySiteID.1$resultUom[Data.AggBySiteID.1$resultUom=="mg/L"] = "ug/L"; Data.AggBySiteID.1 = droplevels(Data.AggBySiteID.1)

#2. Set samples with Max <= LOQ to below LOQ
Data.AggBySiteID.1$AboveLOQ = Data.AggBySiteID.1$resultNumberOfSamples > Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ
Data.AggBySiteID.1$AboveLOQ[Data.AggBySiteID.1$resultMaximumValue<=Data.AggBySiteID.1$procedureLOQValue] = FALSE
summary(Data.AggBySiteID.1$AboveLOQ)

#Remove useless variables from datase
usefull.var = c("monitoringSiteIdentifier","phenomenonTimeReferenceYear", "Site.Year", "resultMeanValue", "resultMaximumValue","CAS","N.Measured","N.Detected","procedureLOQValue", "AboveLOQ", "resultNumberOfSamples", "resultQualityNumberOfSamplesBelowLOQ")
usefull.var1 <- usefull.var[usefull.var %in% names(Data.AggBySiteID.1)]
# Data.AggBySiteID.1 = Data.AggBySiteID.1[,names(Data.AggBySiteID.1)%in%usefull.var]
Data.AggBySiteID.1 = Data.AggBySiteID.1[,..usefull.var1]

# How many entries for the same CAS per Site_Year ##################################################
EntriesPerCAS = aggregate(resultMeanValue ~ CAS + monitoringSiteIdentifier + phenomenonTimeReferenceYear, data =Data.AggBySiteID.1, length)
table(EntriesPerCAS$resultMeanValue)
# 1      2      3      4      5      6      7      8      9
# 638102    211    162    167     54     21      8      5      3
# There are some repeated entries (multiple Averages for the same CAS, station and Year)
# Will take the maximum as worst-case scenario for the Average
EntriesPerCAS = EntriesPerCAS[EntriesPerCAS$resultMeanValue>1,]
unique(EntriesPerCAS$monitoringSiteIdentifier) #49 Sites from Italy with repeated entries

# Aggregate to keep the maximum
# @Zhenglei: is there a better way to keep only the entry with the maximum value for the Mean?
## Data.AggBySiteID.1 <- aggregate(. ~ monitoringSiteIdentifier + phenomenonTimeReferenceYear + CAS, data = Data.AggBySiteID.1, max)
Data.AggBySiteID.1 = Data.AggBySiteID.1 %>% group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,CAS) %>%
 summarise(k=which.max(resultMaximumValue),procedureLOQValue=procedureLOQValue[k],resultMeanValue=resultMeanValue[k],
           resultMaximumValue=resultMaximumValue[k],AboveLOQ=AboveLOQ[k],resultNumberOfSamples=resultNumberOfSamples[k],
           resultQualityNumberOfSamplesBelowLOQ=resultQualityNumberOfSamplesBelowLOQ[k])




###########################################################################
## RESTRICT THE DATASET BASED ON THE N CHEMICALS MEASURED PER STATION

# How many determinads measured per station and Year? ################################################
#variable to keep
var = c("monitoringSiteIdentifier","CAS","phenomenonTimeReferenceYear")
DataForNChem = Data.AggBySiteID.1[,names(Data.AggBySiteID.1)%in%c(var)]
measured.chem = aggregate(CAS ~ monitoringSiteIdentifier + phenomenonTimeReferenceYear, data =DataForNChem, length)
summary(measured.chem$CAS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00    2.00   22.00   31.46   45.00  191.00

# How many determinads measured per station and Year?
detected.chem = aggregate(CAS ~ monitoringSiteIdentifier + phenomenonTimeReferenceYear, data =DataForNChem[Data.AggBySiteID.1$resultNumberOfSamples>Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ,], length)
summary(detected.chem$CAS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   5.000   9.637  12.000 125.000

# Will only keep stations measuring at least 10 chemicals the same year
measured.chem = measured.chem[measured.chem$CAS>9,] #11668 StationID - Year combinations

summary(measured.chem$CAS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 10.00   18.00   36.00   42.94   53.00  271.00
summary(measured.chem$phenomenonTimeReferenceYear)
#1998 - 2018
length(unique(measured.chem$monitoringSiteIdentifier))
# 4609 Station IDs

#Subset data for these selection of Station Year combinations
Data.AggBySiteID.1$Site.Year = paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")
measured.chem$Site.Year = paste(measured.chem$monitoringSiteIdentifier,measured.chem$phenomenonTimeReferenceYear, sep=".")
detected.chem$Site.Year = paste(detected.chem$monitoringSiteIdentifier,detected.chem$phenomenonTimeReferenceYear, sep=".")

# Merge measured & detected
measured.chem = merge(measured.chem, detected.chem, by=c("monitoringSiteIdentifier","phenomenonTimeReferenceYear","Site.Year"), all.x = T)
names(measured.chem)[c(4,5)] = c("N.Measured", "N.Detected")
measured.chem$N.Detected[is.na(measured.chem$N.Detected)] = 0 # NAs are non-detected
summary(measured.chem)
hist(measured.chem$N.Measured)
hist(measured.chem$N.Detected)
measured.chem$Ratio.MD = measured.chem$N.Detected/measured.chem$N.Measured
summary(measured.chem$Ratio.MD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.06809 0.19048 0.24554 0.33333 1.00000

# Some summaries
table(measured.chem$N.Measured>9)
# TRUE
# 11668
table(measured.chem$N.Measured>50)
# FALSE  TRUE
# 8098  3570
table(measured.chem$N.Measured>100)
# FALSE  TRUE
# 10559  1109
table(measured.chem$N.Detected>10)
# FALSE  TRUE
# 7316  4352
table(measured.chem$N.Detected>50)
# FALSE  TRUE
# 11201   467

#Subset the original dataset by merging to measured chem
Data.AggBySiteID.2 = merge(Data.AggBySiteID.1,measured.chem, by=c("monitoringSiteIdentifier","phenomenonTimeReferenceYear","Site.Year"), all.y = T); Data.AggBySiteID.2 = droplevels(Data.AggBySiteID.2)
#571999 observations
length(unique(Data.AggBySiteID.2$monitoringSiteIdentifier)) #4609 Site IDs
length(unique(Data.AggBySiteID.2$CAS)) #327 CAS

#Remove useless variables from datase
usefull.var = c("monitoringSiteIdentifier","phenomenonTimeReferenceYear", "Site.Year", "resultMeanValue", "resultMaximumValue","CAS","N.Measured","N.Detected","procedureLOQValue", "AboveLOQ", "resultNumberOfSamples", "resultQualityNumberOfSamplesBelowLOQ")
Data.AggBySiteID.2 = Data.AggBySiteID.2[,names(Data.AggBySiteID.2)%in%usefull.var]

#########################################################################################
###### MAKE SURE YOU SELECT THE RIGHT LOQ TREATMENT OPTION BELOW BEFORE CONTINUING!!!!
#########################################################################################

# Options for samples below LOQ
# Option 1: BelowLOQ = 0
# Option 2: BelowLOW = 1/2 LOQ

# v2021: 619600     12
length((levels(Data.AggBySiteID.2$monitoringSiteIdentifier)))
length(unique(Data.AggBySiteID.2$monitoringSiteIdentifier))
length(unique(Stations$monitoringSiteIdentifier)) #  34522
stations <- Stations[Stations$monitoringSiteIdentifier %in% unique(Data.AggBySiteID.2$monitoringSiteIdentifier),]
dim(stations) # 4641x22
sum(complete.cases(stations[,c("lon","lat")])) ##4319

which(!complete.cases(stations[,c("lon","lat")]))

## Station summary:
usethis::use_data(stations,overwrite = T)

Res_Max_all <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),
                          Exclusion.CAS = c("Excluded_None","Refined_Metals_PAH","Excluded_Metals","Current_Use","Banned_Listed"),Data.SSD.1=Data.SSD.1,useAllMax=T)
Res_Mean_all <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),
                           Exclusion.CAS = c("Excluded_None","Refined_Metals_PAH","Excluded_Metals","Current_Use","Banned_Listed"),Data.SSD.1=Data.SSD.1,useAllMean=T)

Res <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),
                  Exclusion.CAS = c("Excluded_None","Refined_Metals_PAH","Excluded_Metals","Current_Use","Banned_Listed"),Data.SSD.1=Data.SSD.1)
stations1_Mean <- Res_Mean_all %>% group_by(LOQ.type,Exclusion.CAS,monitoringSiteIdentifier) %>%
  summarize(mean_HQ_Max_Acute=mean(HQ_Max_Acute,na.rm=T),mean_HQ_Max_Chronic=mean(HQ_Max_Chronic,na.rm=T),mean_HI_Acute=mean(HI_Acute,na.rm=T),mean_HI_Chronic=mean(HI_Chronic,na.rm=T),
            mean_N_measured=mean(Nmeasured,na.rm=T),mean_N_detected=mean(Ndetected,na.rm=T),mean_MCR_Chronic=mean(MCR.HC05.Chronic.NOEC,na.rm=T),mean_MCR_Acute=mean(MCR.HC50.Acute.EC50,na.rm=T))%>%
   mutate(Group_AF1_Chronic=get_group(HI=mean_HI_Chronic,HQ=mean_HQ_Max_Chronic,MCR=mean_MCR_Chronic,AF=1)) %>%
  mutate(Group_AF1_Acute=get_group(HI=mean_HI_Acute,HQ=mean_HQ_Max_Acute,MCR=mean_MCR_Acute,AF=1))
stations1_Max <- Res_Max_all %>% group_by(LOQ.type,Exclusion.CAS,monitoringSiteIdentifier) %>%
  summarize(mean_HQ_Max_Acute=mean(HQ_Max_Acute,na.rm=T),mean_HQ_Max_Chronic=mean(HQ_Max_Chronic,na.rm=T),mean_HI_Acute=mean(HI_Acute,na.rm=T),mean_HI_Chronic=mean(HI_Chronic,na.rm=T),
            mean_N_measured=mean(Nmeasured,na.rm=T),mean_N_detected=mean(Ndetected,na.rm=T),mean_MCR_Chronic=mean(MCR.HC05.Chronic.NOEC,na.rm=T),mean_MCR_Acute=mean(MCR.HC50.Acute.EC50,na.rm=T))%>%
  mutate(Group_AF1_Chronic=get_group(HI=mean_HI_Chronic,HQ=mean_HQ_Max_Chronic,MCR=mean_MCR_Chronic,AF=1)) %>%
  mutate(Group_AF1_Acute=get_group(HI=mean_HI_Acute,HQ=mean_HQ_Max_Acute,MCR=mean_MCR_Acute,AF=1))

stations1_Mean <- left_join(stations1_Mean,stations)
stations1_Max <- left_join(stations1_Max,stations)
stations1_Mean$Exclusion.CAS <- factor(stations1$Exclusion.CAS,levels=c("Excluded_None","Refined_Metals_PAH","Excluded_Metals","Banned_Listed","Current_Use"))
stations1_Max$Exclusion.CAS <- factor(stations1$Exclusion.CAS,levels=c("Excluded_None","Refined_Metals_PAH","Excluded_Metals","Banned_Listed","Current_Use"))
usethis::use_data(stations1_Max,overwrite = T)
usethis::use_data(stations1_Mean,overwrite = T)
############################################################
usethis::use_data(stations1,overwrite = T)



Res <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),
                  Exclusion.CAS = c("Excluded_None","Refined_Metals_PAH","Excluded_Metals_PAH","Current_Use"),Data.SSD.1=Data.SSD.1)
stations2 <- Res %>% group_by(LOQ.type,Exclusion.CAS,monitoringSiteIdentifier) %>%
  summarize(mean_HQ_Max_Acute=mean(HQ_Max_Acute,na.rm=T),mean_HQ_Max_Chronic=mean(HQ_Max_Chronic,na.rm=T),mean_HI_Acute=mean(HI_Acute,na.rm=T),mean_HI_Chronic=mean(HI_Chronic,na.rm=T),mean_N_measured=mean(Nmeasured,na.rm=T),mean_MCR_Chronic=mean(MCR.HC05.Chronic.NOEC,na.rm=T),mean_MCR_Acute=mean(MCR.HC50.Acute.EC50,na.rm=T))%>%
  mutate(Group_AF1_Chronic=get_group(HI=mean_HI_Chronic,HQ=mean_HQ_Max_Chronic,MCR=mean_MCR_Chronic,AF=1)) %>%
  mutate(Group_AF1_Acute=get_group(HI=mean_HI_Acute,HQ=mean_HQ_Max_Acute,MCR=mean_MCR_Acute,AF=1))
stations2 <- left_join(stations2,stations)
usethis::use_data(stations2,overwrite = T)

correctMetal <- F

if(correctMetal==T){
  ## Need to refine the concentrations of Metals!!
  metals <- chemclass %>% filter(Chem.Group=="Metal")
  all(metals$CAS %in% Data.AggBySiteID.2$CAS) ## should be TRUE
  Data.AggBySiteID.2 <- Data.AggBySiteID.2%>% left_join(metals[,c("CAS","Reference.Conc")])%>%
    mutate(Reference.Conc=replace_na(Reference.Conc,0)) %>% ##%>% filter(CAS %in% metals$CAS)
    mutate(resultMeanValue=resultMeanValue-Reference.Conc,resultMaximumValue=resultMaximumValue-Reference.Conc) ## %>% summarise(r1=sum(resultMeanValue<0),r2=sum(resultMaxValue<0))
  Data.AggBySiteID.2$resultMeanValue[Data.AggBySiteID.2$resultMeanValue<0] <- 0
  Data.AggBySiteID.2$resultMaximumValue[Data.AggBySiteID.2$resultMaximumValue<0] <- 0
  Data.AggBySiteID.2$Reference.Conc <- NULL
}

usethis::use_data(Data.AggBySiteID.2,overwrite = T)

rm(Data.AggBySiteID)
rm(Data.AggBySiteID.1)
rm(Stations)


donotrun <- T
if(!donotrun){
casinfo_site_year <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS,Site.Year) %>% nest()%>%
  mutate(AboveLOQ=any(data$AboveLOQ),SummaryLOQ=map(data,~summary(.x$procedureLOQValue)))%>%
  mutate(tidied = map(SummaryLOQ, tidy)) %>% unnest(tidied) %>%
  dplyr::select(Site.Year,CAS,AboveLOQ,minimum,q1, median,  mean, q3, maximum)
write.csv(casinfo_site_year,file="inst/manuscript_v3/cas_site_year_LOQ.csv")
}
