#' Function to check sensitivity to Sample and Chemical per year restrictions
#'
#' @param NSample.Year
#' @param Nchemicals.Year
#'
#' @return
#' @export
#'
#' @examples
data_prepare <- function(NSample.Year=5,Nchemicals.Year=10){



  ## ----setup-----------------------------------------------------------------------------------------------------
  library(WaterBase)




  #Define and set paths
  ## this need to be changed by the user if reproducing the base data in the package
  Dir_Data = "~/s3raw/WB-Zhenglei/"


  library(tidyverse)
  # Load Endpoint data
  #L.Posthuma SSDs
  Data.SSD <- read.csv(paste0(Dir_Data,"Copy of etc4373-sup-0002-supmat.csv"))

  #Load data. Basic exploration of stations information

  Stations = read.csv(paste0(Dir_Data,"Waterbase_v2020_1_S_WISE6_SpatialObject_DerivedData.csv")) # v2021 ==> too many missing information use the previous one instead.
  Stations = read.csv(paste0(Dir_Data,"Waterbase_v2019_1_S_WISE6_SpatialObject_DerivedData_repaired.csv")) # v2021
  #How many Stations?
  Stations$monitoringSiteIdentifier = as.factor(Stations$monitoringSiteIdentifier)
  length(levels(Stations$monitoringSiteIdentifier))
  #60776
  #How many have complete lat/long info?
  sum(complete.cases(Stations[names(Stations)%in%c("lon", "lat")]))
  #48563


  ## --------------------------------------------------------------------------------------------------------------
  library(data.table)
  Data.AggBySiteID = fread(paste0(Dir_Data,"Waterbase_v2020_1_T_WISE6_AggregatedData.csv"))
  ## Data.AggBySiteID = read.csv(paste0(Dir_Data,"Waterbase_v2020_1_T_WISE6_AggregatedData.csv")) #2021 #3962763 rows
  dim(Data.AggBySiteID)
  # [1] 3962763      31
  Data.AggBySiteID <- left_join(Data.AggBySiteID,Stations[,c("monitoringSiteIdentifier","waterBodyIdentifier","countryCode")])
  length(unique(Data.AggBySiteID$monitoringSiteIdentifier))
  length(unique(cbind(Data.AggBySiteID$monitoringSiteIdentifier,Data.AggBySiteID$phenomenonTimeReferenceYear)))
  length(unique(Data.AggBySiteID$observedPropertyDeterminandCode))
  data_prep <- data.frame(Cat="Data Selection",Step=1,N.Sample=3962763,N.CAS=552,N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=25974,N.SiteYear=330096)


  ## --------------------------------------------------------------------------------------------------------------
  names(Data.AggBySiteID)[1] = "monitoringSiteIdentifier"

  Data.AggBySiteID$monitoringSiteIdentifier = as.factor(Data.AggBySiteID$monitoringSiteIdentifier)
  # Prepare CAS numbers
  Data.AggBySiteID$CAS = gsub("CAS_|-|_", "", Data.AggBySiteID$observedPropertyDeterminandCode)

  #How many CAS numbers in SSDs?
  length(unique(Data.AggBySiteID$CAS)) #552 determinads



  ## --------------------------------------------------------------------------------------------------------------
  # Keep only "W" to avoid replications
  Data.AggBySiteID = Data.AggBySiteID[Data.AggBySiteID$procedureAnalysedMatrix=="W",]
  dim(Data.AggBySiteID)

  data_prep <- rbind(data_prep,data.frame(Cat="Data Selection",Step=2,N.Sample=3469602,N.CAS=length(unique(Data.AggBySiteID$CAS)),N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID$monitoringSiteIdentifier)),N.SiteYear=length(unique(paste(Data.AggBySiteID$monitoringSiteIdentifier,Data.AggBySiteID$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  table(Data.AggBySiteID$resultObservationStatus)

  #              A
  #3085822  383780
  ## Note in the current dataset there is no "O" result observation status
  ## Data.AggBySiteID =Data.AggBySiteID[Data.AggBySiteID$resultObservationStatus!="O",]

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

  dim(Data.AggBySiteID)

  data_prep <- rbind(data_prep,data.frame(Cat="Data Selection",Step=3,N.Sample=3343206,N.CAS=length(unique(Data.AggBySiteID$CAS)),N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID$monitoringSiteIdentifier)),N.SiteYear=length(unique(paste(Data.AggBySiteID$monitoringSiteIdentifier,Data.AggBySiteID$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  unique(Data.AggBySiteID$metadata_statements)
  table(Data.AggBySiteID$metadata_statements!="NOTE_LEGACY: Measurement has been confirmed by country to be taken from a highly polluted area")
  # FALSE    TRUE
  # 85 3343121
  #Remove 85 because MAF is not to cover hot spot contamination
  Data.AggBySiteID =Data.AggBySiteID[Data.AggBySiteID$metadata_statements!="NOTE_LEGACY: Measurement has been confirmed by country to be taken from a highly polluted area",]


  dim(Data.AggBySiteID)

  data_prep <- rbind(data_prep,data.frame(Cat="Data Selection",Step=4,N.Sample=3343121,N.CAS=length(unique(Data.AggBySiteID$CAS)),N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID$monitoringSiteIdentifier)),N.SiteYear=length(unique(paste(Data.AggBySiteID$monitoringSiteIdentifier,Data.AggBySiteID$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
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

  dim(Data.AggBySiteID)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Selection",Step=5,N.Sample=2989996,N.CAS=length(unique(Data.AggBySiteID$CAS)),N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID$monitoringSiteIdentifier)),N.SiteYear=length(unique(paste(Data.AggBySiteID$monitoringSiteIdentifier,Data.AggBySiteID$phenomenonTimeReferenceYear, sep=".")))))

  # browser()

  ## --------------------------------------------------------------------------------------------------------------
  # How many entries from the existing data we will be able to use?

  id_del <- grepl("EEA_",(Data.AggBySiteID$observedPropertyDeterminandCode),ignore.case = T)
  unique(Data.AggBySiteID$observedPropertyDeterminandCode[id_del])
  Data.AggBySiteID <- Data.AggBySiteID[!id_del,]


  dim(Data.AggBySiteID)

  data_prep <- rbind(data_prep,data.frame(Cat="Data Selection",Step=6,N.Sample=nrow(Data.AggBySiteID),N.CAS=length(unique(Data.AggBySiteID$CAS)),N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID$monitoringSiteIdentifier,Data.AggBySiteID$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------



  ## --------------------------------------------------------------------------------------------------------------

  CAS.SSD = unique(Data.AggBySiteID$CAS)[unique(Data.AggBySiteID$CAS)%in%as.character(unique(Data.SSD$CAS.))]
  length(unique(CAS.SSD)) #373

  Data.AggBySiteID.1 = droplevels(Data.AggBySiteID) #[1] 1816359      32
  Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID$CAS%in%CAS.SSD,]

  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=1,N.Sample=nrow(Data.AggBySiteID.1),N.CAS=length(unique(Data.AggBySiteID$CAS)),N.Country=length(unique(Data.AggBySiteID$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))



  ## --------------------------------------------------------------------------------------------------------------


  Data.SSD.1 = Data.SSD[Data.SSD$CAS.%in%CAS.SSD,]; Data.SSD.1 = droplevels(Data.SSD.1) #373x38

  table(Data.SSD.1$X.Species.Acute.EC50>=3)
  # FALSE  TRUE
  # 22   348
  # Remove chem that have less than 3 sp per SSD
  Data.SSD.1 = Data.SSD.1[Data.SSD.1$X.Species.Acute.EC50>=3,]
  Data.SSD.1 = Data.SSD.1[Data.SSD.1$X.Species.Chronic.NOEC>=3,]
  Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X.Species.Chronic.NOEC),]
  Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X.Species.Acute.EC50),]
  dim(Data.SSD.1)
  ##[1] 346  38

  # Data quality
  table(Data.SSD.1$OrigenQuality.Chronic.NOEC)
  table(Data.SSD.1$OrigenQuality.Acute.EC50)

  # Highest data quality in all cases
  summary(Data.SSD.1$X.Species.Acute.EC50)
  summary(Data.SSD.1$X.Species.Chronic.NOEC)
  summary(Data.SSD.1$X.TaxClass.Acute.EC50)

  # Taxa clasess
  # Different taxa classes, need to evaluate what exactly is this.
  table(Data.SSD.1$X.TaxClass.Acute.EC50>=3)
  # FALSE  TRUE
  # 39   307
  table(Data.SSD.1$X.TaxClass.Chronic.NOEC>=3)
  # FALSE  TRUE
  # 61   285


  # SSD data with Slope data
  Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X10LogSSDSlope.ug.L..SigmaChronic.NOEC),]
  Data.SSD.1 = Data.SSD.1[!is.na(Data.SSD.1$X10LogSSDSlope.ug.L..SigmaAcute.EC50),]


  # Final interim set of chemicals for analysis
  CAS.SSD.1 = CAS.SSD[CAS.SSD%in%Data.SSD.1$CAS.]
  Data.SSD.1 = Data.SSD.1[Data.SSD.1$CAS.%in%CAS.SSD.1,]; Data.SSD.1 = droplevels(Data.SSD.1)


  ## --------------------------------------------------------------------------------------------------------------
  #########################################################################
  # Subset Water-base rivers based on available SSDs

  Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$CAS%in%CAS.SSD.1,];
  Data.AggBySiteID.1 = droplevels(Data.AggBySiteID.1) #1679005 entries [1] 1506907      32

  dim(Data.AggBySiteID.1)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=2,N.Sample=1506907,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  # Units
  table(Data.AggBySiteID.1$resultUom)
  ## mg{NH4}/L   mg{P}/L   mg{S}/L      mg/L      ug/L
  ##    102903    186232        10     97814   1429400

  ## mg{P}/L mg{S}/L    mg/L    ug/L
  ##   91727      10   22611 1392559

  #Remove all but mg/L and ug/L and transform all to mg/L
  Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$resultUom%in%c("mg/L", "ug/L"),] #1469931 # 1415170
  dim(Data.AggBySiteID.1) ## 1415170 32

  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=3,N.Sample=1415170,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  #transform to ug/L
  Data.AggBySiteID.1$resultMinimumValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMinimumValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
  Data.AggBySiteID.1$resultMeanValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMeanValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
  Data.AggBySiteID.1$resultMedianValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMedianValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
  Data.AggBySiteID.1$resultMaximumValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$resultMaximumValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
  # NOTE: Unit label is changed later after removing unreliable LOQs

  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=4,N.Sample=1415170,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  # Subset for complete data for Average and Maximun values
  ## 43 missing mean and 26681 missing maximum
  Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultMeanValue),]
  Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultMaximumValue),]

  dim(Data.AggBySiteID.1)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=5,N.Sample=1388489,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear,sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  # Remove stations do not reporting number of samples #47
  Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultNumberOfSamples),]

  # Stations reporting more than 5 samples per year to report annual aggregates
  Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$resultNumberOfSamples>=NSample.Year,] #962409 v2021


  dim(Data.AggBySiteID.1)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=6,N.Sample=962409,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  ###########################################################################
  # LOQs: Remove unreliable entries based on missing/unreliable LOQs reporting

  # Remove data for which LOQ is not reported ## 203610 were removed
  Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$procedureLOQValue),]

  # Remove data with LOQ <= 0 # 7113 were removed
  Data.AggBySiteID.1 = Data.AggBySiteID.1[Data.AggBySiteID.1$procedureLOQValue> 0,]

  # N Samples below LOQ
  summary(Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ)

  # Remove data that do not report "resultQualityNumberOfSamplesBelowLOQ" ## 114452 were removed
  Data.AggBySiteID.1 = Data.AggBySiteID.1[!is.na(Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ),] #640197 #[1] 637234

  #N samples below LOQ
  summary(Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
  #    0.00    6.00   11.00    9.24   12.00  366.00  283332
  ##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  ##   0.000   6.000  11.000   9.678  12.000 366.000

  # Convert reported LOQs to ug/L
  Data.AggBySiteID.1$procedureLOQValue[Data.AggBySiteID.1$resultUom=="mg/L"] = Data.AggBySiteID.1$procedureLOQValue[Data.AggBySiteID.1$resultUom=="mg/L"]*1000 # to ug/L
  # Change unit label
  Data.AggBySiteID.1$resultUom[Data.AggBySiteID.1$resultUom=="mg/L"] = "ug/L"; Data.AggBySiteID.1 = droplevels(Data.AggBySiteID.1)

  # Number of samples above LOQ
  table(Data.AggBySiteID.1$resultNumberOfSamples > Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ)
  #73.9% of samples below LOQ; 30% actual measurements
  # This is very important to see what to do with BelowLOQ values!

  #if N samples > 1, and N.detects > 2, then max > mean. Therefore, max > mean can be used to separate the samples with detects from those with non-detects
  table(Data.AggBySiteID.1$resultMaximumValue>Data.AggBySiteID.1$resultMeanValue)

  # 73% below LOQ, 37% above LOQ
  # The number should be equal to the value in the previous logic, which is not!




  #2. Set samples with Max <= LOQ to below LOQ
  Data.AggBySiteID.1$AboveLOQ = Data.AggBySiteID.1$resultNumberOfSamples > Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ

  Data.AggBySiteID.1$AboveLOQ[Data.AggBySiteID.1$resultMaximumValue<=Data.AggBySiteID.1$procedureLOQValue] = FALSE

  summary(Data.AggBySiteID.1$AboveLOQ)

  Data.AggBySiteID.1$AboveLOQ[Data.AggBySiteID.1$resultMaximumValue<Data.AggBySiteID.1$resultMeanValue]
  Data.AggBySiteID.1$AboveLOQ[Data.AggBySiteID.1$resultMaximumValue<=Data.AggBySiteID.1$procedureLOQValue]  ## 479812

  Data.AggBySiteID.1$AboveLOQ[Data.AggBySiteID.1$resultNumberOfSamples <= Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ]


  ## --------------------------------------------------------------------------------------------------------------
  dim(Data.AggBySiteID.1)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=7,N.Sample=637234,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------



  # Aggregate to keep the maximum
  # @Zhenglei: is there a better way to keep only the entry with the maximum value for the Mean?
  ## Data.AggBySiteID.1 <- aggregate(. ~ monitoringSiteIdentifier + phenomenonTimeReferenceYear + CAS, data = Data.AggBySiteID.1, max)
  Data.AggBySiteID.1 = Data.AggBySiteID.1 %>% group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,CAS) %>%
    summarise(k=which.max(resultMaximumValue),procedureLOQValue=procedureLOQValue[k],resultMeanValue=resultMeanValue[k],
              resultMaximumValue=resultMaximumValue[k],AboveLOQ=AboveLOQ[k],resultNumberOfSamples=resultNumberOfSamples[k],
              resultQualityNumberOfSamplesBelowLOQ=resultQualityNumberOfSamplesBelowLOQ[k],waterBodyIdentifier=waterBodyIdentifier[k],countryCode=countryCode[k])




  dim(Data.AggBySiteID.1)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=8,N.Sample=635771,N.CAS=length(unique(Data.AggBySiteID.1$CAS)),N.Country=length(unique(Data.AggBySiteID.1$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.1$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.1$monitoringSiteIdentifier)),
                                          N.SiteYear=length(unique(paste(Data.AggBySiteID.1$monitoringSiteIdentifier,Data.AggBySiteID.1$phenomenonTimeReferenceYear, sep=".")))))


  ## --------------------------------------------------------------------------------------------------------------
  ###########################################################################
  ## RESTRICT THE DATASET BASED ON THE N CHEMICALS MEASURED PER STATION

  # How many determinads measured per station and Year? ################################################
  #variable to keep
  var = c("monitoringSiteIdentifier","CAS","phenomenonTimeReferenceYear")
  DataForNChem = Data.AggBySiteID.1[,names(Data.AggBySiteID.1)%in%c(var)]
  # browser()
  measured.chem = aggregate(CAS ~ monitoringSiteIdentifier + phenomenonTimeReferenceYear, data =DataForNChem, length)
  summary(measured.chem$CAS)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 1.00    2.00   22.00   31.46   45.00  191.00
  ## Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  ##   1.00    2.00   22.00   32.62   45.00  192.00
  tmp <- Data.AggBySiteID.1%>% group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear)%>%summarise(Nmeasured=length(CAS))
  summary(tmp$Nmeasured) ## validated Ismael's original code

  # How many determinads measured per station and Year?
  detected.chem = aggregate(CAS ~ monitoringSiteIdentifier + phenomenonTimeReferenceYear, data =DataForNChem[Data.AggBySiteID.1$resultNumberOfSamples>Data.AggBySiteID.1$resultQualityNumberOfSamplesBelowLOQ,], length)
  summary(detected.chem$CAS)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 1.000   1.000   5.000   9.637  12.000 125.000

  # Will only keep stations measuring at least 10 chemicals the same year
  measured.chem = measured.chem[measured.chem$CAS>(Nchemicals.Year-1),] #11668 StationID - Year combinations

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
  dim(Data.AggBySiteID.1)
  #Subset the original dataset by merging to measured chem
  Data.AggBySiteID.2 = merge(Data.AggBySiteID.1,measured.chem, by=c("monitoringSiteIdentifier","phenomenonTimeReferenceYear","Site.Year"), all.y = T); Data.AggBySiteID.2 = droplevels(Data.AggBySiteID.2)

  #571999 observations
  length(unique(Data.AggBySiteID.2$monitoringSiteIdentifier)) #4664 Site IDs
  length(unique(Data.AggBySiteID.2$CAS)) #334 CAS


  dim(Data.AggBySiteID.2)
  data_prep <- rbind(data_prep,data.frame(Cat="Data Curation",Step=9,N.Sample=nrow(Data.AggBySiteID.2),N.CAS=length(unique(Data.AggBySiteID.2$CAS)),N.Country=length(unique(Data.AggBySiteID.2$countryCode))-1,N.WaterBody=length(unique(Data.AggBySiteID.2$waterBodyIdentifier))-1,N.Site=length(unique(Data.AggBySiteID.2$monitoringSiteIdentifier)),
                                          N.SiteYear= length(unique(paste(Data.AggBySiteID.2$monitoringSiteIdentifier,Data.AggBySiteID.2$phenomenonTimeReferenceYear, sep=".")))))

  #Remove useless variables from datase
  usefull.var = c("monitoringSiteIdentifier","phenomenonTimeReferenceYear", "Site.Year", "resultMeanValue", "resultMaximumValue","CAS","N.Measured","N.Detected","procedureLOQValue", "AboveLOQ", "resultNumberOfSamples", "resultQualityNumberOfSamplesBelowLOQ")
  Data.AggBySiteID.2 = Data.AggBySiteID.2[,names(Data.AggBySiteID.2)%in%usefull.var]










  return(list(Data.AggBySiteID.2=Data.AggBySiteID.2,data_prep=data_prep))
}

#' Calculate the CEFIC_MIAT table
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param LOQ.type
#' @param Exclusion.CAS
#'
#' @return
#' @export
#'
#' @examples
getCEFIC_MIAT <- function(Data.AggBySiteID.2,chemclass,LOQ.type="LOQ.T0",Exclusion.CAS="Excluded_Metals" ){
  tmp1 <- Data.AggBySiteID.2 %>% group_by(Site.Year) %>% summarise(N=length(Site.Year)) %>% mutate(n1=mean(N))
  data("Data.SSD")
  #data("chemclass")
  CAS.SSD = unique(Data.AggBySiteID.2$CAS)[unique(Data.AggBySiteID.2$CAS)%in%as.character(unique(Data.SSD$CAS.))]
  Data.SSD.1 = Data.SSD[Data.SSD$CAS.%in%CAS.SSD,]; Data.SSD.1 = droplevels(Data.SSD.1)
  Res_Max <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                        Exclusion.CAS = Exclusion.CAS ,Data.SSD.1=Data.SSD.1,useAllMax=T)
  Res_Mean <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                         Exclusion.CAS = Exclusion.CAS ,Data.SSD.1=Data.SSD.1,useAllMean=T)
  N_CAS_Site_tmp <- Res_Mean %>%
    group_by(LOQ.type,Exclusion.CAS) %>%
    summarise(N.Site=length(unique(monitoringSiteIdentifier)),
              N.Sample=length(monitoringSiteIdentifier)) #%>%
   # mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))

  table_chronic_mean <- Res_Mean%>%
    group_by(LOQ.type,Exclusion.CAS,Group_AF1_Chronic) %>%
    summarise(N=length(monitoringSiteIdentifier))%>%
    ungroup %>%
    left_join(N_CAS_Site_tmp[,c("LOQ.type","Exclusion.CAS","N.Sample")])%>%
    mutate(Perc.Mean=N/N.Sample*100)%>%
    mutate(N.Site.Mean=N)%>%
    select(c(LOQ.type,Exclusion.CAS,Group_AF1_Chronic,N.Sample,N.Site.Mean,Perc.Mean))
  table_chronic_max <- Res_Max%>% group_by(LOQ.type,Exclusion.CAS,Group_AF1_Chronic) %>%
    summarise(N=length(monitoringSiteIdentifier))%>% ungroup %>%
    left_join(N_CAS_Site_tmp[,c("LOQ.type","Exclusion.CAS","N.Sample")])%>% mutate(Perc.Max=N/N.Sample*100)%>% mutate(N.Site.Max=N)%>% select(c(LOQ.type,Exclusion.CAS,Group_AF1_Chronic,N.Sample,N.Site.Max,Perc.Max))

  table_chronic <-left_join(table_chronic_mean,table_chronic_max)
  ## attr(table_chronic$Perc.Max, "format") <- "%.2f"
  ## fdata(table_chronic)
  table_chronic$Perc.Max <- scales::percent(table_chronic$Perc.Max/100, accuracy = 0.1)
  table_chronic$Perc.Mean <- scales::percent(table_chronic$Perc.Mean/100, accuracy = 0.1)
  return(table_chronic)
}
