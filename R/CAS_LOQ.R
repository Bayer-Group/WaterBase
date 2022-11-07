#' select the designed datasets according to LOQ settings and CAS exclusions
#'
#' @param Data.AggBySiteID.2 data stored in the package, after minimum preprocessing
#' @param chemclass chemical information data
#' @param LOQ.type "LOQ.T0" or "LOQ.T1"
#' @param Data.SSD.1 SSD information
#' @param factor.LOQ use 1 or 1/2 as LOQ handling method
#' @param LOQ Always use TRUE
#' @param SSDrefine whether to use refined SSD information,
#' if excluded none, then use original waterbase data with original SSD. otherwise, used PAH, metals refinement
#' @param Exclusion.CAS exclusion
#' @param Exclusion Keep it as TRUE
#' @param HU whether to calculate HU
#' @param PAF whether to calculate PAF
#'
#' @return
#' @export
#'
#' @examples
preprocessWaterbase <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),Data.SSD.1,
                                Exclusion.CAS = c("Excluded_None", "Excluded_Metals","Excluded_Metals_PAH", "Current_Use", "Metals", "Banned_Listed","Refined_HQ_Metals_PAH"),
                                factor.LOQ=1,
                                LOQ=T, ## Do not change this to be FALSE!
                                SSDrefine=T,
                                Exclusion=T,
                                HU=T,
                                PAF=T,
                                useAllMax=F,
                                useAllMean=F

){

  if(Exclusion.CAS=="Excluded_None"){
    ## do not refine the PAH and Metal???
    SSDrefine <- F
  }else{
    ## NEED to refine Metals Concentrations first!!!!
    ## Need to refine the concentrations of Metals!!
    metals <- chemclass %>% filter(Chem.Group=="Metal")
    all(metals$CAS %in% Data.AggBySiteID.2$CAS) ## should be TRUE
    Data.AggBySiteID.2 <- Data.AggBySiteID.2%>% left_join(metals[,c("CAS","Reference.Conc")])%>%
      mutate(Reference.Conc=replace_na(Reference.Conc,0)) %>% ##%>% filter(CAS %in% metals$CAS)
      mutate(resultMeanValue=resultMeanValue-Reference.Conc,resultMaximumValue=resultMaximumValue-Reference.Conc) ## %>% summarise(r1=sum(resultMeanValue<0),r2=sum(resultMaxValue<0))
    Data.AggBySiteID.2$resultMeanValue[Data.AggBySiteID.2$resultMeanValue<0] <- 0
    Data.AggBySiteID.2$resultMaximumValue[Data.AggBySiteID.2$resultMaximumValue<0] <- 0
    Data.AggBySiteID.2$Reference.Conc <- NULL
    SSDrefine <- T
  }
  Data.AggBySiteID.3 <-  Data.AggBySiteID.2

  if(LOQ){
    if(LOQ.type == "LOQ.T0"){
      Data.AggBySiteID.3$resultMaximumValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = 0
      Data.AggBySiteID.3$resultMeanValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = 0
    }else{
      Data.AggBySiteID.3 = Data.AggBySiteID.2
      Data.AggBySiteID.3$resultMaximumValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = factor.LOQ*Data.AggBySiteID.3$procedureLOQValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"]
      Data.AggBySiteID.3$resultMeanValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = factor.LOQ*Data.AggBySiteID.3$procedureLOQValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"]
    }
    #Data.AggBySiteID.3a <- Data.AggBySiteID.3
  }





  ## Read in ssd information
  ssdInfo <- Data.SSD.1[,c("CAS.","X10LogSSDMedianConcentration.ug.L..MuAcute.EC50","X10LogSSDMedianConcentration.ug.L..MuChronic.NOEC", "X10LogSSDSlope.ug.L..SigmaAcute.EC50","X10LogSSDSlope.ug.L..SigmaChronic.NOEC")] %>% mutate(CAS.=as.character(CAS.))
  Data.AggBySiteID.3a <- left_join(Data.AggBySiteID.3,ssdInfo,by=c("CAS"="CAS."))


  #names(Data.AggBySiteID.3a)[(13:16)] = c("SSDLOG10.Mu.Acute.EC50", "SSDLOG10.Mu.chronic.NOEC", "SSDLOG10.Sigma.Acute.EC50","SSDLOG10.Sigma.chronic.NOEC")
  data.table::setnames(Data.AggBySiteID.3a,old =c("X10LogSSDMedianConcentration.ug.L..MuAcute.EC50","X10LogSSDMedianConcentration.ug.L..MuChronic.NOEC", "X10LogSSDSlope.ug.L..SigmaAcute.EC50","X10LogSSDSlope.ug.L..SigmaChronic.NOEC"),
                       new=c("SSDLOG10.Mu.Acute.EC50", "SSDLOG10.Mu.chronic.NOEC", "SSDLOG10.Sigma.Acute.EC50","SSDLOG10.Sigma.chronic.NOEC")
  )
  #Hazard.Index to the midpoint SSD by chem (reference to the HC05)
  #Hazard.Index to the HC05 SSD by chem
  # Taken from: https://edild.github.io/ssd/
  # Will use compoud-specific slopes. Assuming 0.7 SSD slope is quite worst case (as how the data comes out)

  #Compute HC05 by chemical ############################################################################
  #The mu and sig parameters are given in log10 units!
  Data.AggBySiteID.3a$HC50.Acute = 10^qnorm(0.5,Data.AggBySiteID.3a$SSDLOG10.Mu.Acute.EC50, Data.AggBySiteID.3a$SSDLOG10.Sigma.Acute.EC50)
  Data.AggBySiteID.3a$HC05.Chronic = 10^qnorm(0.05,Data.AggBySiteID.3a$SSDLOG10.Mu.chronic.NOEC, Data.AggBySiteID.3a$SSDLOG10.Sigma.chronic.NOEC)
  if(SSDrefine){ ## only if Excluded none, then no refinement on PAH.
    ## Need to refine PAH HC05.Chronic using PNEC!!!
    ## all
    PAHs <- chemclass %>% filter(Chem.Group=="PAH") %>% filter(!is.na(PNEC))
    ##  all(PAHs$CAS %in% Data.AggBySiteID.2$CAS) ## should be TRUE
    Data.AggBySiteID.3a <- Data.AggBySiteID.3a%>% left_join(PAHs[,c("CAS","PNEC")])
    ## plot(density(Data.AggBySiteID.3a$HC05.Chronic[!is.na(Data.AggBySiteID.3a$PNEC)]/Data.AggBySiteID.3a$PNEC[!is.na(Data.AggBySiteID.3a$PNEC)]))
    ##summary((Data.AggBySiteID.3a$HC05.Chronic[!is.na(Data.AggBySiteID.3a$PNEC)]/Data.AggBySiteID.3a$PNEC[!is.na(Data.AggBySiteID.3a$PNEC)]))
    # sum((Data.AggBySiteID.3a$HC05.Chronic[!is.na(Data.AggBySiteID.3a$PNEC)]/Data.AggBySiteID.3a$PNEC[!is.na(Data.AggBySiteID.3a$PNEC)])>1)
    Data.AggBySiteID.3a$HC05.Chronic[!is.na(Data.AggBySiteID.3a$PNEC)] <-Data.AggBySiteID.3a$PNEC[!is.na(Data.AggBySiteID.3a$PNEC)]
    Data.AggBySiteID.3a$PNEC <- NULL
  }
  if(!useAllMean){
    Data.AggBySiteID.3a$HQ.Acute = Data.AggBySiteID.3a$resultMaximumValue/Data.AggBySiteID.3a$HC50.Acute
  }else{
    Data.AggBySiteID.3a$HQ.Acute = Data.AggBySiteID.3a$resultMeanValue/Data.AggBySiteID.3a$HC50.Acute
  }
  if(!useAllMax){
    Data.AggBySiteID.3a$HQ.Chronic = Data.AggBySiteID.3a$resultMeanValue/Data.AggBySiteID.3a$HC05.Chronic
  }else{
    Data.AggBySiteID.3a$HQ.Chronic = Data.AggBySiteID.3a$resultMaximumValue/Data.AggBySiteID.3a$HC05.Chronic
  }

  # PAF by chemical ##################################################################################
  # Calculate HU (Hazard Unit) Anti-log10 are taken on Mu because HUs needs to be calculated in log-non-tramsformed data.
  if(HU){
    if(!useAllMean){
      Data.AggBySiteID.3a$HU.Acute.EC50 = Data.AggBySiteID.3a$resultMaximumValue/(10^Data.AggBySiteID.3a$SSDLOG10.Mu.Acute.EC50)
    }else{
      Data.AggBySiteID.3a$HU.Acute.EC50 = Data.AggBySiteID.3a$resultMeanValue/(10^Data.AggBySiteID.3a$SSDLOG10.Mu.Acute.EC50)

    }
    if(!useAllMax){
      Data.AggBySiteID.3a$HU.Chronic.NOEC = Data.AggBySiteID.3a$resultMeanValue/(10^Data.AggBySiteID.3a$SSDLOG10.Mu.chronic.NOEC)
    }else{
      Data.AggBySiteID.3a$HU.Chronic.NOEC = Data.AggBySiteID.3a$resultMaximumValue/(10^Data.AggBySiteID.3a$SSDLOG10.Mu.chronic.NOEC)

    }
  }
  if(PAF){
    # Calculate PAF
    Data.AggBySiteID.3a$PAF.Acute.EC50 = pnorm(log10(Data.AggBySiteID.3a$HU.Acute.EC50),0,Data.AggBySiteID.3a$SSDLOG10.Sigma.Acute.EC50,log=F) ## Note this is equal to pnorm(log10(Data.AggBySiteID.3a$resultMaximumValue),Data.AggBySiteID.3a$SSDLOG10.Mu.Acute.EC50 SSDLOG10.Sigma.Acute.EC50,)
    # del <- Data.AggBySiteID.3a %>% mutate(PAF.Acute= pnorm(log10(resultMaximumValue),SSDLOG10.Mu.Acute.EC50,SSDLOG10.Sigma.Acute.EC50))
    Data.AggBySiteID.3a$PAF.Chronic.NOEC = pnorm(log10(Data.AggBySiteID.3a$HU.Chronic.NOEC),0,Data.AggBySiteID.3a$SSDLOG10.Sigma.chronic.NOEC,log=F)
  }
  ## Reassign dataset 3.
  Data.AggBySiteID.3 <- Data.AggBySiteID.3a
  ## rm("Data.AggBySiteID.3a")

  if(Exclusion){
    require(tidyverse)
    if(Exclusion.CAS=="Excluded_Metals"){
      cas <- chemclass %>% filter(Chem.Group!="Metal")
    }else if(Exclusion.CAS=="Current_Use"){
      cas <- chemclass %>% filter(class=="Current_Use")
    }else if(Exclusion.CAS=="Metals"){
      cas <- chemclass %>% filter(Chem.Group=="Metal")
    }else if(Exclusion.CAS=="PAH"){
      cas <- chemclass %>% filter(Chem.Group=="PAH")
    }else if(Exclusion.CAS=="Banned_Listed"){
      #cas1 <- chemclass %>% filter(Chem.Group!="Metal" & Approved=="Y" & Chem.Group!="PAH"& Priority=="")
      #cas <-  chemclass %>% filter(Chem.Group!="Metal")%>% filter(!(CAS %in% cas1$CAS))
      cas <-  chemclass %>% filter(class=="Banned_Listed")
    }else if(Exclusion.CAS=="Excluded_Metals_PAH"){
      cas <-  chemclass %>% filter(Chem.Group!="Metal" & Chem.Group!="PAH")
    }else if(Exclusion.CAS=="Refined_HQ_Metals_PAH"){
      cas <-  chemclass
    }else{
      cas <- chemclass
    }
    Data.AggBySiteID.3 <- Data.AggBySiteID.3 %>% filter(CAS %in% cas$CAS)
  }
  print(dim(Data.AggBySiteID.3)) ## [1] 616795     16
  return(Data.AggBySiteID.3)
}

#' Deprecated: select the designed datasets according to LOQ settings and CAS exclusions
#'
#' @param Data.AggBySiteID.2 data stored in the package, after minimum preprocessing
#' @param chemclass chemical info
#' @param LOQ.type T0 or T1
#' @param Exclusion.CAS exclusion CAS choices
#'
#' @return
#' @export
#'
#' @examples
selectdata <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),
                       Exclusion.CAS = c("Excluded_None", "Excluded_Metals","Excluded_Metals_PAH", "Current_Use", "Metals", "Banned_Listed","Refined_HQ_Metals_PAH"),
                       factor.LOQ=1
)
{
  if(LOQ.type == "LOQ.T0"){
    Data.AggBySiteID.3 = Data.AggBySiteID.2
    Data.AggBySiteID.3$resultMaximumValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = 0
    Data.AggBySiteID.3$resultMeanValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = 0
  }else{
    Data.AggBySiteID.3 = Data.AggBySiteID.2
    Data.AggBySiteID.3$resultMaximumValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = factor.LOQ*Data.AggBySiteID.3$procedureLOQValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"]
    Data.AggBySiteID.3$resultMeanValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"] = factor.LOQ*Data.AggBySiteID.3$procedureLOQValue[Data.AggBySiteID.3$AboveLOQ=="FALSE"]
  }
  require(tidyverse)
  if(Exclusion.CAS=="Excluded_Metals"){
    cas <- chemclass %>% filter(Chem.Group!="Metal")
  }else if(Exclusion.CAS=="Current_Use"){
    cas <- chemclass %>% filter(Chem.Group!="Metal" & Approved=="Y" & Chem.Group!="PAH"& Priority=="")
  }else if(Exclusion.CAS=="Metals"){
    cas <- chemclass %>% filter(Chem.Group=="Metal")
  }else if(Exclusion.CAS=="Banned_Listed"){
    #cas1 <- chemclass %>% filter(Chem.Group!="Metal" & Approved=="Y" & Chem.Group!="PAH"& Priority=="")
    #cas <-  chemclass %>% filter(Chem.Group!="Metal")%>% filter(!(CAS %in% cas1$CAS))
    cas <-  chemclass %>% filter(Chem.Group!="Metal")%>% filter(!( Approved=="Y" & Chem.Group!="PAH"& Priority==""))
  }else if(Exclusion.CAS=="Excluded_Metals_PAH"){
    cas <-  chemclass %>% filter(Chem.Group!="Metal" & Chem.Group!="PAH")
  }else if(Exclusion.CAS=="Refined_HQ_Metals_PAH"){
    cas <-  chemclass
  }else{
    cas <- chemclass
  }
  Data.AggBySiteID.3 <- Data.AggBySiteID.3 %>% filter(CAS %in% cas$CAS)
  return(Data.AggBySiteID.3)
}


#' get the class of the chemicals
#'
#' @param Chem.Group chemical group
#' @param Approved whether being approved
#' @param Priority "STC", "WFD", "STC, WFD"
#' @param Candidate.SVHC
#' @param SVHC
#'
#' @return class being either Metal, Current_Use, or Banned_Listed.
#' @export
#'
#' @examples
getclass_deprecated <- function(Chem.Group,Approved,Priority,Candidate.SVHC,SVHC){
  classes <- sapply(1:length(Chem.Group),function(i)
  {
    if(Chem.Group[i]=="Metal") class <- "Metal" else{
      if(Chem.Group[i]!="PAH" & Approved[i]=="Y" & Priority[i] =="" & is.na(Candidate.SVHC[i]) & is.na(SVHC[i])){
        class <- "Current_Use"

      }else{
        class <- "Banned_Listed"
      }
    }
    return(class)
  })
  classes
}


#' get the class of the chemicals
#'
#' @param Chem.Group chemical group
#' @param Approved whether being approved
#' @param Priority "STC", "WFD", "STC, WFD"
#' @param Current_Use Manually identifies
#'
#' @return class being either Metal, Current_Use, or Banned_Listed.
#' @export
#'
#' @examples
getclass <- function(Chem.Group,Approved,Priority,Current_Use){
  classes <- sapply(1:length(Chem.Group),function(i)
  {
    if(Chem.Group[i]=="Metal") class <- "Metal" else{
      if(Current_Use[i]=="YES"){
        class <- "Current_Use"
        if(Chem.Group[i]=="PAH") print("? PAH!")
        if(Priority[i]!="")   print("Priority not empty")
        if(Approved[i]!="Y") {
          print("Approved missing information!")
        }
        # if(Chem.Group[i]!="PAH" & Approved[i]=="Y" & Priority[i] =="" & is.na(Candidate.SVHC[i]) & is.na(SVHC[i])){
        #   class <- "Current_Use"

      }else{
        class <- "Banned_Listed"
      }
    }
    return(class)
  })
  classes
}

#' calculate the number of unique CAS
#'
#' @param Exclusion.CAS exlustion option
#' @param chemclass chemical grouping information data
#'
#' @return number of cas
#' @export
#'
#' @examples
getNcas <- function(Exclusion.CAS,chemclass){
  if(Exclusion.CAS=="Excluded_Metals"){
    cas <- chemclass %>% filter(Chem.Group!="Metal")
  }else if(Exclusion.CAS=="Current_Use"){
    cas <- chemclass %>% filter(class=="Current_Use")
  }else if(Exclusion.CAS=="Metals"){
    cas <- chemclass %>% filter(Chem.Group=="Metal")
  }else if(Exclusion.CAS=="Banned_Listed"){
    #cas1 <- chemclass %>% filter(Chem.Group!="Metal" & Approved=="Y" & Chem.Group!="PAH"& Priority=="")
    #cas <-  chemclass %>% filter(Chem.Group!="Metal")%>% filter(!(CAS %in% cas1$CAS))
    cas <-  chemclass %>% filter(class=="Banned_Listed")
  }else if(Exclusion.CAS=="Excluded_Metals_PAH"){
    cas <-  chemclass %>% filter(Chem.Group!="Metal" & Chem.Group!="PAH")
  }else{
    cas <- chemclass
  }
  ncas <- nrow(cas)
  return(ncas)
}


#' Get the main drivers for each site year combination
#'
#' @param Data.AggBySiteID.2  the WaterBase data
#' @param chemclass  Chemical class data
#' @param LOQ.type  handling of LOQ data
#' @param Exclusion.CAS excluding classes
#' @param Data.SSD.1 SSD information
#' @param threshold take 0.1 as default
#'
#' @return Stacked data frames with LOQ.type and Exclusion.CAS as identifier
#' @export
#'
#' @examples
getDriverData <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),
                          Exclusion.CAS = c("Excluded_None", "Excluded_Metals", "Current_Use", "Metals", "Banned_Listed"),
                          Data.SSD.1,threshold=0.1,include0=F,cutoff=0,exclude.driver=NULL,...){
  Res <- NULL
  for(i in 1:length(LOQ.type)){
    for(j in 1:length(Exclusion.CAS)){
      Res1 <- getDriverData1(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type[i],
                             Exclusion.CAS = Exclusion.CAS[j],
                             Data.SSD.1=Data.SSD.1,threshold = threshold,include0 = include0,cutoff=cutoff,exclude.driver=exclude.driver,...)
      Res1$LOQ.type=LOQ.type[i]
      Res1$Exclusion.CAS = Exclusion.CAS[j]
      Res <- rbind(Res,Res1)
    }
  }
  return(Res)
}


#' helper function to extract driver data
#'
#' @param data
#' @param threshold filter HI by threshold
#' @param include0  whether to include 0 HQ values if the CAS is in the top 3 drivers but there are 0 values in some of the samples.
#'
#' @return
#' @export
#'
#' @examples
driverhelper <- function(data,threshold=0.1,include0=F){
  # Res_nest$data[[1]]
  # # A tibble: 31 x 9
  # procedureLOQVal… resultNumberOfS… resultQualityNu… resultMeanValue resultMaximumVa…
  # <dbl>            <int>            <int>           <dbl>            <dbl>
  #   1            0.015               12               12               0                0
  # 2            0.01                12               12               0                0
  # 3            0.05                12               12               0                0
  # 4            0.001               12               12               0                0
  # 5            0.01                12               12               0                0
  # 6            0.01                12               12               0                0
  # 7            0.01                12               12               0                0
  # 8            0.1                 12               12               0                0
  # 9            0.01                12               12               0                0
  # 10            0.008               12               12               0                0
  # # … with 21 more rows, and 4 more variables: CAS <chr>, N.Measured <int>,
  # #   N.Detected <dbl>, AboveLOQ <lgl>

  data <-data %>% mutate(HQ_Sum_Acute=sum(HQ.Acute,na.rm=T),HQ_Sum_Chronic=sum(HQ.Chronic,na.rm=T))%>%
    mutate(HI_Acute=HQ_Sum_Acute,HI_Chronic=HQ_Sum_Chronic)
  if(data$HI_Acute[1]>threshold){
    data_Acute <- data[order(data$HQ.Acute,decreasing = T),]
    nr <- sum(data_Acute$HQ.Acute>0)
    Acute <- data_Acute[1:nr,] %>% select(c(CAS,HQ.Acute)) %>% mutate(Component=paste("Driver",1:nr))
    if(nr<=3){
      Acute<- rbind(Acute,data.frame(CAS="",HQ.Acute=0,Component="All Other"))
    }else{
      Acute <- Acute[1:3,]
      Acute<- rbind(Acute,data.frame(CAS="",HQ.Acute=sum(data_Acute[4:nr,]$HQ.Acute),Component="All Other"))
    }
    Acute <- Acute %>% rename(HQ=HQ.Acute) %>% mutate(test="Acute",HI=data$HI_Acute[1])
    if(include0){
      added <- data_Acute %>% filter((CAS %in% Acute$CAS) & (HQ.Acute==0)) %>%
        dplyr::select(c(CAS,HQ.Acute))%>% rename(HQ=HQ.Acute)%>% mutate(test="Acute",HI=data$HI_Acute[1])
      added <- left_join(added,Acute[,c("CAS","Component")])
      Acute <- rbind(Acute,added)
    }
  }else{
    # Acute <- data.frame(CAS="",HQ.Acute=data$HI_Acute[1],Component="All Other")
    Acute <- NULL
  }
  if(data$HI_Chronic[1]>threshold){
    data_Chronic <- data[order(data$HQ.Chronic,decreasing = T),]
    nr <- sum(data_Chronic$HQ.Chronic>0)
    Chronic <- data_Chronic[1:nr,] %>% select(c(CAS,HQ.Chronic)) %>% mutate(Component=paste("Driver",1:nr))
    if(nr<=3){
      Chronic<- rbind(Chronic,data.frame(CAS="",HQ.Chronic=0,Component="All Other"))
    }else{
      Chronic <- Chronic[1:3,]
      Chronic<- rbind(Chronic,data.frame(CAS="",HQ.Chronic=sum(data_Chronic[4:nr,]$HQ.Chronic),Component="All Other"))
    }
    Chronic <- Chronic %>% rename(HQ=HQ.Chronic) %>% mutate(test="Chronic",HI=data$HI_Chronic[1])
    if(include0){
      added <- data_Chronic%>% filter((CAS %in% Chronic$CAS) & (HQ.Chronic==0)) %>%
        dplyr::select(c(CAS,HQ.Chronic))%>% rename(HQ=HQ.Chronic) %>% mutate(test="Chronic",HI=data$HI_Chronic[1])
      added <- left_join(added,Chronic[,c("CAS","Component")])
      Chronic <- rbind(Chronic,added)
    }
  }else{
    # Chronic<- data.frame(CAS="",HQ.Chronic=data$HI_Chronic[1],Component="All Other")
    Chronic <- NULL
  }


  return(rbind(Acute,Chronic))
}


#' Get the main drivers for each site year combination
#'
#' @param Data.AggBySiteID.2  the WaterBase data
#' @param chemclass  Chemical class data
#' @param LOQ.type  handling of LOQ data
#' @param Exclusion.CAS excluding classes
#' @param Data.SSD.1 SSD information
#' @param threshold take 0.1 as default
#' @return one data frame with each site year having component 1,2,3,4 with CAS and chemical names.
#' @export
#'
#' @examples
getDriverData1 <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),
                           Exclusion.CAS = c("Excluded_None", "Excluded_Metals",
                                             "Current_Use", "Metals", "Banned_Listed"),
                           Data.SSD.1,
                           threshold=0.1,include0=F,cutoff=0,exclude.driver=NULL,...){
  Data.AggBySiteID.3a <-  preprocessWaterbase(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                                              Data.SSD.1 = Data.SSD.1,Exclusion.CAS=Exclusion.CAS,...)
  sumN <- Data.AggBySiteID.3a %>%  group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>% summarise(Nmeasured=length(CAS))
  Data.AggBySiteID.3a <-left_join(Data.AggBySiteID.3a,sumN) %>%filter(Nmeasured>cutoff)
  if(!is.null(exclude.driver)){
    Data.AggBySiteID.3a <-Data.AggBySiteID.3a %>%filter(!(CAS %in% exclude.driver))
  }
  Res_nest <- Data.AggBySiteID.3a %>%
    group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>% nest()
  Res <- Res_nest %>% mutate(driver=purrr::map(data,.f=driverhelper,threshold=threshold,include0=include0))%>%
    select(-data) %>% unnest(cols = c(driver))

  return(Res)

}

sep_legacy_helper <- function(data,threshold=0.1,include0=F){
  data <-data %>% mutate(HQ_Sum_Acute=sum(HQ.Acute,na.rm=T),HQ_Sum_Chronic=sum(HQ.Chronic,na.rm=T))%>%
    mutate(HI_Acute=HQ_Sum_Acute,HI_Chronic=HQ_Sum_Chronic)
  if(data$HI_Acute[1]>threshold){
    data_Acute <- data[order(data$HQ.Acute,decreasing = T),]
    nr <- sum(data_Acute$HQ.Acute>0) ### Keep all the rows(samples) with HI > 0
    if(nr>0) Acute <- data_Acute[1:nr,] %>% select(c(CAS,HQ.Acute,class)) %>% mutate(Component=class)%>%
      dplyr::group_by(Component)%>%dplyr::summarise(HQ.Acute=sum(HQ.Acute)) else Acute <- NULL

    Acute <- Acute %>% rename(HQ=HQ.Acute) %>% mutate(test="Acute",HI=data$HI_Acute[1])
    if(include0){ ### inlcude HQ=0 CAS samples when the CAS are included in the CASes that appeared in the HI>threshold samples

      ### Note that since for 2 category separation, it does not matter whether include 0 any more, the CAS name
      ### will be lost after joining banned_listed and current_Use.
      ### do nothing!
    }
  }else{
    # Acute <- data.frame(CAS="",HQ.Acute=data$HI_Acute[1],Component="All Other")
    Acute <- NULL
  }
  if(data$HI_Chronic[1]>threshold){
    data_Chronic <- data[order(data$HQ.Chronic,decreasing = T),]
    nr <- sum(data_Chronic$HQ.Chronic>0)
    if(nr>0)  Chronic <- data_Chronic[1:nr,] %>% select(c(CAS,HQ.Chronic,class)) %>% mutate(Component=class)%>%
      dplyr::group_by(Component)%>%dplyr::summarise(HQ.Chronic=sum(HQ.Chronic)) else Chronic <- NULL

    Chronic <- Chronic %>% rename(HQ=HQ.Chronic) %>% mutate(test="Chronic",HI=data$HI_Chronic[1])
    if(include0){
      ### Note that since for 2 category separation, it does not matter whether include 0 any more, the CAS name
      ### will be lost after joining banned_listed and current_Use.
      ### do nothing!
    }
  }else{
    # Chronic<- data.frame(CAS="",HQ.Chronic=data$HI_Chronic[1],Component="All Other")
    Chronic <- NULL
  }


  return(rbind(Acute,Chronic))
}

#' Get Legacy vs Curent Data
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param LOQ.type
#' @param Exclusion.CAS
#' @param Data.SSD.1
#' @param threshold threshold for HI
#' @param include0
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
sep_legacy <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),
                       Exclusion.CAS = c("Excluded_None", "Excluded_Metals",
                                         "Current_Use", "Metals", "Banned_Listed"),
                       Data.SSD.1,
                       threshold=0.1,include0=F,cutoff=0,...){

  Data.AggBySiteID.3a <-  preprocessWaterbase(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                                              Data.SSD.1 = Data.SSD.1,Exclusion.CAS=Exclusion.CAS,...)
  Data.AggBySiteID.3a <- left_join(Data.AggBySiteID.3a,chemclass[,c("CAS","class")])
  if(cutoff>0) {
    sts <- Data.AggBySiteID.3a  %>%  group_by(Site.Year,class) %>% summarise(N.current=length(unique(CAS)))%>%
      filter(class=="Current_Use")%>%filter(N.current>cutoff)
    sts <- sts$Site.Year
    Data.AggBySiteID.3a <- Data.AggBySiteID.3a  %>% filter(Site.Year %in% sts)
  }

  Res_nest <- Data.AggBySiteID.3a %>%  group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>% nest()
  Res <- Res_nest %>% mutate(driver=purrr::map(data,.f=sep_legacy_helper,threshold=threshold,include0=include0))%>% select(-data) %>% unnest(cols = c(driver))

  return(Res)
}


#' Get 6driver vs all other data
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param LOQ.type
#' @param Exclusion.CAS
#' @param Data.SSD.1
#' @param threshold
#' @param include0
#' @param cutoff
#' @param drivers
#'
#' @return
#' @export
#'
#' @examples
sep_6driver <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0"),
                        Exclusion.CAS =  "Excluded_Metals",
                        Data.SSD.1,
                        threshold=0.1,include0=F,cutoff=0,drivers=c("2921882", "15687271", "34123596", "330541", "15307865", "58082"
                        )){

  Data.AggBySiteID.3a <-  preprocessWaterbase(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                                              Data.SSD.1 = Data.SSD.1,Exclusion.CAS=Exclusion.CAS)
  Data.AggBySiteID.3a <- left_join(Data.AggBySiteID.3a,chemclass[,c("CAS","class")])
  if(cutoff>0) {
    sts <- Data.AggBySiteID.3a  %>%  group_by(Site.Year,class) %>% summarise(N.current=length(unique(CAS)))%>%
      filter(class=="Current_Use")%>%filter(N.current>cutoff)
    sts <- sts$Site.Year
    Data.AggBySiteID.3a <- Data.AggBySiteID.3a  %>% filter(Site.Year %in% sts)
  }
  Data.AggBySiteID.3a$class <- "All Other"
  ndriver <- length(drivers)
  Data.AggBySiteID.3a$class[Data.AggBySiteID.3a$CAS %in% drivers] <- paste(ndriver,"Drivers")

  Res_nest <- Data.AggBySiteID.3a %>%  group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>% nest()
  Res <- Res_nest %>% mutate(driver=purrr::map(data,.f=sep_legacy_helper,threshold=threshold,include0=include0))%>% select(-data) %>% unnest(cols = c(driver))
  browser()
  return(Res)
}
#' Generate plot for two component comparison. Legacy vs. Current or 6 drivers and all other. ...
#'
#' @param driverdata3
#' @param nsample3
#' @param sumdriver3
#'
#' @return
#' @export
#'
#' @examples
gen_fig3_2component <- function(driverdata3,nsample3,sumdriver3,plot="boxplot"){
  pbar2 <- ggplot(driverdata3,aes(x=x,y=HQ))+
    geom_bar(aes(fill=Component,col=Component),position=position_stack(reverse = F),stat = "identity")+
    facet_wrap(.~Exclusion.CAS,scale="free")+
    theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
    geom_text(data=nsample3,aes(x=sy,y=hmax,label=paste0("N.Sample = ",sy)),hjust=1)+
    xlab("")+theme(legend.position = "bottom")+geom_hline(yintercept = 1,lty=2)+ylab("HI")




  inset_plot <- NULL
  ## tmp <- expand.grid(LOQ.type = c("LOQ.T0","LOQ.T1"),Exclusion.CAS = c("Excluded_Metals_PAH","Current_Use"))
  k <- length(unique(driverdata3$Exclusion.CAS))
  for(i in 1:1)
    for(j in 1:k){
      inset_plot[[2*(i-1)+j]] <- NA
    }

  names(inset_plot) <- paste0("LOQ.T0","-",nsample3$Exclusion.CAS)

  LOQ <- "LOQ.T0"

  for(Exclusion in nsample3$Exclusion.CAS){
    if(plot=="boxplot") p <- ggplot(data=sumdriver3 %>% filter( Exclusion.CAS==Exclusion),
                                    aes(x=Component,y=Perc,fill=Component))+
        geom_boxplot()  +
        guides(fill=FALSE) +
        theme_bw(base_size=12) +  ## makes everything smaller
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
        theme(panel.background = element_rect(fill="white"),  ## white plot background
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=rel(1),angle = 45,hjust=1,vjust=1), ## tiny axis text
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank())
    if(plot=="violinplot") p <- ggplot(data=sumdriver3 %>% filter( Exclusion.CAS==Exclusion),
                                       aes(x=Component,y=Perc,fill=Component))+
        geom_half_violin() +
        geom_dotplot(binaxis = "y", method="histodot", stackdir="up",dotsize = 0.1,stackratio=0.5)+
        guides(fill=FALSE) +
        theme_bw(base_size=12) +  ## makes everything smaller
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
        theme(panel.background = element_rect(fill="white"),  ## white plot background
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=rel(1),angle = 45,hjust=1,vjust=1), ## tiny axis text
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank())
    inset_plot[[paste0(LOQ,"-",Exclusion)]] <- p
  }




  fig3 <- pbar2+annotation_custom2(grob = ggplotGrob(inset_plot[[1]]), data = data.frame(x=nsample3$sy[1], HQ=nsample3$hmax[1]*0.95,

                                                                                         Exclusion.CAS=nsample3$Exclusion.CAS[1]),
                                   xmin = nsample3$sy[1]*0.15, xmax = nsample3$sy[1]*0.95, ymin = nsample3$hmax[1]*0.15, ymax =nsample3$hmax[1]*0.95)+
    annotation_custom2(grob = ggplotGrob(inset_plot[[2]]), data = data.frame(x=nsample3$sy[2], HQ=nsample3$hmax[2]*0.95,

                                                                             Exclusion.CAS=nsample3$Exclusion.CAS[2]),
                       xmin = nsample3$sy[2]*0.15, xmax = nsample3$sy[2]*0.95, ymin = nsample3$hmax[2]*0.15, ymax =nsample3$hmax[2]*0.95)

  return(fig3)
}

#' Get 1 CAS info w.r.t. LOQ inside the AggBySiteID waterbase data.
#'
#' @param Data.AggBySiteID.2  the WaterBase data
#' @param chemclass  Chemical class data
#' @param Data.SSD.1 SSD information
#'
#' @return data frame with Nmeasured detected, LOQ for each CAS
#' @export
#'
#' @examples
getCasData <- function(Data.AggBySiteID.2,chemclass,Data.SSD.1,Exclusion.CAS="Refined_HQ_Metals_PAH"){

  # casinfo <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS,Site.Year) %>%
  #   dplyr::summarise(AboveLOQ=any(AboveLOQ),LOQ=min(procedureLOQValue))%>% ungroup %>%
  #   dplyr::group_by(CAS)%>%
  #   dplyr::summarise(Nmeasured=length(unique(Site.Year)), Ndetected=sum(AboveLOQ),LOQ=min(LOQ),maxLOQ=max(LOQ))

  casinfoN <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS,Site.Year) %>%
    dplyr::summarise(AboveLOQ=any(AboveLOQ))%>% ungroup %>%
    dplyr::group_by(CAS)%>%
    dplyr::summarise(Nmeasured=length(unique(Site.Year)), Ndetected=sum(AboveLOQ))

  casinfoLOQ <-  Data.AggBySiteID.2 %>% dplyr::group_by(CAS)%>% nest() %>%
    mutate(SummaryLOQ=map(data,~summary(.x$procedureLOQValue)) )%>%
    mutate(tidied = map(SummaryLOQ, tidy)) %>% unnest(tidied) %>%
    dplyr::select(-c(SummaryLOQ,data)) %>%  rename_with(function(x)paste0("LOQ.",x), where(is.numeric))

  casinfo <- left_join(casinfoN,casinfoLOQ)%>% mutate(LOQ=LOQ.minimum)

  casinfoConcMax <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS)%>% nest() %>%
    mutate(SummaryConc=map(data,~summary(.x$resultMaximumValue)) )%>%
    mutate(tidied = map(SummaryConc, tidy)) %>% unnest(tidied) %>%
    dplyr::select(-c(SummaryConc,data)) %>%  rename_with(function(x)paste0("Max.Conc.",x), where(is.numeric))
  casinfo <- left_join(casinfo,casinfoConcMax)

  casinfoConcMean <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS)%>% nest() %>%
    mutate(SummaryConc=map(data,~summary(.x$resultMeanValue)) )%>%
    mutate(tidied = map(SummaryConc, tidy)) %>% unnest(tidied) %>%
    dplyr::select(-c(SummaryConc,data)) %>%  rename_with(function(x)paste0("Mean.Conc.",x), where(is.numeric))
  casinfo <- left_join(casinfo,casinfoConcMax)
  ### Cases refined PAH and Metals.
  if(Exclusion.CAS=="Excluded_None"){
    Data_Excluded_None <- preprocessWaterbase(Data.AggBySiteID.2,chemclass, LOQ.type="LOQ.T0",
                                              Data.SSD.1 = Data.SSD.1,
                                              Exclusion.CAS = "Excluded_None")[,c(1:16,18:20)]%>%
      rename(HQ.Chronic.Mean=HQ.Chronic,HQ.Acute.Mean=HQ.Acute)
    tmp <- preprocessWaterbase(Data.AggBySiteID.2,chemclass, LOQ.type="LOQ.T0",Data.SSD.1 = Data.SSD.1,Exclusion.CAS = "Excluded_None",useAllMax = T)
    Data_Excluded_None <- cbind(Data_Excluded_None,HQ.Chronic.Max=tmp[,20],HQ.Acute.Max=tmp[,19])
    casinfoConcMax <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS)%>% nest() %>%
      mutate(SummaryConc=map(data,~summary(.x$resultMaximumValue)) )%>%
      mutate(tidied = map(SummaryConc, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryConc,data)) %>%  rename_with(function(x)paste0("Max.Conc.",x), where(is.numeric))
    casinfo <- left_join(casinfo,casinfoConcMax)

    casinfoConcMean <- Data.AggBySiteID.2 %>% dplyr::group_by(CAS)%>% nest() %>%
      mutate(SummaryConc=map(data,~summary(.x$resultMeanValue)) )%>%
      mutate(tidied = map(SummaryConc, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryConc,data)) %>%  rename_with(function(x)paste0("Mean.Conc.",x), where(is.numeric))
    casinfo <- left_join(casinfo,casinfoConcMax)



    tmpdat <- Data_Excluded_None %>% dplyr::group_by(CAS)%>% nest()
      hq1<-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Chronic.Mean)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Chronic.Mean.",x), where(is.numeric))
    hq2 <-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Acute.Mean)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Acute.Mean.",x), where(is.numeric))
    hq3<-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Chronic.Max)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Chronic.Max.",x), where(is.numeric))
    hq4 <-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Acute.Max)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Acute.Max.",x), where(is.numeric))
    casinfoHQ <- purrr::reduce(list(hq1,hq2,hq3,hq4), dplyr::left_join)
  }

  if(Exclusion.CAS=="Refined_HQ_Metals_PAH"){
    Data_Refined_HQ_Metals_PAH <- preprocessWaterbase(Data.AggBySiteID.2,chemclass, LOQ.type="LOQ.T0",
                                                   Data.SSD.1 = Data.SSD.1,
                                                   Exclusion.CAS = "Refined_HQ_Metals_PAH")[,c(1:16,18:20)]%>%
      rename(HQ.Chronic.Mean=HQ.Chronic,HQ.Acute.Mean=HQ.Acute)
    tmp <- preprocessWaterbase(Data.AggBySiteID.2,chemclass, LOQ.type="LOQ.T0",Data.SSD.1 = Data.SSD.1,
                               Exclusion.CAS = "Refined_HQ_Metals_PAH",useAllMax = T)
    Data_Refined_HQ_Metals_PAH <- cbind(Data_Refined_HQ_Metals_PAH,HQ.Chronic.Max=tmp[,20],HQ.Acute.Max=tmp[,19])



    casinfoConcMax <- Data_Refined_HQ_Metals_PAH %>% dplyr::group_by(CAS)%>% nest() %>%
      mutate(SummaryConc=map(data,~summary(.x$resultMaximumValue)) )%>%
      mutate(tidied = map(SummaryConc, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryConc,data)) %>%  rename_with(function(x)paste0("Max.Conc.",x), where(is.numeric))
    casinfo <- left_join(casinfo,casinfoConcMax)

    casinfoConcMean <- Data_Refined_HQ_Metals_PAH %>% dplyr::group_by(CAS)%>% nest() %>%
      mutate(SummaryConc=map(data,~summary(.x$resultMeanValue)) )%>%
      mutate(tidied = map(SummaryConc, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryConc,data)) %>%  rename_with(function(x)paste0("Mean.Conc.",x), where(is.numeric))
    casinfo <- left_join(casinfo,casinfoConcMax)



    tmpdat <- Data_Refined_HQ_Metals_PAH %>% dplyr::group_by(CAS)%>% nest()
    hq1<-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Chronic.Mean)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Chronic.Mean.",x), where(is.numeric))
    hq2 <-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Acute.Mean)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Acute.Mean.",x), where(is.numeric))
    hq3<-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Chronic.Max)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Chronic.Max.",x), where(is.numeric))
    hq4 <-tmpdat %>%
      mutate(SummaryHQ=map(data,~summary(.x$HQ.Acute.Max)) )%>%
      mutate(tidied = map(SummaryHQ, tidy)) %>% unnest(tidied) %>%
      dplyr::select(-c(SummaryHQ,data)) %>%  rename_with(function(x)paste0("HQ.Acute.Max.",x), where(is.numeric))
    casinfoHQ <- purrr::reduce(list(hq1,hq2,hq3,hq4), dplyr::left_join)
  }
  casinfo <- left_join(casinfo,casinfoHQ)
  chemclass$CAS <- as.character(chemclass$CAS)
  casdata <- left_join(casinfo,chemclass)
  ssdInfo <- Data.SSD.1[,c("CAS.","X10LogSSDMedianConcentration.ug.L..MuAcute.EC50","X10LogSSDMedianConcentration.ug.L..MuChronic.NOEC", "X10LogSSDSlope.ug.L..SigmaAcute.EC50","X10LogSSDSlope.ug.L..SigmaChronic.NOEC")] %>% mutate(CAS.=as.character(CAS.))
  names(ssdInfo) <- c("CAS","SSDLOG10.Mu.Acute.EC50", "SSDLOG10.Mu.chronic.NOEC", "SSDLOG10.Sigma.Acute.EC50","SSDLOG10.Sigma.chronic.NOEC")

  casdata <- left_join(casdata,ssdInfo) %>%
    mutate(LOQ.Chronic=LOQ/10^qnorm(0.05,SSDLOG10.Mu.chronic.NOEC, SSDLOG10.Sigma.chronic.NOEC),
           LOQ.Acute=LOQ/10^qnorm(0.5,SSDLOG10.Mu.Acute.EC50, SSDLOG10.Sigma.Acute.EC50)
    )%>%
    mutate(HC50.Acute = 10^qnorm(0.5,SSDLOG10.Mu.Acute.EC50, SSDLOG10.Sigma.Acute.EC50),
           HC05.Chronic = 10^qnorm(0.05,SSDLOG10.Mu.chronic.NOEC,SSDLOG10.Sigma.chronic.NOEC))
  return(casdata)
}




#' Get HQ chronic data for driver data
#'
#' @param Data.AggBySiteID.2  the WaterBase data
#' @param chemclass  Chemical class data
#' @param Data.SSD.1 SSD information
#' @param driverdata driverdata
#' @return data frame with Nmeasured detected, LOQ for each CAS
#' @export
#'
#' @examples
getHQData_Chronic <- function(Data.AggBySiteID.2,chemclass,Data.SSD.1,LOQ.type,driverdata,Exclusion.CAS="Excluded_Metals"){
  excludecas <- Exclusion.CAS
  driverdata <- driverdata %>% filter(test=="Chronic" & CAS!="" & Exclusion.CAS==excludecas & Component!="All Other") ## note that filter could postentially not working
  unique(driverdata1$Exclusion.CAS)
  HQ.data <- NULL
  for(i in 1:length(LOQ.type)){
    driverdata1 <- driverdata %>% filter(LOQ.type==LOQ.type[i])
    # Data.AggBySiteID.3a <-  selectdata(Data.AggBySiteID.2,chemclass,LOQ.type =LOQ.type[i] ,Exclusion.CAS = "Excluded_Metals")
    Data.AggBySiteID.3a <-  preprocessWaterbase(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type[i],
                                                Data.SSD.1 = Data.SSD.1,Exclusion.CAS= Exclusion.CAS,
                                                HU=F,PAF=F)

    Data.AggBySiteID.3a<- Data.AggBySiteID.3a %>%  group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>%
      mutate(HI_Acute=sum(HQ.Acute,na.rm=T),HI_Chronic=sum(HQ.Chronic,na.rm=T))
    Data.AggBySiteID.3a <- Data.AggBySiteID.3a %>% dplyr::select(c(Site.Year,HQ.Acute,HQ.Chronic,CAS,HI_Chronic,HI_Acute))%>%
      filter(CAS %in% driverdata1$CAS) %>% mutate(LOQ.type=LOQ.type[i])
    HQ.data <- rbind(HQ.data,Data.AggBySiteID.3a)
  }

  return(HQ.data)
}


#' Get the relevant data w.r.t. Site and Year
#'
#' @param Data.AggBySiteID.2  the WaterBase data
#' @param chemclass  Chemical class data
#' @param LOQ.type  handling of LOQ data
#' @param Exclusion.CAS excluding classes
#' @param Data.SSD.1 SSD information
#'
#' @return Stacked data frames with LOQ.type and Exclusion.CAS as identifier
#' @export
#'
#' @examples
getResData <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),
                       Exclusion.CAS = c("Excluded_None", "Excluded_Metals", "Current_Use", "Metals", "Banned_Listed"),
                       exclude.substances=NULL,
                       Data.SSD.1,...){
  Res <- NULL
  for(i in 1:length(LOQ.type)){
    for(j in 1:length(Exclusion.CAS)){
      Res1 <- getResData1(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type[i],
                          Exclusion.CAS = Exclusion.CAS[j],exclude.substances=exclude.substances,
                          Data.SSD.1=Data.SSD.1,...)
      Res1$LOQ.type=LOQ.type[i]
      Res1$Exclusion.CAS = Exclusion.CAS[j]
      Res <- rbind(Res,Res1)
    }
  }
  Res$Exclusion.CAS <- factor(Res$Exclusion.CAS,levels=Exclusion.CAS)
  return(Res)
}

#' Get 1 dataset
#'
#' @param Data.AggBySiteID.2  the WaterBase data
#' @param chemclass  Chemical class data
#' @param LOQ.type  handling of LOQ data
#' @param Exclusion.CAS excluding classes
#' @param Data.SSD.1 SSD information
#' @param exclude.substances exclude substances CAS
#'
#' @return data frame with MCR, HQ, HI, PAF, and other additional informaion
#' @export
#'
#' @examples
getResData1 <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0", "LOQ.T1"),
                        Exclusion.CAS = c("Excluded_None", "Excluded_Metals", "Current_Use", "Metals", "Banned_Listed","Refined_HQ_Metals_PAH"),
                        exclude.substances=NULL,
                        Data.SSD.1,...){
  Data.AggBySiteID.3a <-  preprocessWaterbase(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                                              Data.SSD.1 = Data.SSD.1,Exclusion.CAS=Exclusion.CAS,...)
  if(!is.null(exclude.substances)){
    Data.AggBySiteID.3a <-  Data.AggBySiteID.3a %>% filter(!(CAS %in% exclude.substances))
  }
  Res <- Data.AggBySiteID.3a %>%  group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>%
    summarise(HU_Sum_Acute=sum(HU.Acute.EC50,na.rm=T),HU_Sum_Chronic=sum(HU.Chronic.NOEC,na.rm=T),
              HQ_Sum_Acute=sum(HQ.Acute,na.rm=T),HQ_Sum_Chronic=sum(HQ.Chronic,na.rm=T),
              HU_Max_Acute=max(HU.Acute.EC50,na.rm=T),HU_Max_Chronic=max(HU.Chronic.NOEC,na.rm=T),
              HQ_Max_Acute=max(HQ.Acute,na.rm=T),HQ_Max_Chronic=max(HQ.Chronic,na.rm=T),
              PAF_Max_Acute=max(PAF.Acute.EC50,na.rm=T),PAF_Max_Chronic=max(PAF.Chronic.NOEC,na.rm=T),
              Nmeasured=length(CAS),Ndetected=sum(resultMaximumValue>0),Ratio.MD=Ndetected/Nmeasured) %>%
    mutate(HI_Acute=HQ_Sum_Acute,HI_Chronic=HQ_Sum_Chronic,
           MaxHQ.Acute=HQ_Max_Acute,MaxHQ.Chronic=HQ_Max_Chronic)%>%
    mutate(msPAF_Sum_Acute=pnorm(log10(HU_Sum_Acute),0,0.7),msPAF_Max_Acute=pnorm(log10(HU_Max_Acute),0,0.7),
           msPAF_Sum_Chronic=pnorm(log10(HU_Sum_Chronic),0,0.7),msPAF_Max_Chronic=pnorm(log10(HU_Max_Chronic),0,0.7))%>%
    mutate(MCR.HC05.Chronic.NOEC = HQ_Sum_Chronic/HQ_Max_Chronic,
           MCR.HC50.Acute.EC50=HQ_Sum_Acute/HQ_Max_Acute)
  Res <- Res %>% ungroup %>% group_by(monitoringSiteIdentifier) %>% mutate(maxHI_site=max(HQ_Sum_Chronic,na.rm=T))  ## Aggregate over "Year"
  Res <- Res %>% ungroup %>% mutate(Group_AF1_Chronic=get_group(HI=HQ_Sum_Chronic,HQ=HQ_Max_Chronic,MCR=MCR.HC05.Chronic.NOEC,AF=1),
                                    #Group_AF3_Chronic=get_group(HI=HQ_Sum_Chronic,HQ=HQ_Max_Chronic,MCR=MCR.HC05.Chronic.NOEC,AF=3),
                                    #Group_AF10_Chronic=get_group(HI=HQ_Sum_Chronic,HQ=HQ_Max_Chronic,MCR=MCR.HC05.Chronic.NOEC,AF=10),
                                    Group_AF1_Acute=get_group(HI=HQ_Sum_Acute,HQ=HQ_Max_Acute,MCR=MCR.HC50.Acute.EC50,AF=1))#,
  #Group_AF3_Acute=get_group(HI=HQ_Sum_Acute,HQ=HQ_Max_Acute,MCR=MCR.HC50.Acute.EC50,AF=3),
  # Group_AF10_Acute=get_group(HI=HQ_Sum_Acute,HQ=HQ_Max_Acute,MCR=MCR.HC50.Acute.EC50,AF=10))
  return(Res)

}



#' Title
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param Exclusion.CAS
#' @param exclude.substances
#' @param Data.SSD.1
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
getResData_LOQ <- function(Data.AggBySiteID.2,chemclass,
                        Exclusion.CAS = c("Excluded_Metals"),
                        exclude.substances=NULL,
                        Data.SSD.1,LOQ.type="LOQ.T1",...){
  Data.AggBySiteID.3a <-  Data.AggBySiteID.2 %>% mutate(resultMaximumValue=procedureLOQValue,
                                                        resultMeanValue=procedureLOQValue,AboveLOQ=TRUE)
  Data.AggBySiteID.3a <-  preprocessWaterbase(Data.AggBySiteID.3a,chemclass,LOQ.type = LOQ.type,
                                              Data.SSD.1 = Data.SSD.1,Exclusion.CAS=Exclusion.CAS,...)
  if(!is.null(exclude.substances)){
    Data.AggBySiteID.3a <-  Data.AggBySiteID.3a %>% filter(!(CAS %in% exclude.substances))
  }
  Res <- Data.AggBySiteID.3a %>%  group_by(monitoringSiteIdentifier,phenomenonTimeReferenceYear,Site.Year) %>%
    summarise(HU_Sum_Acute=sum(HU.Acute.EC50,na.rm=T),HU_Sum_Chronic=sum(HU.Chronic.NOEC,na.rm=T),
              HQ_Sum_Acute=sum(HQ.Acute,na.rm=T),HQ_Sum_Chronic=sum(HQ.Chronic,na.rm=T),
              HU_Max_Acute=max(HU.Acute.EC50,na.rm=T),HU_Max_Chronic=max(HU.Chronic.NOEC,na.rm=T),
              HQ_Max_Acute=max(HQ.Acute,na.rm=T),HQ_Max_Chronic=max(HQ.Chronic,na.rm=T),
              PAF_Max_Acute=max(PAF.Acute.EC50,na.rm=T),PAF_Max_Chronic=max(PAF.Chronic.NOEC,na.rm=T),
              Nmeasured=length(CAS),Ndetected=sum(resultMaximumValue>0),Ratio.MD=Ndetected/Nmeasured) %>%
    mutate(HI_Acute=HQ_Sum_Acute,HI_Chronic=HQ_Sum_Chronic,
           MaxHQ.Acute=HQ_Max_Acute,MaxHQ.Chronic=HQ_Max_Chronic)%>%
    mutate(msPAF_Sum_Acute=pnorm(log10(HU_Sum_Acute),0,0.7),msPAF_Max_Acute=pnorm(log10(HU_Max_Acute),0,0.7),
           msPAF_Sum_Chronic=pnorm(log10(HU_Sum_Chronic),0,0.7),msPAF_Max_Chronic=pnorm(log10(HU_Max_Chronic),0,0.7))%>%
    mutate(MCR.HC05.Chronic.NOEC = HQ_Sum_Chronic/HQ_Max_Chronic,
           MCR.HC50.Acute.EC50=HQ_Sum_Acute/HQ_Max_Acute)
  Res <- Res %>% ungroup %>% group_by(monitoringSiteIdentifier) %>% mutate(maxHI_site=max(HQ_Sum_Chronic,na.rm=T))  ## Aggregate over "Year"
  Res <- Res %>% ungroup %>% mutate(Group_AF1_Chronic=get_group(HI=HQ_Sum_Chronic,HQ=HQ_Max_Chronic,MCR=MCR.HC05.Chronic.NOEC,AF=1),
                                    #Group_AF3_Chronic=get_group(HI=HQ_Sum_Chronic,HQ=HQ_Max_Chronic,MCR=MCR.HC05.Chronic.NOEC,AF=3),
                                    #Group_AF10_Chronic=get_group(HI=HQ_Sum_Chronic,HQ=HQ_Max_Chronic,MCR=MCR.HC05.Chronic.NOEC,AF=10),
                                    Group_AF1_Acute=get_group(HI=HQ_Sum_Acute,HQ=HQ_Max_Acute,MCR=MCR.HC50.Acute.EC50,AF=1))#,
  #Group_AF3_Acute=get_group(HI=HQ_Sum_Acute,HQ=HQ_Max_Acute,MCR=MCR.HC50.Acute.EC50,AF=3),
  # Group_AF10_Acute=get_group(HI=HQ_Sum_Acute,HQ=HQ_Max_Acute,MCR=MCR.HC50.Acute.EC50,AF=10))
  Res$Exclusion.CAS <- Exclusion.CAS
  Res$LOQ.type <- LOQ.type
  return(Res)

}



get_group_1 <- function(HI,HQ_max,MCR,AF=c(1,3,10)){
  cutoff <- 1/AF

  if(HI>cutoff & HQ_max > cutoff ) {
    gr <- "I"
  }else{
    if(HI<=cutoff & HQ_max <= cutoff) gr <-"II" else{
      if(HI > cutoff & HQ_max <= cutoff){
        if(MCR <= 2) gr <- "IIIA" else gr <- "IIIB"
      }else{
        gr <- NA
      }
    }
  }
  return(gr)
}

#'  CEFIC-MIAT decision tree table (Valloton &Proce 2016)
#'
#' @param HIv vector of HI
#' @param HQ_maxv vector of max HQ
#' @param MCRv vector of MCR
#' @param AF assessment factor default being 1
#'
#' @return result of classification of I,II, IIIa, IIIb
#' @export
#'
#' @examples
get_group <- function(HIv,HQ_maxv,MCRv,AF=1){
  sapply(1:length(HIv), function(x)get_group_1(HIv[x],HQ_maxv[x],MCRv[x],AF))
}


#' Function to generate one row of figure 1
#'
#' @param Res summarized data based on waterbase db
#' @param cutoff threshold for number of measured chemicals
#' @param df
#' @param test
#' @param used
#' @param exclude.cas
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
generate_fig1_1 <- function(Res,
                            cutoff=0,df = data.frame(x=1e-6,y=8),
                            test="Chronic",used=c("Mean"),
                            exclude.cas=NULL,chemclass,label3B=TRUE,...){
  ## c("Diclofenac", "Ibuprofen", "Chlorpyrifos", "diuron", "Isoproturon", "Caffeine")
  # if(used=="Max") Res <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),

  #                   Exclusion.CAS = c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),Data.SSD.1,useAllMax=T)
  # if(used=="Mean") Res <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),
  #                   Exclusion.CAS = c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),Data.SSD.1,useAllMean=T)

  ##

  ## define the recting for coloring different areas in the plot
  recting <- rbind(data.frame(xmin=-Inf,xmax=1,ymin=-Inf,ymax=Inf,col="green"),
                   data.frame(xmin=1,xmax=Inf,ymin=0,ymax=2,col="red"),
                   data.frame(xmin=1,xmax=Inf,ymin=2,ymax=Inf,col="yellow"))
  Res <- Res %>% filter(LOQ.type=="LOQ.T0" & Exclusion.CAS != "Excluded_None") %>% droplevels(.)

  if(cutoff>0){
    ## Keep only stations that has N.measured>cuttoff current use chemicals.
    sts <-  Res%>%filter(Exclusion.CAS =="Current_Use" & Nmeasured>cutoff)
    sts <- unique(sts$Site.Year)
    Res <- Res%>%filter(Site.Year %in% sts)
  }


  ## Get number of Cas and number of sites.
  N_CAS_Site <- Res %>% group_by(LOQ.type,Exclusion.CAS) %>%
    summarise(N.Site=length(unique(monitoringSiteIdentifier)),
              N.Sample=length(monitoringSiteIdentifier)) %>%
    mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))


  if(!is.null(exclude.cas)){
    if(!is.na(as.numeric(exclude.cas[1]))){
      chemclass1 <- chemclass %>% filter( !(CAS %in% exclude.cas))

    }else chemclass1 <- chemclass %>% filter( !(Substance %in% exclude.cas))
  }else{chemclass1 <- chemclass}
  N.CAS <- sapply(c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),function(x)getNcas(x,chemclass1))
  ## N_CAS_Site$N.CAS <- rep(N.CAS[unique(N_CAS_Site$Exclusion.CAS)],2) ## When there multiple LOQ types.
  N_CAS_Site$N.CAS <-N.CAS[unique(N_CAS_Site$Exclusion.CAS)]
  N_CAS_Site <- data.frame(df,N_CAS_Site)
  ## use only chronic.
  # N_CAS_Site$y <- 8

  if(test=="Acute"){
    Nclass <- Res %>% group_by(LOQ.type,Exclusion.CAS,Group_AF1_Acute) %>%
      summarise(N=length(monitoringSiteIdentifier))%>%
      left_join(N_CAS_Site[,c("LOQ.type","Exclusion.CAS","N.Sample")])%>%
      mutate(Perc=N/N.Sample*100)

    if(label3B){
      Nclass <- Nclass[,c("LOQ.type","Exclusion.CAS","Perc","Group_AF1_Acute")]%>%
        rename(Group=Group_AF1_Acute)%>%filter(Group!="IIIA")
    }else{
      Nclass1 <- Nclass %>% filter(Group_AF1_Acute=="IIIA" | Group_AF1_Acute =="IIIB")%>%
        group_by(LOQ.type,Exclusion.CAS) %>% summarise(Perc=sum(Perc)) %>%
        mutate(Group="III")

      Nclass <- rbind(Nclass1,Nclass[,c("LOQ.type","Exclusion.CAS","Perc","Group_AF1_Acute")]%>%
                        rename(Group=Group_AF1_Acute)%>%filter(Group=="I"|Group=="II"))
    }
    df1 <- data.frame(x=1e-6,y=7,Label="No concern",Group="II")
    df1 <- rbind(df1,data.frame(x=20,y=5,Label="Single Substance Risk",Group="I"))
    if(label3B){
      df1 <- rbind(df1,data.frame(x=5,y=10,Label="Multiple Chemicals Risk",Group="IIIB"))
    }else df1 <- rbind(df1,data.frame(x=5,y=10,Label="Multiple Chemicals Risk",Group="III"))
    Nclass <- Nclass %>% left_join(df1)
    Res_1 <- Res %>% mutate(MCR=MCR.HC50.Acute.EC50,HI=HI_Acute,Group_AF1=Group_AF1_Acute)
  }

  if(test=="Chronic"){

    Nclass <- Res %>% group_by(LOQ.type,Exclusion.CAS,Group_AF1_Chronic) %>% summarise(N=length(monitoringSiteIdentifier))%>% left_join(N_CAS_Site[,c("LOQ.type","Exclusion.CAS","N.Sample")])%>% mutate(Perc=N/N.Sample*100)

    if(label3B){
      Nclass <- Nclass[,c("LOQ.type","Exclusion.CAS","Perc","Group_AF1_Chronic")]%>%
        rename(Group=Group_AF1_Chronic)%>%filter(Group!="IIIA")
    }else{
      Nclass1 <- Nclass %>% filter(Group_AF1_Chronic=="IIIA" | Group_AF1_Chronic =="IIIB")%>% group_by(LOQ.type,Exclusion.CAS) %>% summarise(Perc=sum(Perc)) %>% mutate(Group="III")

      Nclass <- rbind(Nclass1,Nclass[,c("LOQ.type","Exclusion.CAS","Perc","Group_AF1_Chronic")]%>%rename(Group=Group_AF1_Chronic)%>%filter(Group=="I"|Group=="II"))
    }
    df1 <- data.frame(x=1e-6,y=4,Label="No concern",Group="II")
    df1 <- rbind(df1,data.frame(x=50,y=5,Label="Single Substance Risk",Group="I"))
    if(label3B){
      df1 <- rbind(df1,data.frame(x=5,y=9,Label="Multiple Chemicals Risk",Group="IIIB"))
    }else df1 <- rbind(df1,data.frame(x=5,y=9,Label="Multiple Chemicals Risk",Group="III"))
    Nclass <- Nclass %>% left_join(df1)
    Res_1 <- Res %>% mutate(MCR=MCR.HC05.Chronic.NOEC,HI=HI_Chronic,Group_AF1=Group_AF1_Chronic)
  }

  N_CAS_Site_1 <- N_CAS_Site%>% filter(LOQ.type=="LOQ.T0" & Exclusion.CAS != "Excluded_None")%>% droplevels(.)
  Nclass_1 <- Nclass%>% filter(LOQ.type=="LOQ.T0" & Exclusion.CAS != "Excluded_None")%>% droplevels(.)
  Res_1 <- Res_1 %>%  mutate(test=test,cutoff=cutoff,used=used)
  ##

  p1 <- ggplot(Res_1,aes(x=HI,y=MCR))+
    geom_point(aes(col=Group_AF1))+geom_text(data=N_CAS_Site_1,aes(x=x,y=y,label=paste0("N.Site = ",N.Site,"\nN.CAS = ",N.CAS,"\nN.Sample = ",N.Sample)))+
    geom_text(data=Nclass_1,aes(x=x,y=y,label=paste0(Group,"\n",formatC(Perc,digits = 1,format="f"),"%")),vjust="inward",hjust=0.6)  +
    scale_x_log10()+geom_vline(xintercept = 1)+facet_grid(LOQ.type~Exclusion.CAS,scale="free")+
    annotate("rect", xmin = 1e-10, xmax = 1, ymin = recting[1,3], ymax = recting[1,4],fill="yellow",alpha=0.2)+
    geom_line(data=data.frame(x=seq(1,10,length=100),y=seq(1,10,length=100)),aes(x=x,y=y))+
    scale_y_continuous(breaks=c(2,4,6,8),limits=(c(1,10)),expand = c(0,0))+theme(legend.position = "bottom")
  ##

  return(list(p1=p1,Res_1=Res_1,Nclass_1=Nclass_1,N_CAS_Site_1=N_CAS_Site_1))
}




#' Generate the CEFIC-MIAT plot only for LOQ data.
#'
#' @param Res_LOQ
#' @param cutoff
#' @param df
#' @param test
#' @param used
#' @param exclude.cas
#' @param chemclass
#' @param label3B
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
generate_fig1_LOQ <- function(Res_LOQ,
                            cutoff=5,df = data.frame(x=1e-6,y=8),
                            test="Chronic",used=c("Mean"),
                            exclude.cas=NULL,chemclass,label3B=TRUE,...){
  ## c("Diclofenac", "Ibuprofen", "Chlorpyrifos", "diuron", "Isoproturon", "Caffeine")
  # if(used=="Max") Res <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),

  #                   Exclusion.CAS = c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),Data.SSD.1,useAllMax=T)
  # if(used=="Mean") Res <- getResData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0","LOQ.T1"),
  #                   Exclusion.CAS = c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),Data.SSD.1,useAllMean=T)

  ##
  Res <- Res_LOQ
  ## define the recting for coloring different areas in the plot
  recting <- rbind(data.frame(xmin=-Inf,xmax=1,ymin=-Inf,ymax=Inf,col="green"),
                   data.frame(xmin=1,xmax=Inf,ymin=0,ymax=2,col="red"),
                   data.frame(xmin=1,xmax=Inf,ymin=2,ymax=Inf,col="yellow"))




  ## Get number of Cas and number of sites.
  N_CAS_Site <- Res %>% group_by(LOQ.type,Exclusion.CAS) %>%
    summarise(N.Site=length(unique(monitoringSiteIdentifier)),
              N.Sample=length(monitoringSiteIdentifier)) %>%
    mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))


  if(!is.null(exclude.cas)){
    if(!is.na(as.numeric(exclude.cas[1]))){
      chemclass1 <- chemclass %>% filter( !(CAS %in% exclude.cas))

    }else chemclass1 <- chemclass %>% filter( !(Substance %in% exclude.cas))
  }else{chemclass1 <- chemclass}
  N.CAS <- sapply(c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),function(x)getNcas(x,chemclass1))
  ## N_CAS_Site$N.CAS <- rep(N.CAS[unique(N_CAS_Site$Exclusion.CAS)],2) ## When there multiple LOQ types.
  N_CAS_Site$N.CAS <-N.CAS[unique(N_CAS_Site$Exclusion.CAS)]
  N_CAS_Site <- data.frame(df,N_CAS_Site)
  ## use only chronic.
  # N_CAS_Site$y <- 8



  if(test=="Chronic"){

    Nclass <- Res %>% group_by(LOQ.type,Exclusion.CAS,Group_AF1_Chronic) %>% summarise(N=length(monitoringSiteIdentifier))%>% left_join(N_CAS_Site[,c("LOQ.type","Exclusion.CAS","N.Sample")])%>% mutate(Perc=N/N.Sample*100)

    if(label3B){
      Nclass <- Nclass[,c("LOQ.type","Exclusion.CAS","Perc","Group_AF1_Chronic")]%>%
        rename(Group=Group_AF1_Chronic)%>%filter(Group!="IIIA")
    }else{
      Nclass1 <- Nclass %>% filter(Group_AF1_Chronic=="IIIA" | Group_AF1_Chronic =="IIIB")%>% group_by(LOQ.type,Exclusion.CAS) %>% summarise(Perc=sum(Perc)) %>% mutate(Group="III")

      Nclass <- rbind(Nclass1,Nclass[,c("LOQ.type","Exclusion.CAS","Perc","Group_AF1_Chronic")]%>%rename(Group=Group_AF1_Chronic)%>%filter(Group=="I"|Group=="II"))
    }
    df1 <- data.frame(x=1e-6,y=4,Label="No concern",Group="II")
    df1 <- rbind(df1,data.frame(x=50,y=5,Label="Single Substance Risk",Group="I"))
    if(label3B){
      df1 <- rbind(df1,data.frame(x=5,y=9,Label="Multiple Chemicals Risk",Group="IIIB"))
    }else df1 <- rbind(df1,data.frame(x=5,y=9,Label="Multiple Chemicals Risk",Group="III"))
    Nclass <- Nclass %>% left_join(df1)
    Res_1 <- Res %>% mutate(MCR=MCR.HC05.Chronic.NOEC,HI=HI_Chronic,Group_AF1=Group_AF1_Chronic)
  }

  N_CAS_Site_1 <- N_CAS_Site%>% filter(LOQ.type=="LOQ.T1" & Exclusion.CAS != "Excluded_None")%>% droplevels(.)
  Nclass_1 <- Nclass%>% filter(LOQ.type=="LOQ.T1" & Exclusion.CAS != "Excluded_None")%>% droplevels(.)
  Res_1 <- Res_1 %>%  mutate(test=test,cutoff=cutoff,used=used)
  ##

  p1 <- ggplot(Res_1,aes(x=HI,y=MCR))+
    geom_point(aes(col=Group_AF1))+geom_text(data=N_CAS_Site_1,aes(x=x,y=y,label=paste0("N.Site = ",N.Site,"\nN.CAS = ",N.CAS,"\nN.Sample = ",N.Sample)))+
    geom_text(data=Nclass_1,aes(x=x,y=y,label=paste0(Group,"\n",formatC(Perc,digits = 1,format="f"),"%")),vjust="inward",hjust=0.6)  +
    scale_x_log10()+geom_vline(xintercept = 1)+facet_grid(LOQ.type~Exclusion.CAS,scale="free")+
    annotate("rect", xmin = 1e-10, xmax = 1, ymin = recting[1,3], ymax = recting[1,4],fill="yellow",alpha=0.2)+
    geom_line(data=data.frame(x=seq(1,10,length=100),y=seq(1,10,length=100)),aes(x=x,y=y))+
    scale_y_continuous(limits=(c(1,10)),expand = c(0,0))+theme(legend.position = "bottom")
  ##

  return(list(p1=p1,Res_1=Res_1,Nclass_1=Nclass_1,N_CAS_Site_1=N_CAS_Site_1))
}




#' layerd plot definitions in inset_plot
#'
#' @param grob
#' @param xmin
#' @param xmax
#' @param ymin
#' @param ymax
#' @param data
#'
#' @return
#' @export
#'
#' @examples
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data)
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob,
                                          xmin = xmin, xmax = xmax,
                                          ymin = ymin, ymax = ymax))
}




#' Generate figure 2 (Figure 3 in Manuscript)
#'
#' @param Res processed waterbase data
#' @param cutoff 0, not excluding sites
#' @param df indication of plotting legends position in the figure
#' @param test "Chronic" or "Acute"
#' @param used "Mean" or "Max"
#' @param exclude.cas whether to exclude some chemicals
#' @param chemclass chemical information
#' @param ... other to be passed arguments
#'
#' @return
#' @export
#'
#' @examples
generate_fig2 <- function(Res,
                          cutoff=0,df = data.frame(x=1e-6,y=8),
                          test="Chronic",used=c("Mean"),
                          exclude.cas=NULL,chemclass,...){

  Res <- Res %>% filter(Exclusion.CAS != "Excluded_None") %>% droplevels(.)

  if(cutoff>0){
    ## Keep only stations that has N.measured>cuttoff current use chemicals.
    sts <-  Res%>%filter(Exclusion.CAS =="Current_Use" & Nmeasured>cutoff)
    sts <- unique(sts$Site.Year)
    Res <- Res%>%filter(Site.Year %in% sts)
  }


  ## Get number of Cas and number of sites.
  N_CAS_Site <- Res %>% group_by(LOQ.type,Exclusion.CAS) %>%
    summarise(N.Site=length(unique(monitoringSiteIdentifier)),N.Sample=length(monitoringSiteIdentifier)) %>%
    mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))


  if(!is.null(exclude.cas)){
    chemclass1 <- chemclass %>% filter( !(Substance %in% exclude.cas))
  }else{chemclass1 <- chemclass}
  N.CAS <- sapply(c("Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use"),function(x)getNcas(x,chemclass1))
  ## N_CAS_Site$N.CAS <- rep(N.CAS[unique(N_CAS_Site$Exclusion.CAS)],2) ## When there multiple LOQ types.
  N_CAS_Site$N.CAS <-N.CAS[as.character((N_CAS_Site$Exclusion.CAS))]#rep(N.CAS[unique(N_CAS_Site$Exclusion.CAS)],length(unique(N_CAS_Site$LOQ.type)))
  N_CAS_Site <- data.frame(df,N_CAS_Site)
  ## use only chronic.
  # N_CAS_Site$y <- 8

  if(test=="Acute"){

    Res_1 <- Res %>% mutate(MCR=MCR.HC50.Acute.EC50,HI=HI_Acute,Group_AF1=Group_AF1_Acute)
  }

  if(test=="Chronic"){


    Res_1 <- Res %>% mutate(MCR=MCR.HC05.Chronic.NOEC,HI=HI_Chronic,Group_AF1=Group_AF1_Chronic)
  }
  ## df <- data.frame(x=200,y=1e-6)


  N_CAS_Site_1 <- N_CAS_Site %>% droplevels(.)

  Res_1 <- Res_1 %>%  mutate(test=test,cutoff=cutoff,used=used)
  ######

  tmp <- Res_1[,c("HI","HI_Acute","HI_Chronic","Ndetected","Nmeasured","LOQ.type","Exclusion.CAS")]



  tmp <- tmp %>% pivot_longer(cols = c(Ndetected,Nmeasured),values_to="N",names_to="Ntype") # %>%
  tmp$Exclusion.CAS <- factor(tmp$Exclusion.CAS,levels=(c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))
  tmp$Ntype<- factor(tmp$Ntype,levels=(c("Nmeasured","Ndetected")))
  sumtmp <- tmp %>% group_by(LOQ.type,Exclusion.CAS,Ntype)%>%summarise(Perc=sum(HI<1)/length(HI)*100,x1=max(N))
  sumtmp <- cbind(df,sumtmp)

  N_CAS_Site_1$x <- sumtmp$x[1]
  N_CAS_Site_1$y <- 1e-03
  N_CAS_Site_1$Ntype <- "Ndetected"
  N_CAS_Site_1 <- rbind(N_CAS_Site_1,N_CAS_Site_1%>% mutate(Ntype="Nmeasured"))

  ##
  N_CAS_Site_1$Ntype<- factor(N_CAS_Site_1$Ntype,levels=(c("Nmeasured","Ndetected")))
  N_CAS_Site_1$Exclusion.CAS <- factor(N_CAS_Site_1$Exclusion.CAS,levels=(c("Excluded_None","Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))
  N_CAS_Site_1 <- N_CAS_Site_1%>% filter(Exclusion.CAS=="Excluded_Metals"| Exclusion.CAS =="Current_Use")%>% droplevels(.)
  sumtmp_1 <- sumtmp %>% filter(Exclusion.CAS=="Excluded_Metals"| Exclusion.CAS =="Current_Use")%>% droplevels(.)%>%mutate(used=used)
  tmp_1 <- tmp %>% filter(Exclusion.CAS=="Excluded_Metals"| Exclusion.CAS =="Current_Use")%>% droplevels(.)%>%mutate(used=used)
  N_CAS_Site_1  <- N_CAS_Site_1 %>% filter(Exclusion.CAS=="Excluded_Metals"| Exclusion.CAS =="Current_Use")%>% droplevels(.)%>%mutate(used=used)
  return(list(tmp_1=tmp_1,N_CAS_Site_1=N_CAS_Site_1,sumtmp_1=sumtmp_1))

}

#' wrapper function to generate figure 3
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param LOQ.type
#' @param Exclusion.CAS
#' @param Data.SSD.1
#' @param threshold default 1, HI > 1
#' @param include0  whether to include HI 0?
#' @param threshold_current  the current use data, if we need use a smaller threshold for HI, default is 0.1
#' @param cutoff N measured > cutoff
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
gen_fig3_wrapper <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0"), Exclusion.CAS = c("Refined_HQ_Metals_PAH", "Excluded_Metals","Current_Use"),
                             Data.SSD.1,threshold = 1,include0=F,threshold_current=1,cutoff=10,addTag=T,...){
  driverdata <- getDriverData(Data.AggBySiteID.2,chemclass,LOQ.type = LOQ.type,
                              Exclusion.CAS = Exclusion.CAS, Data.SSD.1=Data.SSD.1,threshold = threshold ,
                              include0=include0,cutoff=cutoff,...)
  sumdriver <- driverdata %>% group_by(LOQ.type,Exclusion.CAS,test,Component) %>% summarise(HQ=sum(HQ)) %>% ungroup %>% group_by(LOQ.type,Exclusion.CAS,test) %>% mutate(sumHQ=sum(HQ),Perc=HQ/sumHQ) %>% mutate(Component=factor(Component,levels=c("Driver 1","Driver 2","Driver 3","All Other" )))
  driverdata1 <- driverdata %>% group_by(LOQ.type,Exclusion.CAS,test) %>% nest()%>%
    mutate(data=purrr::map(data,function(df){
      df1<- df%>%group_by(Site.Year,monitoringSiteIdentifier,phenomenonTimeReferenceYear,HI) %>% nest()

      df <- df1[order(df1$HI,decreasing = T),]%>% add_column(x=1:nrow(df1))%>%unnest(cols = c(data))%>% ungroup
      return(df)}))%>%
    unnest(cols = c(data))%>%mutate(Component=factor(Component,levels=c("Driver 1","Driver 2","Driver 3","All Other" )))

  driverdata1 <- driverdata1 %>% mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))
  #site.year <- driverdata%>%filter(Exclusion.CAS=="Excluded_Metals",test=="Chronic") %>%
   # group_by(LOQ.type)%>% nest()%>%mutate(sy=map(data,function(df)unique(df$Site.Year)))

  site.year <- driverdata%>%filter(Exclusion.CAS=="Excluded_Metals",test=="Chronic") %>%
    group_by(LOQ.type)%>% nest()%>%mutate(sy=purrr::map(data,function(df)unique(df$Site.Year)))

  ## extract driver data with the same site.year.

  ## threshold_current <- 1 ## Do we want re-generate the currentdat????
  if(threshold_current<threshold){
    currentdat <- getDriverData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0"),
                                Exclusion.CAS = c("Current_Use"),Data.SSD.1=Data.SSD.1, threshold=threshold_current,
                                include0=include0,cutoff=cutoff,...) ### different from the original driver data with threshold being 0.1!!!
    currentdat <- currentdat %>% filter(test=="Chronic")
    ## take only the samples same as in Excluded_Metals case!! Note that in this case current_use part would include more samples in the plotted figures compared to the original driverdata. !!!
    tmp1 <- currentdat %>% filter(LOQ.type=="LOQ.T0") %>% filter(Site.Year %in% site.year$sy[[1]])


    currentdat1 <- tmp1 %>%group_by(LOQ.type) %>% nest()%>% mutate(data=purrr::map(data,function(df){
      df1<- df%>%group_by(Site.Year,monitoringSiteIdentifier,phenomenonTimeReferenceYear,HI) %>% nest()
      df <- df1[order(df1$HI,decreasing = T),]%>% add_column(x=1:nrow(df1))%>%unnest(cols = c(data))%>% ungroup
      return(df)}))%>%unnest(cols = c(data))%>%mutate(Component=factor(Component,levels=c("Driver 1","Driver 2","Driver 3","All Other" )))

    driverdata2 <- rbind(driverdata1 %>% filter(Exclusion.CAS!="Current_Use" & test=="Chronic"),currentdat1)
  }else{
    driverdata2 <- driverdata1 %>% filter(test=="Chronic")
  }

  driverdata2 <- driverdata2 %>%
    mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))


  nsample2 <- driverdata2 %>% group_by(LOQ.type,Exclusion.CAS,test)%>% nest()%>%
    mutate(sy=purrr::map(data,function(df) length(unique(df$Site.Year))),hmax=purrr::map(data,function(df) max(df$HI)))%>%
    select(-data) %>% unnest(cols=c(sy,hmax)) %>%
    mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))


  sumdriver2 <- driverdata2 %>% group_by(LOQ.type,Exclusion.CAS,test,Site.Year) %>%
    mutate(sumHQ=sum(HQ))%>% ungroup %>%
    group_by(LOQ.type,Exclusion.CAS,test,Component,Site.Year) %>%
    summarise(Perc=HQ/sumHQ) %>%
    mutate(Component=factor(Component,levels=c("Driver 1","Driver 2","Driver 3","All Other" )))

  fig3 <- generate_fig3(driverdata2,nsample2,sumdriver2 = sumdriver2)

  if(addTag){
    fig3_tagdata <- data.frame(LOQ.type=LOQ.type,Exclusion.CAS=Exclusion.CAS)%>%
      mutate(Exclusion.CAS=factor(Exclusion.CAS,levels=c("Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")))
    fig3_tagdata$x = nsample2$sy*0.05
    fig3_tagdata$y <- nsample2$hmax*0.99
    fig3_tagdata$Tag <- LETTERS[1:3]

    fig3 <- fig3+geom_text(data=fig3_tagdata,aes(x=x,y=y,label=Tag),size=6)
  }
  return(fig3)
}

#' Generate Annotated Figure 3
#'
#' @param driverdata2 driver data
#' @param nsample2 summarised driver data
#' @param LOQ.type default LOQ.T0
#'
#' @return ggplot object
#' @export
#'
#' @examples
generate_fig3 <- function(driverdata2,nsample2,LOQ.type="LOQ.T0",sumdriver2){
  pbar2 <- ggplot(driverdata2%>%filter(test=="Chronic"),aes(x=x,y=HQ))+geom_bar(aes(fill=Component,col=Component),position=position_stack(reverse = T),stat = "identity")+facet_wrap(LOQ.type~Exclusion.CAS,scale="free")+ theme(axis.text.x = element_blank(), axis.ticks = element_blank())+geom_text(data=nsample2 %>% filter(test=="Chronic"),aes(x=sy,y=hmax,label=paste0("N.Sample = ",sy)),hjust=1)+xlab("")+theme(legend.position = "bottom")+geom_hline(yintercept = 1,lty=2)+ylab("HI")




  inset_plot <- NULL
  ## tmp <- expand.grid(LOQ.type = c("LOQ.T0","LOQ.T1"),Exclusion.CAS = c("Excluded_Metals_PAH","Current_Use"))
  for(i in 1:1)
    for(j in 1:3){
      inset_plot[[2*(i-1)+j]] <- NA
    }

  names(inset_plot) <- paste0(nsample2$LOQ.type,"-",nsample2$Exclusion.CAS)

  for (LOQ in LOQ.type) {
    for(Exclusion in c("Refined_HQ_Metals_PAH","Excluded_Metals","Current_Use")){
      p <- ggplot(data=sumdriver2 %>% filter(test=="Chronic" & Exclusion.CAS==Exclusion & LOQ.type==LOQ),
                  aes(x=Component,y=Perc,fill=Component))+
        geom_boxplot()  +
        guides(fill=FALSE) +
        theme_bw(base_size=12) +  ## makes everything smaller
        scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
        theme(panel.background = element_rect(fill="white"),  ## white plot background
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=rel(1),angle = 45,hjust=1,vjust=1), ## tiny axis text
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank())

      inset_plot[[paste0(LOQ,"-",Exclusion)]] <- p
    }
  }



  fig3 <- pbar2+annotation_custom2(grob = ggplotGrob(inset_plot[[1]]), data = data.frame(x=nsample2$sy[1], HQ=nsample2$hmax[1]*0.95,
                                                                                         LOQ.type=nsample2$LOQ.type[1],
                                                                                         Exclusion.CAS=nsample2$Exclusion.CAS[1]),
                                   xmin = nsample2$sy[1]*0.15, xmax = nsample2$sy[1]*0.95, ymin = nsample2$hmax[1]*0.15, ymax =nsample2$hmax[1]*0.95)+
    annotation_custom2(grob = ggplotGrob(inset_plot[[2]]), data = data.frame(x=nsample2$sy[2], HQ=nsample2$hmax[2]*0.95,                                        LOQ.type=nsample2$LOQ.type[2],
                                                                             Exclusion.CAS=nsample2$Exclusion.CAS[2]),
                       xmin = nsample2$sy[2]*0.15, xmax = nsample2$sy[2]*0.95, ymin = nsample2$hmax[2]*0.15, ymax =nsample2$hmax[2]*0.95)+
    annotation_custom2(grob = ggplotGrob(inset_plot[[3]]), data = data.frame(x=nsample2$sy[3], HQ=nsample2$hmax[3]*0.95,                                        LOQ.type=nsample2$LOQ.type[3],
                                                                             Exclusion.CAS=nsample2$Exclusion.CAS[3]),
                       xmin = nsample2$sy[3]*0.15, xmax = nsample2$sy[3]*0.95, ymin = nsample2$hmax[3]*0.15, ymax =nsample2$hmax[3]*0.95)


  return(fig3)
}




#' Wrapper function to generate figure 4a
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param LOQ.type
#' @param Exclusion.CAS
#' @param Data.SSD.1
#' @param threshold
#' @param saveplot
#' @param writetable
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
gen_fig4a_wrapper <- function(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0"),
                              Exclusion.CAS = c("Refined_HQ_Metals_PAH","Excluded_Metals"),
                              Data.SSD.1,threshold = 1,saveplot=FALSE,writetable=FALSE,...){
  driverdata <- getDriverData(Data.AggBySiteID.2,chemclass,LOQ.type =LOQ.type,
                              Exclusion.CAS = Exclusion.CAS, Data.SSD.1=Data.SSD.1,threshold = threshold,
                              addcaption="",...)

  sumCAS <- driverdata %>% filter(HI>1) %>% group_by(LOQ.type,Exclusion.CAS,test)%>%mutate(nSample=length(unique(Site.Year))) %>% group_by(LOQ.type,Exclusion.CAS,test,CAS,nSample) %>%nest() %>% mutate(nappear=purrr::map(data,nrow))%>% unnest(cols=c(nappear))%>% mutate(freq=nappear/nSample ) %>% select(-data) %>% group_by(LOQ.type,Exclusion.CAS,test)%>% nest() %>% mutate(data=purrr::map(data,function(df){
    df <- df[order(df$freq,decreasing = T),]
    df <- df%>%add_column(x=1:nrow(df))
    return(df)
  }))%>%unnest(cols = c(data))%>%ungroup #%>% mutate(x=paste0(x,"-",CAS))

  sumCAS <- sumCAS %>% left_join(chemclass[,c("CAS","Substance","Chem.Group","class")]%>%mutate(CAS=as.character(CAS)))%>%filter(CAS!="")%>% filter(x<=41)

  fig4 <- list()
  for(excludecas in Exclusion.CAS){
    fig4a <- generate_fig4a(sumCAS = sumCAS, Exclusion.CAS = excludecas,chemclass=chemclass)
    fig4a <- fig4a+theme(strip.text.x = element_text(size=12),axis.text=element_text(size=12,hjust=1),legend.direction = "horizontal")
    if(saveplot){
      fig4a
      ggsave(paste0("inst/manuscript_2022/Figure4a-Chronic-",excludecas,addcaption,".png"),width = 10,height=8,dpi=300)
      cairo_ps(paste0("inst/manuscript_2022/Figure4a-Chronic-",excludecas,addcaption,".eps"),width = 10,height=8)
      fig4a
      dev.off()
    }
    if(writetable) write.csv(sumCAS%>% mutate(Sub1=paste0(x,freq,"*",Substance)) %>% mutate(Sub1=factor(Sub1,levels=rev(unique(Sub1))))%>%filter(Exclusion.CAS==excludecas &test=="Chronic")%>% select(c(test,LOQ.type,Exclusion.CAS,freq,Substance)),
                             file=paste0("Driver_",excludecas,addcaption,".csv"))
  fig4[[excludecas]] <- fig4a
  }
  return(fig4)

}
#' Function to generate figure 4a
#'
#' @param sumCAS
#' @param Exclusion.CAS
#' @param chemclass
#'
#' @return
#' @export
#'
#' @examples
generate_fig4a <- function(sumCAS,Exclusion.CAS="Refined_HQ_Metals_PAH",chemclass,facet="LOQ.type"){
  sumCAS <- sumCAS%>% mutate(Sub1=paste0(x,freq,"*",Substance)) %>% mutate(Sub1=factor(Sub1,levels=rev(unique(Sub1))))


  color = c(Current_Use = "#2596be", Banned_Listed = "#f28e2b", Metal = "#e15759")
  library(ggtext)
  library(glue)
  labels <- data.frame(Substance=str_split(levels(sumCAS$Sub1),"\\*",simplify=T)[,2]) %>% left_join(chemclass[,c("Substance","class")]) %>%mutate(Sub2= glue("<i style='color: {color[class]}'>{Substance}</i>"))


  library(ggtext)
  excludecas <- Exclusion.CAS
  if(facet=="LOQ.type") ggplot(sumCAS%>%filter(Exclusion.CAS==excludecas &test=="Chronic"),aes(x=Sub1,y=freq,fill=Chem.Group))+geom_bar(stat="identity")+ coord_flip()+facet_wrap(~LOQ.type,scale="free",drop=TRUE)+scale_x_discrete(breaks=levels(sumCAS$Sub1),label=labels$Sub2)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))+xlab("")+  theme(axis.text.x = element_markdown(),axis.text.y = element_markdown(),legend.position = "bottom")+ylab("Frequency in top 3 drivers") else{
    sumCAS$facet <- facet
    ggplot(sumCAS%>%filter(Exclusion.CAS==excludecas &test=="Chronic"),aes(x=Sub1,y=freq,fill=Chem.Group))+geom_bar(stat="identity")+ coord_flip()+facet_wrap(~facet,scale="free",drop=TRUE)+scale_x_discrete(breaks=levels(sumCAS$Sub1),label=labels$Sub2)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))+xlab("")+  theme(axis.text.x = element_markdown(),axis.text.y = element_markdown(),legend.position = "bottom")+ylab("Frequency in top 3 drivers")

  }
}



#' Function to generate figure 4b
#'
#' @param Data.AggBySiteID.2
#' @param chemclass
#' @param Data.SSD.1
#' @param LOQ.type
#' @param driverdata
#' @param Exclusion.CAS
#'
#' @return
#' @export
#'
#' @examples
generate_fig4b <- function(Data.AggBySiteID.2,chemclass,Data.SSD.1,LOQ.type= c("LOQ.T0"),driverdata,Exclusion.CAS="Refined_HQ_Metals_PAH"){
  HQdata <- getHQData_Chronic(Data.AggBySiteID.2,chemclass=chemclass,Data.SSD.1=Data.SSD.1,LOQ.type= c("LOQ.T0"),driverdata,Exclusion.CAS=Exclusion.CAS)
  sumCAS <- HQdata %>% mutate(test="Chronic") %>%group_by(LOQ.type,test)%>%mutate(nSample=length(unique(Site.Year))) %>% group_by(LOQ.type,test,CAS,nSample) %>%nest() %>% mutate(nappear=purrr::map(data,nrow),QI=purrr::map(data,function(df) quantile(df$HQ.Chronic,c(0.25,0.5,0.75,0.95,1))))%>% unnest(cols=c(nappear,QI))%>% mutate(freq=nappear/nSample ) %>% select(-data) %>% mutate(Quantile=c("25%","50%","75%","95%","100%"))%>% group_by(LOQ.type,test,CAS)%>% mutate(q95=QI[Quantile=="95%"])%>% filter(q95>0.01)%>% group_by(LOQ.type,test)%>%nest()%>%mutate(data=purrr::map(data,function(df){
    df <- df %>% group_by(CAS,q95) %>% nest()
    df <- df[order(df$q95,decreasing=TRUE),]
    df <- df %>% add_column(x=1:nrow(df))
    df <- df %>% unnest(cols=c(data))
    return(df)
  }))%>%unnest(cols = c(data))%>%ungroup #%>% mutate(x=paste0(x,"-",CAS))

  sumCAS <- sumCAS %>% left_join(casdata[,c("CAS","Substance","Chem.Group","Nmeasured","Ndetected","LOQ","LOQ.Chronic","LOQ.Acute")]%>%mutate(CAS=as.character(CAS)))%>%filter(CAS!="") %>% filter(x<=41)

  write.csv(sumCAS,file=paste0("inst/manuscript_v3/table_fig4b_",Exclusion.CAS,".csv"))
  sumCAS<- sumCAS%>%filter(Quantile!="100%")
  sumCAS <- sumCAS%>%mutate(Sub1=paste0(x,"*",Substance)) %>% mutate(Sub1=factor(Sub1,levels=rev(unique(Sub1))),Sub2=paste0(Sub1,"|",Nmeasured,"|",Ndetected))
  ## To add min(LOQ) per CAS to the plot with a different symbol


  ## https://github.com/wilkelab/ggtext/issues/27
  # color = c(Current_Use = "#006400", Banned_Listed = "#696969", Metal = "#e15759")
  color = c(Current_Use = "#2596be", Banned_Listed = "#f28e2b", Metal = "#e15759")
  library(ggtext)
  labels <- data.frame(Substance=str_split(levels(sumCAS$Sub1),"\\*",simplify=T)[,2]) %>% left_join(casdata[,c("Substance","class","Nmeasured","Ndetected")])%>%mutate(Substance=paste0(Substance,"|",Nmeasured,"|",Ndetected)) %>%mutate(Sub2= glue::glue("<i style='background-color: {color[class]}'>{Substance}</i>")) %>% mutate(Sub3= glue::glue("<i style='color: {color[class]}'>{Substance}</i>"))


  fig4b <- ggplot(sumCAS%>%filter( test=="Chronic" & LOQ.type=="LOQ.T0" ),aes(x=Sub1,y=QI,fill=Chem.Group))+geom_point(aes(col=Chem.Group,pch=Quantile))+ coord_flip()+facet_wrap(LOQ.type~test,scale="free",drop=TRUE)+scale_x_discrete(breaks=levels(sumCAS$Sub1),label=labels$Sub3)+xlab("")+  theme(axis.text.x = element_markdown(size=12),axis.text.y = element_markdown(size=12,hjust = 1),legend.position = "bottom")+scale_y_log10()+ylab("HQ")##+geom_point(aes(y=LOQ.Chronic),pch="l")
  # if(interactive()){
  #   ggsave(paste0("inst/manuscript_v3/Figure4b_",Exclusion.CAS,".png"),width = 12,height=9,dpi=300)
  # }
  return(fig4b)
}
