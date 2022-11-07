# chemclass <- read.csv("data-raw/WB-WQ-CAS-CLASIFICATION.csv")
chemclass <- read.csv("data-raw/WB-WQ-CAS-CLASIFICATION-08122021.csv") ## Corrected classification information

names(chemclass)
table(chemclass$Chem.Group)
chemclass <- chemclass %>% mutate(class=getclass(Chem.Group,Approved,Priority,Current_Use))
chemclass$Approved[chemclass$Current_Use=="YES"] <- "Y"

chemclass$CAS <- as.character(chemclass$CAS)



## Note that SVHC
## Candidate List of substances of very high concern for Authorisation

PAH_PNEC <- read.csv("data-raw/PAHs_PNECs_Petrotox2021.csv")
which(PAH_PNEC$PAHs[1:2] %in% chemclass$Substance)
which(!(PAH_PNEC$CAS %in% chemclass$CAS.y))
which(chemclass$Substance==PAH_PNEC$PAHs[1])

which(chemclass$Substance==PAH_PNEC$PAHs[2])

chemclass <- chemclass%>% left_join(PAH_PNEC,by=c("CAS.y"="CAS")) %>% dplyr::select(-c("Source","Ref"))
chemclass[c(179,163),"PNEC"] <- PAH_PNEC[1:2,"PNEC"]

chemclass%>% filter(Chem.Group=="PAH")%>%select(c(Substance,PNEC)) %>% filter(is.na(PNEC))
which(chemclass$Substance=="1,2,4-Trimethylbenzene")
which(chemclass$Substance=="2-Methylnaphthalene")
chemclass[c(121,161),]
chemclass[242,"Type"] <- "I"
usethis::use_data(chemclass,overwrite = T)


Metal_ref <- read.csv("data-raw/Metals_Background_FOREG_07292021.csv")
usethis::use_data(Metal_ref,overwrite = T)
all(Metal_ref$CAS %in% chemclass$CAS)
all(Metal_ref$CAS.y %in% chemclass$CAS.y)
all(Metal_ref$Substance %in% chemclass$Substance)
Metal_ref$CAS <- as.character(Metal_ref$CAS)
chemclass <- chemclass%>% left_join(Metal_ref%>%dplyr::select(-c("Unit","Media","Source","Ref")))
chemclass%>% filter(Chem.Group=="Metal")
chemclass$class1 <- chemclass$class
chemclass$class1[chemclass$Chem.Group=="PAH"] <- "PAH"
usethis::use_data(chemclass,overwrite = T)



PAHs <- chemclass %>% filter(Chem.Group=="PAH") %>% filter(!is.na(PNEC))
ssdInfo <- Data.SSD.1[,c("CAS.","X10LogSSDMedianConcentration.ug.L..MuAcute.EC50","X10LogSSDMedianConcentration.ug.L..MuChronic.NOEC", "X10LogSSDSlope.ug.L..SigmaAcute.EC50","X10LogSSDSlope.ug.L..SigmaChronic.NOEC")] %>% mutate(CAS.=as.character(CAS.))
PAHs <- PAHs %>% left_join(ssdInfo,by=c("CAS"="CAS."))
data.table::setnames(PAHs,old =c("X10LogSSDMedianConcentration.ug.L..MuAcute.EC50","X10LogSSDMedianConcentration.ug.L..MuChronic.NOEC", "X10LogSSDSlope.ug.L..SigmaAcute.EC50","X10LogSSDSlope.ug.L..SigmaChronic.NOEC"),
new=c("SSDLOG10.Mu.Acute.EC50", "SSDLOG10.Mu.chronic.NOEC", "SSDLOG10.Sigma.Acute.EC50","SSDLOG10.Sigma.chronic.NOEC")
)
PAHs$HC05.Chronic <- 10^qnorm(0.05,PAHs$SSDLOG10.Mu.chronic.NOEC, PAHs$SSDLOG10.Sigma.chronic.NOEC)
usethis::use_data(PAHs,overwrite = T)

# All data given in ug/L
# For Metals and PAHs
# PAHs = replace HC05 by PNEC
# Metals= subtract background from given concentration


# chemclass with all LOQ information and HC05.
## Add information LOQ summaryN and also HQ (Chronic).
library(broom)
casdata <- getCasData(Data.AggBySiteID.2,chemclass,Data.SSD.1,Exclusion.CAS = "Excluded_None")
casdata$Type <- factor(casdata$Type,levels=c("I","H","F",""))
casdata$Chem.Group <- factor(casdata$Chem.Group,levels=c("Industrial", "Metal", "PAH", "Pest", "Pharma"))

usethis::use_data(casdata,overwrite = T)
names(casdata)
write.csv(casdata[,c(41,42,2:40,64,65)],"inst/SM/SM3_casdata_rawConc.csv")
