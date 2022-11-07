driverdata <- getDriverData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0"), Exclusion.CAS = c("Refined_Metals_PAH","Excluded_Metals"), Data.SSD.1,threshold = 1)

sumCAS <- driverdata %>% filter(HI>1) %>% group_by(LOQ.type,Exclusion.CAS,test)%>%mutate(nSample=length(unique(Site.Year))) %>% group_by(LOQ.type,Exclusion.CAS,test,CAS,nSample) %>%nest() %>% mutate(nappear=purrr::map(data,nrow))%>% unnest(cols=c(nappear))%>% mutate(freq=nappear/nSample ) %>% select(-data) %>% group_by(LOQ.type,Exclusion.CAS,test)%>% nest() %>% mutate(data=purrr::map(data,function(df){
  df <- df[order(df$freq,decreasing = T),]
  df <- df%>%add_column(x=1:nrow(df))
  return(df)
}))%>%unnest(cols = c(data))%>%ungroup #%>% mutate(x=paste0(x,"-",CAS))

sumCAS <- sumCAS %>% left_join(chemclass[,c("CAS","Substance","Chem.Group","class")]%>%mutate(CAS=as.character(CAS)))%>%filter(CAS!="")%>% filter(x<=41)

sumCAS_mean <- sumCAS

driverdata_max <- getDriverData(Data.AggBySiteID.2,chemclass,LOQ.type = c("LOQ.T0"), Exclusion.CAS = c("Refined_Metals_PAH","Excluded_Metals"), Data.SSD.1,threshold = 1,useAllMax=T)

sumCAS <- driverdata_max %>% filter(HI>1) %>% group_by(LOQ.type,Exclusion.CAS,test)%>%mutate(nSample=length(unique(Site.Year))) %>% group_by(LOQ.type,Exclusion.CAS,test,CAS,nSample) %>%nest() %>% mutate(nappear=purrr::map(data,nrow))%>% unnest(cols=c(nappear))%>% mutate(freq=nappear/nSample ) %>% select(-data) %>% group_by(LOQ.type,Exclusion.CAS,test)%>% nest() %>% mutate(data=purrr::map(data,function(df){
  df <- df[order(df$freq,decreasing = T),]
  df <- df%>%add_column(x=1:nrow(df))
  return(df)
}))%>%unnest(cols = c(data))%>%ungroup #%>% mutate(x=paste0(x,"-",CAS))

sumCAS <- sumCAS %>% left_join(chemclass[,c("CAS","Substance","Chem.Group","class")]%>%mutate(CAS=as.character(CAS)))%>%filter(CAS!="")%>% filter(x<=41)
sumCAS_max <- sumCAS

d1 <- sumCAS_mean%>%filter(freq>0.1 & test=="Chronic" & Exclusion.CAS=="Excluded_Metals")
d2 <- sumCAS_max%>%filter(freq>0.1& test=="Chronic" & Exclusion.CAS=="Excluded_Metals")
#sumCAS_mean$CAS
driver10 <- unique(c(d1$CAS,d2$CAS))
usethis::use_data(driver10)
d1 <- sumCAS_mean%>%filter(freq>0.05& test=="Chronic" & Exclusion.CAS=="Excluded_Metals")
d2 <- sumCAS_max%>%filter(freq>0.05& test=="Chronic" & Exclusion.CAS=="Excluded_Metals")
#sumCAS_mean$CAS
driver5 <- unique(c(d1$CAS,d2$CAS))
usethis::use_data(driver5)
