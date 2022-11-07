
#' Simulation function
#'
#' @param cvs A vector of CV
#' @param mean_sim  mean of the mixture components.
#' @param p mixing proportions, added up to 1
#'
#' @return a vector of MDRs
#' @export
#'
#' @examples
simCVn <- function(cvs,mean_sim=1,p=c(0.5,0.5),Nsim=1000,synergy=F,reduction=0.5){
  ##nCV <- length(cvs)
  res <- lapply(cvs,function(x) simCV1(cv=x,mean_sim=mean_sim,p=p,Nsim=Nsim,synergy=synergy,reduction=reduction))
  names(res) <- cvs
  res <- as.data.frame(res)
  require(tidyverse)
  res <- tidyr::pivot_longer(res,cols=everything(),names_to=c("CV"),values_to="MDR") %>%
    mutate(CV=gsub("X","",CV))
  return(res)
}

simCV1 <- function(cv,mean_sim=1,p=c(0.5,0.5),Nsim=1000,synergy=F,reduction=0.5){
## should be vectorized against cv.
  sd_sim <- cv*mean_sim
  v_sim <- sd_sim^2
  mu <- log(mean_sim/sqrt(1+v_sim/mean_sim^2))
  sigma <- sqrt(log(1+v_sim/mean_sim^2))
  if(synergy){
    mean_mix <- mean_sim*(1-reduction)
    mu_mix <- log(mean_mix/sqrt(1+v_sim/mean_mix^2))
    sigma_mix <- sqrt(log(1+v_sim/mean_mix^2))
  }

  ## Values were generated for each toxicity test that would be
  ##  conducted in a mixture study.  F
  ## For example, for a binary mixture AB, a value would be generated for compound A,
  ## compound B, and mixture AB.
  nmix <- length(p)
  MDRs <- sapply(1:Nsim,function(i){
    if(synergy==F) ECx <- rlnorm(nmix+1,meanlog = mu,sdlog=sigma) else{
      ECx <- c(rlnorm(nmix,meanlog = mu,sdlog=sigma),rlnorm(1,meanlog = mu_mix,sdlog=sigma_mix))
    }
    ECx_mix <- 1/(sum(p/ECx[1:nmix]))
    MDR <- ECx_mix/(ECx[1+nmix])
    MDR
  })
  return(MDRs)
}

#' Simulate MDR matrix
#'
#' @param cvs
#' @param mean_sim
#' @param p
#' @param Nsim
#' @param M
#'
#' @return
#' @export
#'
#' @examples
simCVmat <- function(cvs,mean_sim=1,p=c(0.5,0.5),Nsim=1000,M=100,synergy=F,reduction=0.5){
  res <- lapply(1:M,function(i){
    res1 <- simCVn(cvs=cvs,mean_sim=mean_sim,p=p,Nsim=Nsim,synergy=synergy,reduction=reduction)
  })
  names(res) <- paste0("Rep_",1:M)
  res2 <- lapply(res,function(x)x$MDR)
  df <- as.data.frame(res2)
  df$CV <- res[[1]]$CV
  return(df)
}



#' Simulate the confusion matrix under different assumptions
#'
#' @param cvs vector of cv
#' @param mdrs vector of mdr
#' @param Nsim number of simulation
#' @param synergy logical
#' @param reduction reduction factor
#'
#' @return the confusion matrix
#' @export
#'
#' @examples
simConfusion <- function(cvs=c(0.3,0.6,1,1.4,3),
                         mdrs=c(1.25,2,3,5),
                         Nsim=12000,
                         synergy=F,reduction=0.5,q=0.95){

  res1 <- simCVn(cvs=cvs,mean_sim=1,p=1,Nsim=Nsim,synergy = synergy,reduction=reduction)

  c1 <-res1%>%group_by(CV)%>% nest() %>% mutate(confusion=map(data,function(df) sapply(mdrs,function(x) sum(df>x)/Nsim)))%>% mutate(Mixture=1) %>% dplyr::select(-data) %>% unnest(cols=c(confusion)) %>% ungroup %>% mutate(cutoff=rep(mdrs,length(cvs))) %>% pivot_wider(names_from = cutoff,names_prefix="MDR>",values_from=confusion)

  r1 <- res1%>%group_by(CV)%>% summarise(mean=mean(MDR),Q=quantile(MDR,q)) %>% mutate(Mixture=1)


  res2 <- simCVn(cvs=cvs,mean_sim=1,p=c(0.5,0.5),Nsim=Nsim,synergy = synergy,reduction=reduction)

  c2 <-res2%>%group_by(CV)%>% nest() %>% mutate(confusion=map(data,function(df) sapply(mdrs,function(x) sum(df>x)/Nsim)))%>% mutate(Mixture=2) %>% dplyr::select(-data) %>% unnest(cols=c(confusion)) %>% ungroup %>% mutate(cutoff=rep(mdrs,length(cvs))) %>% pivot_wider(names_from = cutoff,names_prefix="MDR>",values_from=confusion)

  r2 <- res2%>%group_by(CV)%>% summarise(mean=mean(MDR),Q=quantile(MDR,q))%>% mutate(Mixture=2)

  res3 <- simCVn(cvs=cvs,mean_sim=1,p=c(1/3,1/3,1/3),Nsim=Nsim,synergy = synergy,reduction=reduction)


  c3 <-res3%>%group_by(CV)%>% nest() %>% mutate(confusion=map(data,function(df) sapply(mdrs,function(x) sum(df>x)/Nsim)))%>% mutate(Mixture=3) %>% dplyr::select(-data) %>% unnest(cols=c(confusion)) %>% ungroup %>% mutate(cutoff=rep(mdrs,length(cvs))) %>% pivot_wider(names_from = cutoff,names_prefix="MDR>",values_from=confusion)

  r3 <- res3%>%group_by(CV)%>% summarise(mean=mean(MDR),Q=quantile(MDR,q))%>% mutate(Mixture=3)
  confusion <- rbind(c1,c2,c3)
  confusion <- confusion %>% pivot_longer(cols = starts_with("MDR"))
  #ggplot(confusion,aes(x=name,y=value,col=CV))+geom_point()+geom_line(aes(group=CV))+facet_grid(.~Mixture)+theme(axis.text.x = element_text(angle = 90))+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))+geom_hline(yintercept = 0.05,lty=2)+xlab("MDR Cutoff")
  res <- rbind(r1,r2,r3)
  return(list(confusion=confusion,res=res))

}



#' Simulate the confusion matrix under different assumptions
#'
#' @param cvs vector of cv
#' @param mdrs vector of mdr
#' @param Nsim number of simulation
#' @param synergy logical
#' @param reduction reduction factor
#'
#' @return the confusion matrix
#' @export
#'
#' @examples
simConfusion_2mixture <- function(cvs=c(0.3,0.6,1,1.4,3),
                         mdrs=c(1.25,2,3,5),
                         Nsim=12000,
                         synergy=F,reduction=0.5,q=0.95){



  res2 <- simCVn(cvs=cvs,mean_sim=1,p=c(0.5,0.5),Nsim=Nsim,synergy = synergy,reduction=reduction)

  c2 <-res2%>%group_by(CV)%>% nest() %>% mutate(confusion=map(data,function(df) sapply(mdrs,function(x) sum(df>x)/Nsim)))%>% mutate(Mixture=2) %>% dplyr::select(-data) %>% unnest(cols=c(confusion)) %>% ungroup %>% mutate(cutoff=rep(mdrs,length(cvs))) %>% pivot_wider(names_from = cutoff,names_prefix="MDR>",values_from=confusion)

  r2 <- res2%>%group_by(CV)%>% summarise(mean=mean(MDR),Q=quantile(MDR,q))%>% mutate(Mixture=2)

  confusion <- c2
  confusion <- confusion %>% pivot_longer(cols = starts_with("MDR")) %>% mutate(MDR = as.numeric(gsub("MDR>","",name)))
  #ggplot(confusion,aes(x=name,y=value,col=CV))+geom_point()+geom_line(aes(group=CV))+facet_grid(.~Mixture)+theme(axis.text.x = element_text(angle = 90))+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))+geom_hline(yintercept = 0.05,lty=2)+xlab("MDR Cutoff")
  res <- r2
  return(list(confusion=confusion,res=res))

}
