library(tidyverse)
library(coxme)
load("data_Methyl_clean/eQTM.data.RData")
results <- data.frame()

for(gene.OI in unique(rna.overlap$gene)){
  print(gene.OI) 
  #Filter expression to gene of interest
  rna.temp <- rna.overlap %>% 
    filter(gene == gene.OI)
  #### DMP model ####
  DMP.temp <- DMP %>% 
    filter(gene == gene.OI)
  
  result.DMP <- data.frame()
  #if gene has DMP
  if(nrow(DMP.temp)>0){
    for(probe.OI in unique(DMP.temp$probeID)){
      DMP.M.temp <- DMP.M.overlap %>% 
        filter(probeID == probe.OI)
      
      #Add to RNA data
      dat.temp <- DMP.M.temp %>% 
        inner_join(rna.temp, by = "FULLIDNO")
      
      ####MEDIA model####
      dat.media <- dat.temp %>% filter(condition=="MEDIA")
      
      fit <- lmekin(expression ~ methyl + KCHCA_AGE_YR_CURRENT + 
                      M0_KCVSEX + (1|FULLIDNO),
                    data=dat.media, varlist=as.matrix(kin.overlap))
      beta <- fit$coefficients$fixed
      nvar <- length(beta)
      nfrail <- nrow(fit$var) - nvar
      se <- sqrt(diag(fit$var)[nfrail + 1:nvar])
      t <- beta/se
      p <- 1 - pchisq((beta/se)^2, 1)
      sigma <- fit$sigma
      
      #####TB model####
      dat.tb <- dat.temp %>% filter(condition=="TB")
      
      fit2 <- lmekin(expression ~ methyl + KCHCA_AGE_YR_CURRENT + 
                       M0_KCVSEX + (1|FULLIDNO),
                     data=dat.tb, varlist=as.matrix(kin.overlap))
      beta2 <- fit2$coefficients$fixed
      nvar <- length(beta2)
      nfrail <- nrow(fit2$var) - nvar
      se <- sqrt(diag(fit2$var)[nfrail + 1:nvar])
      t2 <- beta2/se
      p2 <- 1 - pchisq((beta2/se)^2, 1)
      sigma2 <- fit2$sigma
      
      ####Results####
      result.DMP <- data.frame(
        epigen = rep(probe.OI,length(p)+length(p2)),
        model = c(rep("MEDIA",length(p)),rep("TB",length(p))),
        variable=c(names(p),names(p2)),
        pval = c(p,p2),
        sigma = c(rep(sigma,length(p)),rep(sigma2,length(p))),
        t = c(t,t2),
        beta=c(beta,beta2)) %>% 
        mutate(gene=gene.OI) %>% 
        bind_rows(result.DMP)
      
    }} else{
      result.DMP <- NULL
    }
  
  #### DMR model ####
  DMR.temp <- DMR %>% 
    filter(gene == gene.OI)
  
  result.DMR <- data.frame()
  #if gene has DMR
  if(nrow(DMR.temp)>0){
    for(DMR.OI in unique(DMR.temp$DMR)){
      DMR.M.temp <- DMR.M.overlap %>% 
        filter(DMR ==DMR.OI)
      
      #Add to RNA data
      dat.temp <- DMR.M.temp %>% 
        inner_join(rna.temp, by = "FULLIDNO")
      
      ####MEDIA model####
      dat.media <- dat.temp %>% filter(condition=="MEDIA")
      
      fit <- lmekin(expression ~ methyl*Sample_Group +
                      KCHCA_AGE_YR_CURRENT + 
                      M0_KCVSEX + (1|FULLIDNO),
                    data=dat.media, varlist=as.matrix(kin.overlap))
      beta <- fit$coefficients$fixed
      nvar <- length(beta)
      nfrail <- nrow(fit$var) - nvar
      se <- sqrt(diag(fit$var)[nfrail + 1:nvar])
      t <- beta/se
      p <- 1 - pchisq((beta/se)^2, 1)
      sigma <- fit$sigma
      
      ####TB model####
      dat.tb <- dat.temp %>% filter(condition=="TB")
      
      fit2 <- lmekin(expression ~ methyl*Sample_Group +
                       KCHCA_AGE_YR_CURRENT + 
                       M0_KCVSEX + (1|FULLIDNO),
                     data=dat.tb, varlist=as.matrix(kin.overlap))
      beta2 <- fit2$coefficients$fixed
      nvar <- length(beta2)
      nfrail <- nrow(fit2$var) - nvar
      se <- sqrt(diag(fit2$var)[nfrail + 1:nvar])
      t2 <- beta2/se
      p2 <- 1 - pchisq((beta2/se)^2, 1)
      sigma2 <- fit2$sigma
      
      ####results####
      result.DMR <- data.frame(
        epigen = rep(DMR.OI,length(p)+length(p2)),
        model = c(rep("MEDIA",length(p)),rep("TB",length(p))),
        variable=c(names(p),names(p2)),
        pval = c(p,p2),
        sigma = c(rep(sigma,length(p)),rep(sigma2,length(p))),
        t = c(t,t2),
        beta=c(beta,beta2)) %>% 
        mutate(gene=gene.OI) %>% 
        bind_rows(result.DMR)
    }} else {
      result.DMR <- NULL
    }
  
  results <- bind_rows(result.DMP,result.DMR,results)
  
}

### Save
results %>% 
  group_by(variable, model) %>% 
  mutate(FDR=p.adjust(pval)) %>% 
  ungroup() %>% 
  
  write_csv(file="results/RSTR.eQTM.csv")
