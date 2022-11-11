##### Packages #####
# install.packages(c("tidyverse","devtools"))
# This version must be installed in order for lmekin to work. Newer versions of kimma use lmerel
# devtools::install_github("BIGslu/kimma", "272ede4") 

library(tidyverse)
library(kimma)
set.seed(589)

#### Data ####
#log2 ratios of beta (M-values)
RSTR.M <- read_csv("data/RSTR_clean_M.csv.gz") %>% 
  column_to_rownames("probeID")

#Probe features
probes <- read_csv("data/RSTR_probe_features.csv.gz")

#Complete metadata for methylation samples
meta <- read_csv("data/2020.11.20RSTR_Hawn_metadata.csv", 
                 na=c("",NA)) %>% 
  #Keep only methylation samples
  filter(methylation == "PASS") %>% 
  #Make family and name variables characters
  mutate_at(c("family.3","family.2","family.1","Sample_Name"),
            ~as.character(.)) %>% 
  arrange(FULLIDNO)

#Kinship matrix
kin <- read_csv("data/kinship_Hawn_all.csv") %>% 
  #Keep methyl samples
  pivot_longer(-rowname) %>% 
  filter(rowname %in% meta$FULLIDNO & name %in% meta$FULLIDNO) %>% 
  pivot_wider() %>% 
  column_to_rownames()

## Check data order
identical(rownames(kin),colnames(kin))

meta.kin <- meta %>% 
  filter(FULLIDNO %in% rownames(kin))
identical(rownames(kin),meta.kin$FULLIDNO)

#### Loop setup ####
bin <- 100000
## Short names for outputs
name.key <- data.frame(
  model = c("~ Sample_Group", 
            "~ Sample_Group + (1|family.3)",
            "~ Sample_Group + (1|family.2)",
            "~ Sample_Group + (1|family.1)",
            "~ Sample_Group + (1|FULLIDNO)",
            "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + avgBMI + (1|FULLIDNO)",
            "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + RISK_SCORE + KCB_BCGSCAR + (1|FULLIDNO)",
            "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + (1|FULLIDNO)"
  ),
  name = c("base","family3","family2","family1",
           "kinship","kinship_demo","kinship_risk","kinship_final")
)
#
#### Models without kinship ####
models <- c("~ Sample_Group", 
            "~ Sample_Group + (1|family.3)",
            "~ Sample_Group + (1|family.2)",
            "~ Sample_Group + (1|family.1)")

result.ls <- NULL

for(m in models){
  print(m)
  result.ls <- NULL
  
  for(i in seq(1,nrow(RSTR.M),bin)){
    print(i)
    start <- Sys.time()
    #Recode end to not go beyond data
    if(i+bin-1 > nrow(RSTR.M)){ j <- nrow(RSTR.M) } else { j <- i+bin-1 }
    
    #Set model type
    if(grepl("\\|",m)){ 
      model.lme <- TRUE
      model.lm <- FALSE
    } else{
      model.lme <- FALSE
      model.lm <- TRUE
    }
    
    #Run model
    model_temp <- kmFit(counts = RSTR.M, meta = meta,
                        model = m,
                        run.lm = model.lm, run.lme = model.lme,
                        metrics = TRUE,
                        libraryID = "Sample_Name", 
                        patientID = "FULLIDNO", 
                        subset.genes = rownames(RSTR.M)[i:j])
    
    #Concatenate results
    if(i==1){
      result.ls <- model_temp
    } else if(i!=1 & grepl("\\|",m)){
      result.ls[["lme"]] <- result.ls[["lme"]] %>% 
        bind_rows(model_temp[["lme"]])
      result.ls[["lme.fit"]] <- result.ls[["lme.fit"]] %>% 
        bind_rows(model_temp[["lme.fit"]])
    } else {
      result.ls[["lm"]] <- result.ls[["lm"]] %>% 
        bind_rows(model_temp[["lm"]])
      result.ls[["lm.fit"]] <- result.ls[["lm.fit"]] %>% 
        bind_rows(model_temp[["lm.fit"]])
    }
    print(Sys.time()-start)
  }
  
  #Recalc FDR across all probes in each model
  if(grepl("\\|",m)){
    result.ls[["lme"]] <- result.ls[["lme"]] %>% 
      group_by(model, variable) %>%
      mutate(FDR=p.adjust(pval, method="BH")) %>%
      ungroup()
  } else {
    result.ls[["lm"]] <- result.ls[["lm"]] %>% 
      group_by(model, variable) %>%
      mutate(FDR=p.adjust(pval, method="BH")) %>%
      ungroup()
  }
  
  #Make nice names
  result.name <- name.key %>% 
    filter(model == m) %>% 
    pull(name)
  assign(result.name, result.ls)
  
  #Save
  save(list=result.name, 
       file = paste0("results/model_selection/",result.name,"_model.RData"))
}

#### Models with kinship ####
models2 <- c(
  "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + (1|FULLIDNO)",
  "~ Sample_Group + (1|FULLIDNO)",
  "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + avgBMI + (1|FULLIDNO)",
  "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + RISK_SCORE + KCB_BCGSCAR + (1|FULLIDNO)"
  )

result.ls <- NULL

for(m in models2){
  print(m)
  result.ls <- NULL
  
  for(i in seq(1,nrow(RSTR.M),bin)){
    print(i)
    start <- Sys.time()
    #Recode end to not go beyond data 
    if(i+bin-1 > nrow(RSTR.M)){ j <- nrow(RSTR.M) } else { j <- i+bin-1 }
    
    #Run model
    model_temp <- kmFit(counts = RSTR.M, meta = meta, kin = kin,
                        model = m,
                        run.lmerel = TRUE, # run.lmekin = TRUE,
                        metrics = TRUE,
                        libraryID = "Sample_Name", 
                        patientID = "FULLIDNO", 
                        subset.genes = rownames(RSTR.M)[i:j])
    
    #Concatenate results
    if(i==1){
      result.ls <- model_temp
    } else {
      result.ls[["lmekin"]] <- result.ls[["lmekin"]] %>% 
        bind_rows(model_temp[["lmekin"]])
      result.ls[["lmekin.fit"]] <- result.ls[["lmekin.fit"]] %>% 
        bind_rows(model_temp[["lmekin.fit"]])
    } 
    print(Sys.time()-start)
  }
  
  #Recalculate FDR across all probes in a model
  result.ls[["lmekin"]] <- result.ls[["lmekin"]] %>% 
    group_by(model, variable) %>%
    mutate(FDR=p.adjust(pval, method="BH")) %>%
    ungroup()
  
  #Make nice names
  result.name <- name.key %>% 
    filter(model == m) %>% 
    pull(name)
  assign(result.name, result.ls)
  
  #Save
  save(list=result.name, 
       file = paste0("results/model_selection/",result.name,"_model.RData"))
}
