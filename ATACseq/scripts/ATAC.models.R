library(tidyverse)
library(kimma)
set.seed(589)

###### Load data #######
load("data_ATACseq_clean/ATACseq_data.RData")

#Rename all to FULLIDNO
counts <- as.data.frame(voomQW.Nfree.abund.norm$E) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "RS_SUB_ACCESSION_NO") %>% 
  left_join(select(voomQW.Nfree.abund.norm$targets, RS_SUB_ACCESSION_NO, FULLIDNO)) %>% 
  select(-RS_SUB_ACCESSION_NO) %>% 
  arrange(FULLIDNO) %>% 
  pivot_wider(names_from = FULLIDNO) %>% 
  column_to_rownames()

kin.format <- kin %>% 
  #Keep methyl samples
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  left_join(select(voomQW.Nfree.abund.norm$targets, RS_SUB_ACCESSION_NO, FULLIDNO),
            by = c("rowname"="RS_SUB_ACCESSION_NO")) %>% 
  rename(rowname.new = FULLIDNO) %>% 
  left_join(select(voomQW.Nfree.abund.norm$targets, RS_SUB_ACCESSION_NO, FULLIDNO),
            by = c("name"="RS_SUB_ACCESSION_NO")) %>% 
  rename(name.new = FULLIDNO) %>% 
  select(-rowname, -name) %>% 
  arrange(name.new) %>% 
  pivot_wider(names_from = name.new) %>% 
  arrange(rowname.new) %>% 
  column_to_rownames("rowname.new")

#### Loop setup ####
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

  #Set model type
  if(grepl("\\|",m)){ 
    model.lme <- TRUE
    model.lm <- FALSE
  } else{
    model.lme <- FALSE
    model.lm <- TRUE
  }
  
  start <- Sys.time()
  #Run model
  result.ls <- kmFit(counts = counts, meta = voomQW.Nfree.abund.norm$targets,
                     kin = kin.format,
                     model = m,
                     run.lm = model.lm, run.lme = model.lme,
                     metrics = TRUE,
                     libraryID = "FULLIDNO", patientID = "FULLIDNO")
  print(Sys.time()-start)

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
  "~ Sample_Group + (1|FULLIDNO)",
  "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + avgBMI + (1|FULLIDNO)",
  "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + RISK_SCORE + KCB_BCGSCAR + (1|FULLIDNO)",
  "~ Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + (1|FULLIDNO)"
)

result.ls <- NULL

for(m in models2){
  print(m)
  result.ls <- NULL
  
    start <- Sys.time()
    #Run model
    result.ls <- kmFit(counts = counts, meta = voomQW.Nfree.abund.norm$targets,
                       kin = kin.format,
                       model = m,
                       run.lmerel = TRUE,
                       metrics = TRUE,
                       libraryID = "FULLIDNO", patientID = "FULLIDNO")
    print(Sys.time()-start)
  
  #Make nice names
  result.name <- name.key %>% 
    filter(model == m) %>% 
    pull(name)
  assign(result.name, result.ls)
  
  #Save
  save(list=result.name, 
       file = paste0("results/model_selection/",result.name,"_model.RData"))
}
