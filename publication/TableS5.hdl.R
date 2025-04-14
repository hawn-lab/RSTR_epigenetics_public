library(tidyverse)
library(readxl)

#### HDL in RSTR v LTBI ####
fdr <- read_csv("lipoprotein/results/lipoprot_efflux_model_results.csv") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  filter(gene != "J774.abca1") %>% 
  filter(! gene %in% c("msuHDL","ssHDL")) %>% 
  mutate(gene = recode(gene, 
                       "J774.ind"="J774_induced", 
                       "J774.basal"="J774_basal",
                       "ABCA1.ind"="BHK_induced", 
                       "ABCA1.basal"="BHK_basal",
                       "ABCA1.spec"="BHK_ABCA1_specific")) %>%
  mutate(group = case_when(grepl("sz", gene)~"HDL_size",
                           grepl("J774|BHK",gene)~"HDL_efflux",
                           TRUE~"HDL")) %>% 
  mutate(FDR = ifelse(grepl("HDL",gene), FDR, NA)) %>% 
  rename(outcome=gene) %>% 
  # arrange(group, FDR) %>% 
  select(group, outcome, variable, estimate, pval, FDR) %>% 
  #remove proportion results
  filter(!grepl("_pct$",outcome))

#### CEC vs HDL ####
fdr2 <- read_csv("lipoprotein/results/lipoprot_efflux-HDL_model_results.csv") %>% 
  #remove proportion results
  filter(!grepl("_pct$",model) & gene=="ABCA1.spec") %>% 
  mutate(group = "CEC_HDL") %>% 
  rename(outcome=gene) %>% 
  select(group, outcome, variable, estimate, pval, FDR)


#### Combine and save ####
library(openxlsx)

dfs <- list("HDL~RSTR"=fdr, "CEC~HDL"=fdr2)
write.xlsx(dfs, file = "publication/TableS5.HDL.efflux.xlsx")
