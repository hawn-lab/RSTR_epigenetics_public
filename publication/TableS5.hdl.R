library(tidyverse)
library(readxl)

#### Stats ####
fdr <- read_csv("lipoprotein/results/lipoprot_efflux_model_results.csv") %>% 
  select(gene, variable, estimate, pval, FDR.new) %>% 
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
  mutate(FDR.new = ifelse(grepl("HDL",gene), FDR.new, NA)) %>% 
  rename(FDR=FDR.new, outcome=gene) %>% 
  # arrange(group, FDR) %>% 
  select(group, everything())

#### Combine and save ####
write_csv(fdr, file = "publication/TableS5.HDL.efflux.csv")
