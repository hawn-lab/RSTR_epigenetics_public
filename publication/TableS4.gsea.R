library(tidyverse)
library(openxlsx)

load("results/enrichment/MDM.IFN.gsea.RData")
load("results/enrichment/MTB.RSTR.gsea.RData")

df.ls <- list()

df.ls[["MDM-IFN-GSEA"]] <- gsea_IFN_c5_signif %>% 
  select(group, pathway, pval, NES, leadingEdge) %>% 
  mutate(leadingEdge = as.character(leadingEdge)) %>% 
  arrange(pathway, group)

df.ls[["MTB-RSTR-GSEA"]] <- gsea_MTB_c5_signif %>% 
  select(group, pathway, pval, NES, leadingEdge) %>% 
  mutate(leadingEdge = as.character(leadingEdge)) %>% 
  arrange(pathway, group)

write.xlsx(df.ls, file = "../publication/TableS4.gsea.xlsx")
