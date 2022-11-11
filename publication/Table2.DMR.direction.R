library(tidyverse)

#### DATA ####
#Genes in signif enrich terms
attach("Methyl/results/enrichment/RSTR.DMPR.enrichment.RData")
genes.OI <- enrich.results %>% 
  filter(FDR <= 0.2 & group_in_pathway > 1) %>% 
  unnest(genes) %>% 
  pull(genes) %>% unique()

#Corr with expressions
E <- read_csv("Methyl/results/eQTM_correlation.csv") %>% 
  mutate(corr.group = case_when(pval < 0.2 ~ R),
         corr.group = case_when(corr.group>0 ~ "+",
                                corr.group<0 ~ "-")) %>% 
  select(condition,gene,DMR, corr.group) %>% 
  pivot_wider(names_from = condition, values_from = corr.group)



DMR.OI <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz") %>% 
  #Filter DMR assoc with signif genes
  mutate(DMR_genes = str_split(DMR_genes, "/")) %>% 
  unnest(DMR_genes) %>% 
  filter(DMR_genes %in% genes.OI) %>% 
  #Fold change direction for probes in region
  mutate(FC = RSTR.M.ave-LTBI.M.ave,
         FC.group = case_when(FC>0 ~ "+",
                              FC<0 ~ "-")) %>% 
  #Collapse per DMR
  group_by(DMR, DMR_genes, no.cpgs) %>% 
  arrange(desc(FC.group),desc(annotation.group)) %>% 
  summarise(annotation.group = paste(unique(annotation.group), collapse="/"),
            FC.groups = paste(unique(FC.group), collapse="/"),
            .groups = "drop") %>% 
  full_join(E) %>% 
  select(DMR_genes, annotation.group, DMR,  no.cpgs,  FC.groups, MEDIA, TB) %>% 
  mutate(DMR = gsub("_"," ",DMR)) %>% 
  arrange(FC.groups, DMR_genes)
DMR.OI


write_csv(DMR.OI, file = "publication/Table2.DMR.direction.csv", na = "n.s.")
