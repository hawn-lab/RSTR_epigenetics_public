library(tidyverse)
library(GO.db)

##### Enrichment #####
load("Methyl/results/enrichment/RSTR.DMPR.enrichment.RData")

# add GO ID
goterms <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename(GO_path = `Term(GOTERM)`)


enrich.results %>% 
  filter(FDR <= 0.2 & `k/K` > 0.04) %>%
  mutate(GO_path = gsub("GOBP_","",pathway),
         GO_path = gsub("_"," ",GO_path),
         GO_path = tolower(GO_path),
         GO_path = gsub("high density","high-density",GO_path),
         GO_path = gsub("protein containing","protein-containing",GO_path),
         GO_path = gsub("protein lipid","protein-lipid",GO_path)) %>%
  left_join(goterms) %>% 
  mutate(genes = as.character(genes)) %>% 
  #Add groups from Fig2
  mutate(GO_group = case_when(grepl("HIPPO", pathway) ~ "Hippo signaling",
                              grepl("FATTY_ACID", pathway) ~ "Fatty acid biosynthesis",
                              grepl("LIPID_EXPORT", pathway) ~ "Lipid export",
                              grepl("HIGH_DENSITY_LIPOPROTEIN|COMPLEX", pathway) ~ "HDL remodeling",
                              TRUE~"Other")) %>% 
  dplyr::select(gs_cat, gs_subcat, pathway, GOID, GO_group,
                group_in_pathway, size_pathway, `k/K`, FDR, genes) %>% 
  write_csv("publication/TableS2.enrichment.csv")
