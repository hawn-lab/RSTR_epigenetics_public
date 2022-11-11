library(tidyverse)
library(GO.db)

##### Enrichment #####
load("Methyl/results/enrichment/RSTR.DMPR.enrichment.RData")

# add GO ID
goterms <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename(GO_path = `Term(GOTERM)`)


enrich.results %>% 
  filter(FDR <= 0.2 & group_in_pathway > 1) %>%
  mutate(GO_path = gsub("GOBP_","",pathway),
         GO_path = gsub("_"," ",GO_path),
         GO_path = tolower(GO_path),
         GO_path = gsub("high density","high-density",GO_path),
         GO_path = gsub("protein containing","protein-containing",GO_path),
         GO_path = gsub("protein lipid","protein-lipid",GO_path),
         GO_path = gsub("anterior posterior","anterior/posterior",GO_path),
         GO_path = gsub("wnt","Wnt",GO_path),
         GO_path = gsub("cell cell","cell-cell",GO_path),
         GO_path = gsub("hormone mediated","hormone-mediated",GO_path),
         GO_path = gsub("g protein coupled","G protein-coupled",GO_path),
         GO_path = gsub("adenylate cyclase activating","adenylate cyclase-activating", GO_path),
         GO_path = gsub("adenylate cyclase modulating","adenylate cyclase-modulating", GO_path),
         GO_path = recode(GO_path,
                          "cell-cell signaling by Wnt"="cell-cell signaling by wnt")) %>%
  left_join(goterms) %>% 
  dplyr::select(gs_cat, gs_subcat, pathway, GOID,
         group_in_pathway, size_pathway, FDR, `k/K`, genes) %>% 
  mutate(genes = as.character(genes)) %>% 
  
  write_csv("publication/TableS2.enrichment.csv")
