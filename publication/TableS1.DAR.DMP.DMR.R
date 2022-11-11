library(tidyverse)

#### DAR ####
DAR <- read_csv("ATACseq/results/RSTR.DAR.results.csv",
                col_types = list(CHR="c")) %>% 
  filter(DAR_FDR <= 0.2) %>% 
  rename(site=peak, start=peak_start, end=peak_end, FDR=DAR_FDR, 
         gene=DAR_genes, LTBI_mean=LTBI.voom.mean, RSTR_mean=RSTR.voom.mean) %>% 
  mutate(group="DAR")

#### DMP ####
DMP <- read_csv("Methyl/results/RSTR.DMP.results.csv.gz",
                col_types = list(seqnames_hg38="c")) %>% 
  filter(FDR <= 0.2) %>% 
  mutate(gene = ifelse(is.na(gene), gene_HGNC_hg38, gene)) %>% 
  select(probeID, FDR, seqnames_hg38, start_hg38, end_hg38, strand_hg38, 
         gene, feature, LTBI.M.ave, RSTR.M.ave) %>% 
  rename(site=probeID, CHR=seqnames_hg38, start=start_hg38,
        end=end_hg38, strand=strand_hg38,
        annotation=feature,
        RSTR_mean=RSTR.M.ave, LTBI_mean=LTBI.M.ave) %>% 
  mutate(annotation.group = ifelse(grepl("TSS",annotation), "TSS1500",
                            ifelse(annotation %in% 
                                     c("Body","ExonBnd","1stExon"), "body",
                                   annotation))) %>%
  #fix no gene annotation group
  mutate(annotation.group = ifelse(is.na(gene), "IGR",
                                   annotation.group)) %>% 
  mutate(CHR=gsub("chr", "",CHR)) %>% 
  mutate(group="DMP") %>% 
  arrange(site)

#### DMR ####
DMR <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz",
                col_types = list(CHR="c")) %>% 
  filter(DMR_FDR <= 0.05) %>% 
  rename(site=DMR, FDR=DMR_FDR, start=DMR_start, end=DMR_end,
         gene=DMR_genes, annotation=feature, 
         RSTR_mean=RSTR.M.ave, LTBI_mean=LTBI.M.ave) %>% 
  distinct(site, probeID, FDR, CHR, start, end, annotation, annotation.group, 
           gene, RSTR_mean, LTBI_mean) %>% 
  #fix no gene annotation group
  mutate(annotation.group = ifelse(is.na(gene), "IGR", annotation.group)) %>% 
  mutate(group="DMR") %>% 
  arrange(DMR, start)

#### Combine and save ####
bind_rows(DAR,DMP,DMR) %>%
  select(group, site, probeID, FDR, CHR:end, 
         gene, annotation.group,
         LTBI_mean, RSTR_mean) %>% 
  arrange(group, FDR) %>% 
  write_csv(file="publication/TableS1.DAR.DMP.DMR.csv")
