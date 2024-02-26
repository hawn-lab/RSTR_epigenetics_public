library(tidyverse)

#### DAR ####
DAR <- read_csv("ATACseq/results/RSTR.DAR.results.csv",
                col_types = list(CHR="c")) %>% 
  filter(DAR_FDR <= 0.2) %>% 
  rename(site=peak, start=peak_start, end=peak_end, FDR=DAR_FDR, 
         gene=DAR_genes, LTBI_mean=LTBI.voom.mean, RSTR_mean=RSTR.voom.mean) %>% 
  mutate(group="DAR") %>%
  select(group, site, estimate, FDR, CHR:end, 
         gene, annotation.group,
         LTBI_mean, RSTR_mean)

#### DMP ####
DMP_raw <- read_csv("Methyl/results/RSTR.DMP.results.csv.gz",
                    col_types = list(seqnames_hg38="c")) 
DMP <- DMP_raw %>% 
  filter(FDR <= 0.2) %>% 
  mutate(gene = ifelse(is.na(gene), gene_HGNC_hg38, gene)) %>% 
  select(probeID, estimate, FDR, seqnames_hg38, start_hg38, end_hg38, strand_hg38, 
         gene, feature, LTBI.M.ave, RSTR.M.ave) %>% 
  rename(site=probeID, CHR=seqnames_hg38, start=start_hg38,
         end=end_hg38, strand=strand_hg38,
         annotation=feature,
         RSTR_mean=RSTR.M.ave, LTBI_mean=LTBI.M.ave) %>% 
  mutate(annotation.group = case_when(grepl("TSS",annotation)~"TSS1500",
                                      annotation %in% 
                                        c("ExonBnd","1stExon")~"exon",
                                      annotation %in% c("Body","body") ~ "intron",
                                      TRUE~annotation)) %>%
  #fix no gene annotation group
  mutate(annotation.group = ifelse(is.na(gene), "IGR",
                                   annotation.group)) %>% 
  mutate(CHR=gsub("chr", "",CHR)) %>% 
  mutate(group="DMP") %>%
  select(group, site, estimate, FDR, CHR:end, 
         gene, annotation.group,
         LTBI_mean, RSTR_mean) %>% 
  arrange(site)

#### DMR ####
DMR <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz",
                col_types = list(CHR="c")) %>% 
  rename(site=DMR, FDR=DMR_FDR, start=DMR_start, end=DMR_end,
         gene=DMR_genes, 
         RSTR_mean=RSTR.M.ave, LTBI_mean=LTBI.M.ave,
         total_probes=no.cpgs) %>% 
  mutate(group="DMR") %>% 
  group_by(group, site, total_probes, FDR, 
           CHR, start, end, gene) %>% 
  summarise(annotation.group = paste(unique(annotation.group), 
                                     collapse="/"),
            .groups = "drop") %>% 
  #add mean diff
  left_join(
    read_csv("Methyl/results/RSTR.DMR.results.csv") %>% 
      select(DMR, meandiff) %>% 
      rename(site=DMR)) %>% 
  rename(estimate=meandiff) %>% 
  select(group:total_probes, estimate, everything())

#### DMR probes ####
DMR.probe <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz",
                      col_types = list(CHR="c")) %>% 
  rename(site=DMR, start=probe_start, gene=DMR_genes, 
         RSTR_mean=RSTR.M.ave, LTBI_mean=LTBI.M.ave)  %>% 
  mutate(end=start+1)

DMP.in.DMR <- DMP_raw %>% 
  filter(probeID %in% DMR.probe$probeID) %>% 
  select(probeID, estimate, FDR)

DMR.probe <- DMR.probe %>% 
  #add probe FDR
  left_join(DMP.in.DMR) %>% 
  mutate(group="DMR_probes") %>% 
  select(group, site, probeID, estimate, FDR,
         CHR, start,end, LTBI_mean, RSTR_mean)

#### add min/max to DMR ####
DMR.summ <- DMR.probe %>% 
  group_by(site) %>% 
  summarise(FC_min = min(estimate),
            FC_max = max(estimate),
            .groups = "drop")

DMR2 <- DMR %>% 
  left_join(DMR.summ) %>% 
  mutate(annotation.group = gsub("body","intron",annotation.group))

#### Combine and save ####
library(openxlsx)

dfs <- list("DAR"=DAR, "DMP"=DMP, "DMR"=DMR2,
            "DMR_probe"=DMR.probe)
write.xlsx(dfs, file = "publication/TableS1.DAR.DMP.DMR.xlsx")
