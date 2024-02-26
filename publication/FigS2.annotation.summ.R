library(tidyverse)

#### DAR ####
DAR <- read_csv("ATACseq/results/RSTR.DAR.results.csv",
                col_types = list(CHR="c")) %>% 
  filter(DAR_FDR <= 0.2) %>% 
  rename(site=peak) %>% 
  distinct(site, annotation.group) %>% 
  mutate(group="DAR")

#### DMP ####
DMP <- read_csv("Methyl/results/RSTR.DMP.results.csv.gz",
                col_types = list(seqnames_hg38="c")) %>% 
  filter(FDR <= 0.2) %>% 
  rename(site=probeID) %>% 
  mutate(annotation.group = ifelse(grepl("TSS",feature), "TSS1500",
                                   ifelse(feature == "Body", "intron",
                                          ifelse(feature %in% 
                                                   c("ExonBnd","1stExon"), "exon",
                                                 feature)))) %>%
  distinct(site, annotation.group) %>% 
  mutate(group="DMP")

#### DMR ####
DMR <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz",
                col_types = list(CHR="c")) %>% 
  rename(site=DMR) %>% 
  #Fill annotation if IGR if gene is NA
  mutate(annotation.group = ifelse(is.na(DMR_genes), "IGR", annotation.group)) %>% 
  distinct(site, annotation.group) %>% 
  mutate(group="DMR") %>% 
  #collapse multi anno
  mutate(annotation.group = gsub("body","intron",annotation.group)) %>% 
  group_by(site,group) %>% 
  summarise(annotation.group = paste(sort(unique(annotation.group), decreasing=TRUE), collapse=" & "),
            .groups="drop") %>% 
  mutate(annotation.group = case_when(grepl("TSS",annotation.group)~"TSS1500",
                                      grepl("exon",annotation.group)~"exon",
                                      grepl("5'UTR",annotation.group)~"5'UTR",
                                      grepl("3'UTR",annotation.group)~"3'UTR",
                                      grepl("body",annotation.group)~"intron",
                                      TRUE~annotation.group))

#### Combine ####
dat <- bind_rows(DAR,DMP,DMR) %>%
  count(group, annotation.group) %>% 
  mutate(annotation.group = factor(annotation.group, 
                                   levels=c("TSS1500","exon","5'UTR","3'UTR",
                                            "intron","IGR")))
#total sites
dat.total <- dat %>% 
  group_by(group) %>% 
  summarise(total=sum(n), .groups = "drop")

#### Plot #### 
plot <- dat %>% 
  ggplot() +
  geom_bar(aes(x=group, y=n, fill=annotation.group),
           position = "fill",stat = "identity") + 
  theme_classic() +
  labs(x="",y="Proportion of significant sites",fill="Annotation") +
  geom_text(data=dat.total, aes(label=total, x=group, y=1.01), 
            position=position_dodge(width=0.9), vjust=-0.25)


ggsave(filename="publication/FigS2.annotation.summ.png", plot,
       height=4, width=3)
