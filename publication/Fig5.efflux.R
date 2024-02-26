library(tidyverse)
library(patchwork)
library(readxl)
library(ggpubr)

#### Efflux data ####
#Patient metadata
meta <- read_csv("data_refs/2021.09.08RSTR_Hawn_metadata.csv") %>% 
  distinct(FULLIDNO, RS_SUB_ACCESSION_NO, Sample_Group, M0_KCVSEX, KCHCA_AGE_YR_CURRENT,
           avgBMI, RISK_SCORE) %>% 
  drop_na(RS_SUB_ACCESSION_NO) %>% 
  arrange(FULLIDNO)

#Combine with lipoprotein and efflux data
dat <- read_excel("lipoprotein/data_raw/Study191_DMA_Efflux_ForKim211012.xlsx", 
                  sheet="DMA.data") %>% 
  select(AssayID, SUB_ACCESSION_NO) %>% #AssayID key only on sheet 1
  full_join(read_excel("lipoprotein/data_raw/Study191_DMA_Efflux_ForKim211012.xlsx", 
                       sheet="Efflux.data"), by = "AssayID") %>% 
  select(-Group) %>% 
  left_join(meta, by=c("SUB_ACCESSION_NO"="RS_SUB_ACCESSION_NO"))

# Filter to samples with kinship
kin <- read_csv("data_refs/kinship_Hawn_all.csv")

dat.kin <- dat %>% 
  filter(FULLIDNO %in% intersect(kin$rowname, dat$FULLIDNO)) %>% 
  select(AssayID, SUB_ACCESSION_NO, J774.ind:RISK_SCORE) %>% 
  select(-c(J774.basal,J774.abca1,ABCA1.basal,ABCA1.ind)) %>% 
  pivot_longer(J774.ind:ABCA1.spec, names_to = "gene") %>% 
  mutate(facet_lab = fct_recode(gene, 
                                "Total CEC"="J774.ind",
                                "ABCA1 CEC"="ABCA1.spec"
                                ))

#Filter to significant hits
dat.signif <- read_csv("lipoprotein/results/lipoprot_efflux_model_results.csv") %>% 
  filter(variable == "Sample_Group") %>% 
  right_join(dat.kin) %>% 
  rowwise() %>% 
  mutate(facet_lab = paste0(facet_lab,"\nP = ",round(pval,digits=2)),
         facet_lab = factor(facet_lab,
                            levels=c("Total CEC\nP = 0.39",
                                     "ABCA1 CEC\nP = 0.02")))

### Efflux summary ####
# dat.signif %>% 
#   group_by(gene,Sample_Group) %>% 
#   summarise(m = mean(value),
#             sd = sd(value)) %>% 
#   View()
# 
# dat.signif %>% 
#   group_by(gene) %>% 
#   summarise(m = mean(value),
#             sd = sd(value)) %>% 
#   View()

#### Efflux J774 ####
plot1 <- dat.signif %>%
  filter(grepl("J774",gene)) %>% 
  
  ggplot(aes(x=Sample_Group, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Sample_Group), width=0.2, height=0, size=1) +
  theme_classic() +
  labs(x="",y="J774 CEC (%)") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  theme(legend.position = "none") +
  facet_wrap(~facet_lab)
# plot1

#### Efflux BHK ####
plot3 <- dat.signif %>%
  filter(gene %in% c("ABCA1.spec")) %>% 
  
  ggplot(aes(x=Sample_Group, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Sample_Group), width=0.2, height=0, size=1) +
  theme_classic() +
  labs(x="",y=expression(`BHK CEC (`*Delta*`%)`)) +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  theme(legend.position = "none") +
  facet_wrap(~facet_lab, scales="free")
# plot3

#### RNAseq data ####
load("RNAseq/RSTR_RNAseq_data_for_eQTM.RData")
stat.dat <- data.frame(
  group1 = c("Media\nLTBI","Media\nRSTR","+Mtb\nLTBI"),
  group2 = c("+Mtb\nLTBI","+Mtb\nRSTR","+Mtb\nRSTR"),
  y.position  = c(10.6,11.2,11.8),
  symbol = c("1.3E-14","< 1E-16","0.016")
)

#### RNAseq plot ####
plot4 <- as.data.frame(dat.combined.voom$E) %>% 
  rownames_to_column("geneName") %>% 
  filter(geneName == "ABCA1") %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets) %>% 
  mutate(condition = recode(condition, "MEDIA"="Media",
                            "TB"="+Mtb")) %>% 
  mutate(x=paste(condition, Sample_Group, sep="\n"),
         x=factor(x, levels=c("Media\nLTBI", "Media\nRSTR", 
                              "+Mtb\nLTBI", "+Mtb\nRSTR"
         ))) %>% 
  
  ggplot(aes(x=x, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, aes(color=Sample_Group), size=1) +
  stat_pvalue_manual(stat.dat, label="symbol", size=2.5
                     ) +
  theme_classic() +
  labs(x="", y="ABCA1\nlog2 normalized expression") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  theme(legend.position = "none") +
  lims(y=c(6,12))
# plot4

#### Save ####
lo <- "
ABCCC
"
plot_all <- plot1 + plot3 +plot4 + 
  plot_annotation(tag_levels = "A") +
  plot_layout(design = lo)

# plot_all

ggsave(plot_all, file="publication/Fig5.efflux.png",
       width=6.1, height=3)
ggsave(plot_all, file="publication/Fig5.efflux.pdf",
       width=6.1, height=3)
