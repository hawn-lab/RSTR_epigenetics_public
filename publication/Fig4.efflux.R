library(tidyverse)
library(patchwork)
library(readxl)

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
  select(AssayID, SUB_ACCESSION_NO, tHDL:xlHDL.sz) %>% 
  full_join(read_excel("lipoprotein/data_raw/Study191_DMA_Efflux_ForKim211012.xlsx", 
                       sheet="Efflux.data"), by = "AssayID") %>% 
  select(-Group) %>% 
  left_join(meta, by=c("SUB_ACCESSION_NO"="RS_SUB_ACCESSION_NO"))

# Filter to samples with kinship
kin <- read_csv("data_refs/kinship_Hawn_all.csv")

dat.kin <- dat %>% 
  filter(FULLIDNO %in% intersect(kin$rowname, dat$FULLIDNO)) %>% 
  pivot_longer(tHDL:ABCA1.spec, names_to = "gene")

#Filter to significant hits
dat.signif <- read_csv("lipoprotein/results/lipoprot_efflux_model_results.csv") %>% 
  filter(variable == "Sample_Group") %>% 
  left_join(dat.kin)

#### Efflux plot ####

plot1 <- dat.signif %>% 
  filter(gene == "ABCA1.spec" & pval < 0.05) %>% 
  mutate(gene = "ABCA1 specific efflux") %>% 
  #Create facet label
  mutate(facet.lab = paste0(gene, "\nP = ", signif(pval, digits=2))) %>% 
  ggplot(aes(x=Sample_Group, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Sample_Group), width=0.2, height=0) +
  theme_classic() +
  facet_wrap(~facet.lab, scales = "free") +
  labs(x="",y="Cholesterol efflux (%)") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  theme(legend.position = "none")
plot1

#### Save ####
ggsave(plot1, file="publication/Fig4.efflux.png",
       width=2, height=4)
