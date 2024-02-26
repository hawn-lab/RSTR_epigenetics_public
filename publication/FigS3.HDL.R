library(tidyverse)
library(readxl)
library(patchwork)

#### Data ####
meta <- read_csv("data_refs/2021.09.08RSTR_Hawn_metadata.csv") %>%
  distinct(FULLIDNO, RS_SUB_ACCESSION_NO, Sample_Group) %>%
  drop_na() %>%
  rename(SUB_ACCESSION_NO=RS_SUB_ACCESSION_NO)

dat1 <- read_excel("lipoprotein/data_raw/Study191_DMA_Efflux_ForKim211012.xlsx", 
                   sheet="DMA.data") %>% 
  left_join(meta) %>% 
  select(AssayID, FULLIDNO, Sample_Group, tHDL:xlHDL, -ssHDL, -msuHDL)

dat2 <- read_excel("lipoprotein/data_raw/Study191_DMA_Efflux_ForKim211012.xlsx", 
                   sheet="DMA.data") %>% 
  select(AssayID, xsHDL.sz:xlHDL.sz)

fdr <- read_csv("lipoprotein/results/lipoprot_efflux_model_results.csv") %>% 
  select(gene, variable, estimate, pval, FDR.new) %>% 
  filter(gene != "J774.abca1") %>% 
  filter(variable=="Sample_Group") %>% 
  mutate(group = case_when(grepl("sz", gene)~"HDL_size",
                           TRUE~"HDL")) %>% 
  mutate(signif = ifelse(grepl("HDL",gene), FDR.new, pval)) %>% 
  select(group, gene, signif)

#### Combine data ####
dat <- full_join(dat1, dat2) %>% 
  select(-AssayID) %>% 
  pivot_longer(-c(FULLIDNO:Sample_Group), names_to = "gene") %>%
  inner_join(fdr) %>% 
  mutate(signif_lab = signif(signif, digits=2)) %>% 
  mutate(facet_lab = case_when(group=="HDL_efflux"~paste0(gene,"\nP = ", signif_lab),
                               TRUE ~ paste0(gene,"\nFDR = ", signif_lab))) %>% 
  mutate(facet_lab = gsub(".sz"," size",facet_lab),
         facet_lab = gsub("_"," ",facet_lab)) %>% 
  # order size
  mutate(hdl_ord = factor(gene, levels=c("xsHDL","xsHDL.sz",
                                         "msuHDL",
                                         "ssHDL","sHDL","sHDL.sz",
                                         "mHDL","mHDL.sz",
                                         "mlHDL","mlHDL.sz",
                                         "lHDL","lHDL.sz",
                                         "xlHDL","xlHDL.sz",
                                         "tHDL")),
         hdl_ord = as.numeric(hdl_ord),
         hdl_ord = ifelse(is.na(hdl_ord), 16, hdl_ord),
         facet_lab_ord = fct_reorder(facet_lab, hdl_ord))


          
#### HDL quant ####
p1 <- dat %>% 
  filter(grepl("HDL",gene) & !grepl("sz",gene)) %>% 
  ggplot(aes(x=Sample_Group, y=value))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Sample_Group), width=0.2, height=0, size=1) +
  theme_classic() +
  labs(x="",y="HDL (uM)",
       title = "HDL concentration") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  theme(legend.position = "none") +
  facet_wrap(~facet_lab_ord, scales = "free_y", ncol=4)
# p1

#### HDL size ####
p2 <- dat %>% 
  filter(grepl("sz",gene))  %>% 
  ggplot(aes(x=Sample_Group, y=value))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=Sample_Group), width=0.2, height=0, size=1) +
  theme_classic() +
  labs(x="",y="Particle diameter (nm)",
       title = "HDL size") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  theme(legend.position = "none") +
  facet_wrap(~facet_lab_ord, scales = "free_y", ncol=3)
# p2

### Save ####
p_all <- p1+p2 + plot_layout(widths = c(4,3.5), nrow = 1) +
  plot_annotation(tag_levels = "A")
# p_all

ggsave(filename = "publication/FigureS3.HDL.png", p_all, height=5, width=8)
ggsave(filename = "publication/FigureS3.HDL.pdf", p_all, height=5, width=8)
