library(tidyverse)
library(BIGpicture)
library(patchwork)

#### GSEA - MTB ####
load("Methyl/results/enrichment/MTB.RSTR.gsea.RData")

gsea_format <- gsea_MTB_c5_signif  %>% 
  mutate(signif.col = case_when(pval < 0.05 ~ "P < 0.05",
                                TRUE  ~ "NS")) %>% 
  mutate(signif.col = fct_relevel(signif.col, "NS", after=Inf)) %>% 
  mutate(group = fct_recode(group,
                            "+Mtb vs media in LTBI"="TB LTBI - MEDIA LTBI",
                            "+Mtb vs media in RSTR"="TB RSTR - MEDIA RSTR",
                            "RSTR vs LTBI in media"="MEDIA RSTR - MEDIA LTBI",
                            "RSTR vs LTBI in +Mtb"="TB RSTR - TB LTBI")) %>% 
  mutate(facet_lab = "Monocyte\nMtb:RSTR")

#### GSEA - MDM ####
load("Methyl/results/enrichment/MDM.IFN.gsea.RData")

gsea_format2 <- gsea_IFN_c5_signif  %>% 
  mutate(signif.col = case_when(pval < 0.05 ~ "P < 0.05",
                                TRUE  ~ "NS")) %>% 
  mutate(signif.col = fct_relevel(signif.col, "NS", after=Inf)) %>% 
  filter(group != "IFNL1") %>% 
  mutate(group = fct_recode(group,
                            "IFNa8"="IFNA8")) %>% 
  mutate(facet_lab = "MDM IFN\nstimulation")

##### GSEA plot ####
#Combine data
gsea_format_all <- bind_rows(gsea_format,gsea_format2) %>% 
  mutate(pathway = gsub("GOBP_","",pathway),
         pathway = gsub("_"," ",pathway),
         pathway = tolower(pathway),
         pathway = gsub("high density","high-density",pathway),
         pathway = gsub("protein containing","protein-containing",pathway),
         pathway = gsub("protein lipid","protein-lipid",pathway),
         pathway = gsub("hippo","Hippo",pathway)) %>% 
  mutate(group = factor(group, levels = c("RSTR vs LTBI in +Mtb",
                                          "RSTR vs LTBI in media",
                                          "+Mtb vs media in RSTR",
                                          "+Mtb vs media in LTBI",
                                          "IFNg","IFNe",
                                          "IFNb","IFNa8")),
         facet_lab = factor(facet_lab, levels=c("Monocyte\nMtb:RSTR",
                                                "MDM IFN\nstimulation")))

#significant pw
gsea.signif <- gsea_format_all %>% 
  group_by(pathway) %>% 
  summarise(minP = min(pval)) %>% 
  filter(minP < 0.05) %>% 
  pull(pathway)

p1 <- gsea_format_all %>% 
  filter(pathway %in% gsea.signif) %>% 
  mutate(pathway = factor(pathway, levels = rev(gsea.signif))) %>% 
  ggplot() +
  aes(x=group,y=NES, fill=signif.col) +
  geom_segment(aes(group, xend = group, y = 0, yend = NES)) +
  geom_point(shape=21, size=3, stroke=1) +
  coord_flip() +
  scale_fill_manual(values = c("P < 0.05"="darkred",
                               "NS"="grey80"),
                    drop = FALSE) +
  geom_hline(yintercept = 0)+
  theme_classic() +
  labs(fill="Significance", x="",
       y="Normalized enrichment score") +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  facet_grid(facet_lab~pathway, scales = "free_y")
# p1

#### STRING ####
attach("RNAseq/MDM-IFN_voom.RData")
gsea_format_hgnc <- gsea_format %>% 
  filter(pathway %in% gsea.signif & pval < 0.05) %>% 
  unnest(leadingEdge) %>% 
  mutate(group = recode(group,
                        "+Mtb vs media in LTBI"="Mtb",
                        "+Mtb vs media in RSTR"="Mtb",
                        "RSTR v LTBI in +Mtb"="RSTR")) %>% 
  group_by(group,gs_cat, gs_subcat, pathway) %>% 
  summarise(leadingEdge = list(unique(leadingEdge))) %>% 
  ungroup() %>% 
  #new group to color by IFN instead of pathway
  mutate(group3=pathway, pathway=fct_rev(group)) %>% 
  arrange(pathway)

gsea_format2_hgnc <- gsea_format2 %>% 
  filter(pathway %in% gsea.signif2) %>% 
  unnest(leadingEdge) %>% 
  mutate(group = recode(group,
                        "IFNG"="IFN other",
                        "IFNB"="IFN other",
                        "IFNE"="IFN other")) %>% 
  left_join(dat.voom$genes, by=c("leadingEdge"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  group_by(group,gs_cat, gs_subcat, pathway) %>% 
  summarise(leadingEdge = list(symbol)) %>% 
  ungroup() %>% 
  #new group to color by IFN instead of pathway
  mutate(group2=pathway, pathway=fct_rev(group)) %>% 
  arrange(pathway)

gsea_format3_hgnc <- gsea_format2 %>% 
  filter(pathway %in% gsea.signif2) %>% 
  unnest(leadingEdge) %>% 
  mutate(group = recode(group,
                        "IFNA8"="IFN other",
                        "IFNB"="IFN other",
                        "IFNE"="IFN other")) %>% 
  left_join(dat.voom$genes, by=c("leadingEdge"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  group_by(group,gs_cat, gs_subcat, pathway) %>% 
  summarise(leadingEdge = list(symbol)) %>% 
  ungroup() %>% 
  #new group to color by IFN instead of pathway
  mutate(group2=pathway, pathway=fct_rev(group)) %>% 
  arrange(pathway)

#### IFNA lipids plot ####
#all genes in gene set
library(msigdbr)
col.vec <- c("#117733","#44AA99", "#CC6677","#AA4499",
             "#88CCEE","#DDCC77","grey80")

go.db <- msigdbr("human", "C5","GO:BP")

# All genes in pw
genesA <- go.db %>% 
  filter(gs_name == "GOBP_LIPID_EXPORT_FROM_CELL") %>% 
  pull(human_gene_symbol) %>% unique()

gsea_format_hgncA <- bind_rows(gsea_format_hgnc,gsea_format2_hgnc) %>% 
  add_row(group="DMR", pathway = "DMR gene", 
          leadingEdge=list("PLA2G3","KCNQ1")) %>% 
  mutate(NES=1, pval=0, FDR=0) %>% 
  mutate(group = factor(group, levels = c("Mtb","RSTR","DMR gene", 
                                          "IFNA8","IFN other",
                                          "none")))

mapA <- map_string(genesA, score_threshold = 400)

plotA <- plot_string(mapA, enrichment = gsea_format_hgncA, 
                     fdr_cutoff = 1, layout="lgl",
                     node_size = 1.7, text_size = 3,
                     edge_min = 1,
                     colors = col.vec[-4]) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.border = element_rect(fill=NA))

plotA2 <- plot_string(mapA, enrichment = gsea_format_hgncA,
                      fdr_cutoff = 1, layout="grid",
                      node_size = 1.7, text_size = 3,
                      edge_max =0,
                      colors = col.vec[-c(4,6)]) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.border = element_rect(fill=NA))

# plotA + plotA2

#### IFNG Hippo ####
genesG <- go.db %>% 
  filter(gs_name == "GOBP_HIPPO_SIGNALING") %>% 
  pull(human_gene_symbol) %>% unique()

# genesG <- gsea_format_hgnc %>% 
#   filter(group == "IFNG" & pval < 0.05) %>% 
#   unnest(leadingEdge) %>% 
#   pull(leadingEdge) %>% unique()
# 
# #add methyl genes
# genesG <- c(genesG, "SHANK2","CIT")

gsea_format_hgncG <- bind_rows(gsea_format_hgnc,gsea_format3_hgnc) %>% 
  add_row(group="DMR", pathway = "DMR gene", 
          leadingEdge=list("SHANK2","CIT")) %>% 
  mutate(NES=1, pval=0, FDR=0) %>% 
  mutate(group = factor(group, levels = c("Mtb","RSTR","DMR gene", 
                                          "IFNG","IFN other",
                                          "none")))

mapG <- map_string(genesG, score_threshold = 400)

plotG <- plot_string(mapG, enrichment = gsea_format_hgncG, 
                     fdr_cutoff = 1, 
                     layout = "kk", node_size = 1.7, text_size = 3,
                     edge_min = 1,
                     colors=col.vec[c(2,5,4,6,7)])+
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.border = element_rect(fill=NA))

# plotG

plotG2 <- plot_string(mapG, enrichment = gsea_format_hgncG, 
                      fdr_cutoff = 1, 
                      layout = "grid", node_size = 1.7, text_size = 3,
                      edge_max = 0,
                      colors=col.vec[c(5,4,7)])+
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.border = element_rect(fill=NA))

# plotG + plotG2

#### Save ####
ggsave("publication/Fig4.RNAseq.GSEA.AB.pdf",
       p1, width=8, height=2.5)

p2 <- plotA + plotA2 + plotG + plotG2+
  plot_layout() +
  plot_annotation(tag_levels = "A")

ggsave("publication/Fig4.RNAseq.GSEA.C.pdf",
       p2, width=8, height=18)
