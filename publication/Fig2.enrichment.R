library(tidyverse)
library(patchwork)
library(rrvgo)
library(GO.db)
library(gt)
library(gtExtras)

#### Enrichment data ####
attach("Methyl/results/enrichment/RSTR.DMPR.enrichment.RData")

# add GO ID
goterms <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename(GO_path = `Term(GOTERM)`)

enrich.anno <- enrich.results %>% 
  filter(FDR <= 0.2 & `k/K`>0.04) %>%
  mutate(GO_path = gsub("GOBP_","",pathway),
         GO_path = gsub("_"," ",GO_path),
         GO_path = tolower(GO_path),
         GO_path = gsub("high density","high-density",GO_path),
         GO_path = gsub("protein containing","protein-containing",GO_path),
         GO_path = gsub("protein lipid","protein-lipid",GO_path),
         GO_path = recode(GO_path, 
                          "protein-lipid complex subunit organization"=
                            "protein-lipid complex organization")) %>%
  left_join(goterms)

#### Bubble plot ####
# Group like terms based on gene content.
simMatrix <- calculateSimMatrix(enrich.anno$GOID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(enrich.anno$FDR), enrich.anno$GOID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.5,
                                orgdb="org.Hs.eg.db")

# Adapted from rrvgo::scatterPlot to label all points
set.seed(42)
x <- cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)

df <- cbind(as.data.frame(x$points),
            reducedTerms[match(rownames(x$points), reducedTerms$go),
                         c("term", "parent", "parentTerm", "size")]) %>% 
  left_join(enrich.anno, by=c("term"="GO_path")) %>% 
  #recode small groups
  mutate(GO_group = recode(parentTerm,
                           "protein-lipid complex organization"="HDL remodeling",
                           "biological process involved in intraspecies interaction between organisms"="Other",
                           "negative regulation of hippo signaling"="Hippo signaling",
                           "high-density lipoprotein particle remodeling"="HDL remodeling",
                           "lipid export from cell"="Lipid export",
                           "rhythmic behavior"="Other",
                           "cochlea development"="Other",
                           "regulation of fatty acid biosynthetic process"="Fatty acid\nbiosynthesis"))

col.vec <- c("Fatty acid\nbiosynthesis"="#CEB125",
             "HDL remodeling"="#AA4499",
             "Hippo signaling"="#1EA1E2",
             "Lipid export"="#117733", 
             "Other"="grey60",
             "none"="grey90")

df_lab <- df %>% 
  group_by(GO_group) %>% 
  mutate(V1=mean(V1),V2=mean(V2),
         V1 = case_when(GO_group == "Hippo signaling" ~ V1+0.3,
                        GO_group == "HDL remodeling" ~ V1-0.3,
                        GO_group == "Lipid export" ~ V1-0.15,
                        GO_group == "Other" ~ V1-0.15,
                        TRUE ~ V1),
         V2 = case_when(GO_group == "Other" ~ V2+0.06,
                        GO_group == "Fatty acid\nbiosynthesis" ~ V2-0.08,
                        GO_group == "Lipid export" ~ V2-0.06,
                        GO_group == "HDL remodeling" ~ V2-0.04,
                        TRUE ~ V2))

p1 <- df %>% 
  arrange(desc(GO_group)) %>% 
  ggplot(aes(x=V1, y=V2, color=GO_group)) +
  geom_point(aes(size=`k/K`*100), alpha=.5) +
  scale_color_discrete(guide="none") +
  scale_size_continuous(range=c(1, 10), breaks = c(5,10,15)) +
  labs(size="Percent\nenrichment", title = "A") +
  theme_bw() +
  theme(axis.text   = element_blank(),
        axis.title  = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.83, 0.635),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.box.background = element_rect(colour = "black", size=1.5)
  ) + 
  geom_text(data=df_lab, aes(label=GO_group)) +
  scale_color_manual(values=col.vec ,
                     guide="none") +
  lims(x=c(-0.53,0.53),
       y=c(-0.35,0.43)) +
  coord_fixed()

# p1

#### Save ####
ggsave(p1, filename = "publication/Fig2A.enrichment.png",
       width=3.1, height=2.5)
ggsave(p1, filename = "publication/Fig2A.enrichment.pdf",
       width=3.1, height=2.5)


#### "heatmap" of terms ####
genes.OI <- df %>%
  unnest(genes) %>% 
  pull(genes) %>% unique()

enrich.hm <- df %>% 
  unnest(genes) %>% 
  distinct(genes, GO_group) %>% 
  mutate(GO_group = recode(GO_group,
                           "Hippo signaling"="Hippo",
                           "HDL remodeling"="HDL", 
                           "Fatty acid\nbiosynthesis"="FA", 
                           "Lipid export"="Lipid")) %>% 
  arrange(GO_group) %>% 
  #Make levels for each pathway group
  mutate(value = as.numeric(factor(GO_group))) %>% 
  #Make levels as space blanks
  rowwise() %>% 
  mutate(value = paste0(rep(" ", value), collapse="")) %>% 
  #wide format
  pivot_wider(names_from = GO_group) %>% 
  #make pretty
  replace(is.na(.), "") %>% 
  arrange(genes) %>% 
  dplyr::rename(Gene=genes)

blank.colnames <- c(" ","  ","   ","    ","     ")
# colnames(enrich.hm) <- c("Gene", blank.colnames)

#### DMR annotation and FC ####
DMR <- read_csv("Methyl/results/RSTR.DMR.results.csv")

DMR.OI <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz") %>% 
  #Filter DMR assoc with signif genes
  mutate(DMR_genes = str_split(DMR_genes, "/")) %>% 
  unnest(DMR_genes) %>% 
  # filter(DMR_genes %in% genes.OI) %>% 
  filter(annotation.group != "IGR") %>% 
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
  mutate(annotation.group = recode(annotation.group,
                                   "body/TSS1500"="TSS1500/body",
                                   "exon/5'UTR"="5'UTR/exon",
                                   "body/5'UTR"="5'UTR/body",
                                   "exon/body"="exon",
                                   "5'UTR/TSS1500"="TSS1500/5'UTR")) %>% 
  #add ave DMR info
  left_join(distinct(DMR, DMR, meandiff)) %>% 
  mutate(meandiff = case_when(meandiff>0 ~ "(+)",
                              meandiff<0 ~ "(-)")) %>% 
  rowwise() %>% 
  mutate(FC.groups = ifelse(FC.groups == "-/+",paste(FC.groups,meandiff),
                            FC.groups)) %>% 
  # full_join(E) %>% 
  #Make pretty
  dplyr::select(DMR_genes, annotation.group, DMR,  no.cpgs,  FC.groups #, MEDIA, TB
  ) %>% 
  mutate(annotation.group = gsub("body","intron",annotation.group)) %>% 
  mutate(DMR = gsub("_"," ",DMR)) %>% 
  arrange(FC.groups, DMR_genes) %>% 
  dplyr::rename(Gene=DMR_genes, Annotation=annotation.group,
                "Total probes"=no.cpgs, "Log2 M fold change RSTR - LTBI"=FC.groups) %>% 
  drop_na(Gene)
# DMR.OI

#### Render in dt ####

dt_all <- inner_join(enrich.hm, DMR.OI) %>% 
  arrange(`Log2 M fold change RSTR - LTBI`, Gene) %>% 
  gt() %>% 
  #Color boxes
  data_color(
    columns = FA:Other,
    colors = scales::col_factor(
      palette = c("white","#CEB125","#AA4499","#1EA1E2","#117733","grey60"),
      domain = c("",blank.colnames)),
    apply_to = "fill",
    autocolor_text = FALSE) %>% 
  #Color text for + / - fold change
  data_color(
    columns = `Log2 M fold change RSTR - LTBI`,
    colors = scales::col_factor(
      palette = c("blue","red","blue","red"),
      domain = c("-","+","+/- (-)","+/- (+)")),
    apply_to = "text",
    autocolor_text = FALSE) %>% 
  #Column width
  cols_width(FA:Other ~ px(50)) %>% 
  cols_width(`Log2 M fold change RSTR - LTBI` ~ px(160)) %>% 
  #Text align
  cols_align(align = "left",columns = Gene:DMR) %>% 
  cols_align(align = "center",
             columns = `Total probes`:`Log2 M fold change RSTR - LTBI`) %>% 
  #row breaks
  tab_options(column_labels.border.top.color = "black",
              column_labels.border.bottom.color = "black",
              table_body.border.bottom.color = "black") %>% 
  tab_style(
    style = cell_borders(sides = c("top", "bottom"),color = "white"),
    locations = cells_body(columns = everything(),rows = everything())) %>% 
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(2)),
    locations = cells_body(columns = everything(),rows = c(2,6))) %>% 
  tab_style(
    style = cell_borders(sides = "top"),
    locations = cells_body(columns = everything(),rows = 1)) %>% 
  tab_header("B") %>% 
  opt_align_table_header(align = "left")
# dt_all

#### save ####
gtsave(dt_all, "publication/Fig2B.DMR.png")

DMR.OI %>% 
  separate(`Log2 M fold change RSTR - LTBI`, 
           into=c("Probe log2 M fold change RSTR - LTBI",
                  "Mean DMR log2 M fold change RSTR - LTBI"), sep=" ") %>% 
  mutate(`Mean DMR log2 M fold change RSTR - LTBI` = ifelse(
    is.na(`Mean DMR log2 M fold change RSTR - LTBI`), 
    paste0("(",`Probe log2 M fold change RSTR - LTBI`,")"),
    `Mean DMR log2 M fold change RSTR - LTBI`)) %>% 
  arrange(Gene) %>% 
  write_csv(file = "publication/TableS3.DMR.direction.csv")
