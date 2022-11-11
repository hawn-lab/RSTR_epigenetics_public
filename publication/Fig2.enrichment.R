library(tidyverse)
library(patchwork)
library(rrvgo)
library(GO.db)
library(BIGpicture)

#### Enrichment ####
attach("Methyl/results/enrichment/RSTR.DMPR.enrichment.RData")

# add GO ID
goterms <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename(GO_path = `Term(GOTERM)`)

enrich.anno <- enrich.results %>% 
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
                                threshold=0.7,
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
                           "biological process involved in intraspecies interaction between organisms"="Other",
                           "negative regulation of hippo signaling"="Hippo\nsignaling",
                           "regulation of cellular component size"="Cellular\ncomponent size",
                           "high-density lipoprotein particle remodeling"="HDL\nremodeling",
                           "lipid export from cell"="Lipid export",
                           "rhythmic behavior"="Other",
                           "cochlea development"="Other",
                           "regulation of fatty acid biosynthetic process"="Fatty acid\nbiosynthesis",
                           "positive regulation of cell projection organization"="Cell projection\norganization"))

col.vec <- c("Cell projection\norganization"="#5D3A9B",
             "Cellular\ncomponent size"="#009E73",
             "Fatty acid\nbiosynthesis"="#56B4E9",
             "HDL\nremodeling"="#E69F00",
             "Hippo\nsignaling"="#CC79A7",
             "Lipid export"="#D55E00", 
             "Other"="grey60",
             "none"="grey90")

df_lab <- df %>% 
  group_by(GO_group) %>% 
  mutate(V1=mean(V1),V2=mean(V2),
         V1 = case_when(GO_group == "Hippo\nsignaling" ~ V1-0.07,
                        GO_group == "HDL\nremodeling" ~ V1+0.03,
                        TRUE ~ V1),
         V2 = ifelse(!GO_group %in% c("Other"),
                     V2+0.07, V2-0.07)) 

p1 <- df %>% 
  arrange(desc(GO_group)) %>% 
  ggplot(aes(x=V1, y=V2, color=GO_group)) +
  geom_point(aes(size=`k/K`*100), alpha=.5) +
  scale_color_discrete(guide="none") +
  scale_size_continuous(range=c(1, 15)) +
  labs(size="Percent enrichment") +
  theme_bw() +
  theme(axis.text   = element_blank(),
        axis.title  = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.82),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        legend.box.background = element_rect(colour = "black", size=1.5)) + 
  geom_text(data=df_lab, aes(label=GO_group)) +
  scale_color_manual(values=col.vec ,
                     guide="none") 
p1

#### STRING plot ####
fake.enrichment <- df %>% 
  dplyr::select(GO_group, genes) %>% 
  unnest(genes) %>% 
  distinct() %>% 
  group_by(GO_group) %>% 
  summarise(genes=list(genes)) %>% 
  dplyr::rename(pathway = GO_group) %>% 
  mutate(FDR = 0.1, group_in_pathway = 2, `k/K` = 1) %>% 
  mutate(pathway= sub("\n"," ", pathway))

col.vec2 <- c("Cell projection organization"="#AE9CCD",
              "Cellular component size"="#7FCEB9",
              "Fatty acid biosynthesis"="#AAD9F4",
              "HDL remodeling"="#F2CF7F",
              "Hippo signaling"="#E5BCD3",
              "Lipid export"="#E7A571", 
              "Other"="#B3B3B3",
              "none"="grey90")
#All DMR genes
DMR.signif <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz") %>% 
  filter( annotation.group != "IGR") %>% 
  distinct(DMR_genes) %>% 
  mutate(DMR_genes = str_split(DMR_genes, "/")) %>% 
  unnest(DMR_genes)  %>% 
  drop_na(DMR_genes) %>% 
  pull(DMR_genes) %>% unique()

map <- map_string(genes = DMR.signif)

p2 <- plot_string(map, enrichment = fake.enrichment, colors = col.vec2, 
                  text_size = 3, node_size = 1.5)
# p2

#### Save ####
p_all <- p1 + p2 +
  plot_annotation(tag_levels="A")
# p_all

ggsave(p_all, filename = "publication/Fig2.enrichment.pdf",
       width=15, height=6)
ggsave(p_all, filename = "publication/Fig2.enrichment.png",
       width=15, height=6)
