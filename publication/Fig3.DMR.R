library(tidyverse)
library(patchwork)
library(biomaRt)

#### DMP ####
DMP <- read_csv("Methyl/results/RSTR.DMP.results.csv.gz")

#### DMR ####
#List DMR of interest
genes.OI <- c("APOC3","PLA2G3","KCNQ1", "SHANK2", "CIT")

DMR_all <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz") 
DMR.OI2 <- DMR_all %>% 
  filter( annotation.group != "IGR") %>% 
  mutate(DMR_genes = str_split(DMR_genes, "/")) %>% 
  unnest(DMR_genes)  %>% 
  drop_na(DMR_genes) %>% 
  filter(DMR_genes %in% genes.OI)

#### Gene annotation ####
#Get reference genome
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                   host = "https://useast.ensembl.org")
# searchDatasets(mart = ensembl, pattern = "hsapiens")


gene.POS <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id',
                                 'entrezgene_id',
                                 'chromosome_name',
                                 'strand',
                                 "start_position", "end_position"),
                  mart = ensembl) %>% 
  #Rename
  dplyr::rename(CHR=chromosome_name) %>% 
  # genes of interest
  filter(hgnc_symbol %in% genes.OI) %>% 
  filter(!grepl("^CHR", CHR))

#DMRs in order
DMR.vec <- c("DMR_25","DMR_61","DMR_49","DMR_62","DMR_20","DMR_59")

#### Methyl plots ####
plot.ls <- list()

for(d in DMR.vec) {
  print(d)
  ##### Subset data #####
  #Get DMR info
  DMR_sub <- DMR.OI2 %>% dplyr::filter(DMR==d) %>% 
    mutate(DMR = gsub("_"," ",DMR))
  
  
  DMR_surround <- DMP %>%
    filter(probeID %in% DMR_sub$probeID) %>%
    dplyr::select(probeID, CHR, start_hg38, FDR,
                  RSTR.M.ave, LTBI.M.ave, RSTR.M.sd, LTBI.M.sd) %>%
    #label for DMP
    mutate(FDR=ifelse(FDR<0.2,"*",NA)) %>% 
    mutate(FC = ifelse(RSTR.M.ave > LTBI.M.ave, "RSTR","LTBI")) %>% 
    pivot_longer(RSTR.M.ave:LTBI.M.sd) %>% 
    #keep only 1 label per site
    rowwise() %>% 
    mutate(FDR = ifelse(grepl(FC, name), FDR, NA)) %>% 
    separate(name, into=c("group","measure","name"), sep="[.]") %>% 
    pivot_wider() 
  
  #### plot ####
  plot.title <- paste0(unique(DMR_sub$DMR), " (",
         paste(unique(DMR_sub$annotation.group),
               collapse="/"), ")")
  plot.title <- gsub("body","intron", plot.title)
  
  DMR_M_plot <- 
    DMR_surround %>% 
    
    ggplot(aes(x=start_hg38, y=ave)) +
    geom_ribbon(aes(ymin = ave-sd, ymax = ave+sd, fill = group),
                alpha = 0.3) +
    geom_point(size=1, aes(color = group)) +
    geom_line(aes(color = group)) +
    #signif DMP
    geom_text(aes(label=FDR, y = ave+sd+0.1), size=10, color="black") +
    theme_bw() +
    labs(x="", y="", title = plot.title) +
    scale_color_manual("", labels=c("LTBI.M.ave"="LTBI", 
                                    "RSTR.M.ave"="RSTR"),
                       values = c("#f1a340","#998ec3")) +
    scale_fill_manual("", labels=c("LTBI.M.ave"="LTBI", 
                                   "RSTR.M.ave"="RSTR"),
                      values = c("#f1a340","#998ec3")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # DMR_M_plot
  #Remove legend except last plot
  if(d != DMR.vec[3] &
     d != DMR.vec[length(DMR.vec)]){
    DMR_M_plot <- DMR_M_plot +
      theme(legend.position = "none")
  }
  #Add tag to first plot
  if(d == DMR.vec[1]){
    DMR_M_plot <- DMR_M_plot +
      labs(tag="A", y="Probe log2 M value")
  }
  if(d == DMR.vec[4]){
    DMR_M_plot <- DMR_M_plot +
      labs(tag="C", y="Probe log2 M value")
  }
  ##### Combo #####
  plot.ls[[d]] <- DMR_M_plot
}
# wrap_plots(plot.ls)

#### Gene plots ####
plot.ls2 <- list()

for(d in DMR.vec[-6]) {
  print(d)
  #get second DMR for SHANK2
  if(d == "DMR_20"){
    d <- c(d, "DMR_59")
  }
  #Get DMR info
  DMR_sub <- DMR.OI2 %>% dplyr::filter(DMR %in% d) %>% 
    mutate(DMR = gsub("_"," ",DMR)) %>% 
    distinct(DMR, CHR, DMR_start, DMR_end, DMR_genes)
  
  #Get gene info
  gene_sub <- gene.POS %>% 
    filter(hgnc_symbol %in% DMR_sub$DMR_genes) %>% 
    full_join(DMR_sub, by=c("hgnc_symbol"="DMR_genes"))
  
  #labels
  x_lab <- paste("Chromosome", unique(DMR_sub$CHR), "position")
  
  #set DMR label position
  if(unique(gene_sub$strand) == 1){
    gene_sub_format <- gene_sub %>% 
      filter(CHR.x==CHR.y) %>% 
      mutate(start = start_position,
             end = end_position,
             g_lab = DMR_start,
             g_hjust = -0.2)
  } else{
    gene_sub_format <- gene_sub %>% 
      filter(CHR.x==CHR.y) %>% 
      mutate(end = start_position,
             start = end_position,
             g_lab = DMR_end,
             g_hjust = 1.4)
  }
  
  
  ##### Gene plot #####
  gene_plot <- gene_sub_format %>%
    ggplot() +
    #DMR
    geom_rect(aes(xmin = DMR_start, xmax = DMR_end, ymin = 0.9, ymax = 1.1),
              fill="black", color="black") +
    #Gene
    geom_segment(aes(x=start, xend=end), 
                 y=1, yend=1, color="grey70", size=2,
                 arrow=arrow(length = unit(0.2, "npc"))) +
    geom_text(aes(x=(start+end)/2, label = hgnc_symbol),
              y=0.95, hjust=0.5, vjust=-1) +
    
    theme_classic() +
    lims(y=c(0.9,1.1)) +
    labs(x=x_lab) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.margin = unit(c(0,1,0,0.5), "cm"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_x_continuous(breaks=c(gene_sub_format$start_position, 
                                gene_sub_format$end_position))
  
  #Add tag to first plot
  if(any(d == DMR.vec[1])){
    gene_plot <- gene_plot +
      labs(tag="B")
  }
  if(any(d == DMR.vec[4])){
    gene_plot <- gene_plot +
      labs(tag="D")
  }
  # gene_plot
  
  ##### Save #####
  plot.ls2[[d[1]]] <- gene_plot
}

#### Expression  ####
load("RNAseq/RSTR_RNAseq_data_for_eQTM.RData")

geneOI <- c("APOC3", "PLA2G3", "KCNQ1", #lipid/hdl
            "CIT","SHANK2")             #hippo

# Lmerel
library(kimma)
dat.select <- dat.combined.voom[geneOI,]
kin <- read_csv("data_refs/kinship_Hawn_all.csv")

overlap <- intersect(kin$rowname, dat.select$targets$FULLIDNO)

kin <- kin %>% 
  filter(rowname %in% overlap) %>% 
  dplyr::select(rowname, all_of(overlap)) %>% 
  arrange(rowname) %>% 
  column_to_rownames()

m <- kmFit(dat=dat.select, kin=kin, patientID = "FULLIDNO",
           model = "~condition*Sample_Group+KCHCA_AGE_YR_CURRENT+M0_KCVSEX+experiment+(1|FULLIDNO)",
           run_lmerel = TRUE, 
           run_contrast = TRUE, contrast_var = "condition:Sample_Group")

#No signif interaction or RSTR main
summarize_kmFit(m$lmerel, fdr_cutoff = 0.2)

#Significant Mtb genes
fdr <- m$lmerel %>%
  filter(variable == "condition") %>% 
  filter(FDR < 0.2) %>% 
  mutate(FDR = signif(FDR, digits=2))

dat <- as.data.frame(dat.combined.voom$E) %>% 
  rownames_to_column("geneName") %>% 
  filter(geneName %in% geneOI) %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets, by="libID") %>% 
  mutate(condition = fct_recode(condition, "Media"="MEDIA",
                                "+Mtb"="TB")) %>% 
  mutate(geneName = factor(geneName, levels=geneOI)) %>% 
  #add FDR 
  left_join(fdr, by=c("geneName"="gene")) %>% 
  mutate(facet.lab = paste(geneName, "\nFDR = ", FDR)) %>% 
  mutate(facet.lab2 = case_when(
    geneName %in% c("APOC3", "PLA2G3", "KCNQ1")~"Fatty acids, lipids, HDL",
    geneName %in% c("SHANK2", "CIT")~"Hippo",
    TRUE~"Other"),
    geneName = factor(geneName, levels=c("APOC3", "PLA2G3", "KCNQ1",
                                         "CIT", "SHANK2" ))) %>% 
  arrange(facet.lab2, geneName)

plot3a <- dat %>% 
  filter(geneName %in%c("APOC3", "PLA2G3", "KCNQ1")) %>% 
  ggplot(aes(x=condition, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, aes(color=Sample_Group)) +
  # stat_pvalue_manual(stat.dat, label="symbol") +
  theme_classic() +
  labs(x="", y="Log2 normalized expression",
       tag="E") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  facet_wrap(~facet.lab, ncol=5, scales="free") +
  theme(legend.position = "none")

plot3b <- dat %>% 
  filter(geneName %in%c("CIT","SHANK2")) %>% 
  ggplot(aes(x=condition, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, aes(color=Sample_Group)) +
  # stat_pvalue_manual(stat.dat, label="symbol") +
  theme_classic() +
  labs(x="", y="Log2 normalized expression",
       tag="F") +
  scale_color_manual("", values = c("#f1a340","#998ec3")) +
  facet_wrap(~facet.lab, ncol=5, scales="free")

#### Save ####
layout <- "
AABBCC
DDEEFF
GGHHII
JJKKKK
LLLMMM
"
p_all <- plot.ls[[1]] + plot.ls[[2]] + plot.ls[[3]] +
  plot.ls2[[1]] + plot.ls2[[2]] + plot.ls2[[3]] +
  plot.ls[[4]] + plot.ls[[5]] + plot.ls[[6]] +
  plot.ls2[[4]] + plot.ls2[[5]] +
  plot3a + plot3b +
  plot_layout(heights = c(3,0.7,3,0.7,3), design = layout)
# p_all

ggsave(p_all, filename = "publication/Fig3.DMR.png", width=10, height=10)
ggsave(p_all, filename = "publication/Fig3.DMR.pdf", width=10, height=10)


#### Check DMP in DMR ####
DMP.signif <- DMP %>% 
  filter(FDR<0.2)
DMR_all %>% 
  mutate(DMP = ifelse(probeID %in% DMP.signif$probeID, "Y","N")) %>% 
  group_by(DMR) %>% 
  count(DMP) %>% 
  ungroup() %>% 
  pivot_wider(names_from = DMP, values_from = n) %>% 
  filter(!is.na(Y)) %>% 
  nrow()