library(tidyverse)
library(patchwork)
library(biomaRt)

#### DMP ####
DMP <- read_csv("Methyl/results/RSTR.DMP.results.csv.gz")

#### DMR ####
#List DMR of interest
genes.OI <- c("APOC3","PLA2G3","KCNQ1")

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
searchDatasets(mart = ensembl, pattern = "hsapiens")


gene.POS <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id',
                                 'entrezgene_id',
                                 'chromosome_name',
                                 'start_position','end_position',
                                 'strand'),
                  mart = ensembl) %>% 
  #Rename
  dplyr::rename(CHR=chromosome_name) %>% 
  # genes of interest
  filter(hgnc_symbol %in% genes.OI) %>% 
  filter(!grepl("^CHR", CHR))

#DMRs in order
DMR.vec <- c("DMR_49","DMR_25","DMR_61")

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
    dplyr::select(probeID, CHR, start_hg38,
                  RSTR.M.ave, LTBI.M.ave, RSTR.M.sd, LTBI.M.sd) %>%
    pivot_longer(RSTR.M.ave:LTBI.M.sd) %>% 
    separate(name, into=c("group","measure","name"), sep="[.]") %>% 
    pivot_wider() 
  #### plot ####
  DMR_M_plot <- DMR_surround %>% 
    
    ggplot(aes(x=start_hg38, y=ave)) +
    geom_errorbar(aes(ymin=ave-sd, ymax=ave+sd, color=group, 
                      width=(max(start_hg38)-min(start_hg38))/50)) +
    geom_point(aes(color=group), size=2) +
    theme_bw() +
    labs(x="", y="", title = unique(DMR_sub$DMR)) +
    scale_color_manual("", labels=c("LTBI.M.ave"="LTBI", 
                                    "RSTR.M.ave"="RSTR"),
                       values = c("#f1a340","#998ec3")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  # DMR_M_plot
  #Remove legend except last plot
  if(d != DMR.vec[length(DMR.vec)]){
    DMR_M_plot <- DMR_M_plot +
      theme(legend.position = "none")
  }
  #Add tag to first plot
  if(d == DMR.vec[1]){
    DMR_M_plot <- DMR_M_plot +
      labs(tag="A", y="Probe log2 M value")
  }
  ##### Combo #####
  plot.ls[[d]] <- DMR_M_plot
}
# wrap_plots(plot.ls)

#### Gene plots ####
plot.ls2 <- list()

for(d in DMR.vec) {
  print(d)
  #Get DMR info
  DMR_sub <- DMR.OI2 %>% dplyr::filter(DMR==d) %>% 
    mutate(DMR = gsub("_"," ",DMR)) %>% 
    distinct(DMR, CHR, DMR_start, DMR_end, DMR_genes)
  
  #Get gene info
  gene_sub <- gene.POS %>% 
    filter(hgnc_symbol %in% DMR_sub$DMR_genes) %>% 
    full_join(DMR_sub, by=c("hgnc_symbol"="DMR_genes"))
  
  #labels
  x_lab <- paste("Chromosome", unique(DMR_sub$CHR), "position")
  
  #set DMR label position
  if(gene_sub$strand == 1){
    gene_sub_format <- gene_sub %>% 
      mutate(start = start_position,
             end = end_position,
             g_lab = DMR_start,
             g_hjust = -0.2)
  } else{
    gene_sub_format <- gene_sub %>% 
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
              fill="grey70") +
    # geom_text(aes(x=g_lab, label=DMR, hjust=g_hjust), y=0.95) +
    #Gene
    geom_segment(aes(x=start, xend=end), 
                 y=1, yend=1, color="black", size=2,
                 arrow=arrow(length = unit(0.2, "npc"))) +
    geom_text(aes(x=(start+end)/2, label = hgnc_symbol),
              y=1, hjust=0.5, vjust=-1) +
    
    theme_classic() +
    lims(y=c(0.9,1.1)) +
    labs(x=x_lab) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.margin = unit(c(0,1,0,0.5), "cm"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  #Add tag to first plot
  if(d == DMR.vec[1]){
    gene_plot <- gene_plot +
      labs(tag="B")
  }
  # gene_plot
  
  ##### Save #####
  plot.ls2[[d]] <- gene_plot
}

#### Expression correlation ####
attach("Methyl/data_Methyl_clean/eQTM.data.RData")

corr.dat <- DMR.M.overlap %>% 
  filter(DMR %in% DMR.vec) %>% 
  inner_join(rna.overlap) %>% 
  inner_join(distinct(DMR.OI, DMR, DMR_genes),
             by=c("DMR","gene"="DMR_genes")) %>% 
  left_join(dplyr::select(meta.overlap,-Sample_Name)) %>%
  mutate(DMR = factor(gsub("_"," ",DMR), levels=gsub("_"," ",DMR.vec)))

corr.dat.summ <- corr.dat %>% 
  group_by(DMR) %>% 
  summarise(minM=max(methyl)-0.1,
            minE=min(expression),
            .groups = "drop") %>%
  full_join(distinct(corr.dat, DMR, condition)) %>% 
  mutate(minE = ifelse(condition=="TB",minE-0.2, minE))

corr.result <- read_csv("Methyl/results/eQTM_correlation.csv") %>% 
  filter(DMR %in% DMR.vec) %>% 
  mutate(DMR = factor(gsub("_"," ",DMR), levels=gsub("_"," ",DMR.vec))) %>% 
  mutate(lab = paste0("R = ",round(R,2),", p = ",round(pval,2))) %>% 
  full_join(corr.dat.summ)

plot.ls3 <- list()
for(d in DMR.vec){
  corr.sub <- corr.dat %>% 
    filter(DMR == gsub("_"," ", d))
  corr.result.sub <- corr.result %>% 
    filter(DMR == gsub("_"," ", d))
  
  corr.plot <- corr.sub %>% 
    ggplot() +
    aes(x = methyl, y=expression, color = condition) +
    geom_point(size=2) +
    geom_smooth(method = "lm", aes(fill=condition), formula = y ~ x,
                show.legend = FALSE, alpha=0.1, se=FALSE, size=0.5) +
    #Correlation
    geom_text(data=corr.result.sub,
              aes(label=lab, x=minM, y=minE), show.legend = FALSE) +
    theme_bw() +
    #Recolor
    scale_color_manual("", labels=c("MEDIA"="Media", 
                                    "TB"="+Mtb"),
                       values = c("#0571b0","#ca0020")) +
    scale_fill_manual("", labels=c("MEDIA"="Media", 
                                   "TB"="+Mtb"),
                      values = c("#0571b0","#ca0020")) +
    #labels etc
    labs(x=paste(corr.result.sub$DMR, "mean log2 M value", sep="\n"), 
         y=paste(corr.result.sub$gene,"log2 gene expression", sep="\n")) +
    theme(plot.margin = unit(c(0,0,0,0.5), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  #Remove legend except last plot
  if(d != DMR.vec[length(DMR.vec)]){
    corr.plot <- corr.plot +
      theme(legend.position = "none")
  }
  #Add tag to first plot
  if(d == DMR.vec[1]){
    corr.plot <- corr.plot +
      labs(tag="C")
  }
  plot.ls3[[d]] <- corr.plot
}

#### Save ####
p_all <- wrap_plots(plot.ls) / wrap_plots(plot.ls2) / wrap_plots(plot.ls3) +
  plot_layout(heights = c(8,2,8))
# p_all

ggsave(p_all, filename = "publication/Fig3.DMR.png", width=18, height=9)
ggsave(p_all, filename = "publication/Fig3.DMR.pdf", width=18, height=9)

#### DMR probes are DMP? ####
DMP.probe <- DMP %>% 
  filter(FDR <= 0.2)  %>% 
  pull(probeID) %>% unique()

is.dmp <- DMR_all %>% 
  mutate(DMP = ifelse(probeID %in% DMP.probe, "Y","N")) %>% 
  # select(DMR,probeID,DMR_genes,DMP) %>% 
  count(DMR,DMR_genes,DMP) %>% 
  pivot_wider(names_from = DMP, values_from = n)
