library(tidyverse)
library(facetscales) #(GitHub package)
# Manhattan plots
library(ggman) #(GitHub package)
set.seed(4389)

#### Load data ####
## ATAC
DAR <- read_csv("ATACseq/results/RSTR.DAR.results.csv",
                col_types = cols(CHR="c")) %>% 
  distinct(peak, CHR, peak_start, peak_end, DAR_FDR) %>% 
  #Format for combining
  dplyr::rename(ID=peak, start=peak_start, end=peak_end, FDR=DAR_FDR) %>% 
  mutate(group="Differentially accessible\nregions (DAR)")

## METHYL
## Probes
DMP <- read_csv("Methyl/results/RSTR.DMP.results.csv.gz",
                col_types = cols(CHR="c")) %>% 
  distinct(probeID, CHR, start_hg38, end_hg38, FDR) %>% 
  #Format for combining
  dplyr::rename(ID=probeID, start=start_hg38, end=end_hg38)%>% 
  mutate(group="Differentially methylated\nprobes (DMP)")

### Regions
DMR <- read_csv("Methyl/results/RSTR.DMR.results.cpgs.csv.gz",
                col_types = cols(CHR="c")) %>% 
  distinct(DMR, CHR, DMR_start, DMR_end, DMR_FDR) %>% 
  #Format for combining
  dplyr::rename(ID=DMR, start=DMR_start, end=DMR_end, FDR=DMR_FDR)%>% 
  mutate(group="Differentially methylated\nregions (DMR)")

#### Combine data ####
man.dat <- DAR %>% 
  bind_rows(DMP) %>% 
  bind_rows(DMR) %>% 
  #Remove XY and alt chromosomes
  filter(CHR %in% c(1:22)) 

#### Plot param ####
#Y scales
scales_y <- list(
  `Differentially accessible\nregions (DAR)` = scale_y_continuous(
    limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)),
  `Differentially methylated\nprobes (DMP)` = scale_y_continuous(
    limits = c(0, 3), breaks = seq(0, 3, 1)),
  `Differentially methylated\nregions (DMR)` = scale_y_continuous(
    limits = c(0, 300), breaks = seq(0, 300, 100)))

#FDR cutoff lines
hlines <- data.frame(
  group = c("Differentially accessible\nregions (DAR)",
          "Differentially methylated\nprobes (DMP)",
          "Differentially methylated\nregions (DMR)"),
  intercept = c(0.2,0.2,1E-70),
  label = c("FDR = 0.2","FDR = 0.2","FDR = 1E-70")
)

#### Plot ####
man.plot <- ggman(man.dat, snp="ID", bp="start",
                  chrom="CHR", pval="FDR", pointSize=1, logTransform=TRUE, 
                  sigLine=NA, relative.positions=TRUE, ymin = 0) +
  #FDR cutoff lines
  geom_hline(data=hlines, aes(yintercept= -log10(intercept),
                 linetype = label)) +
  scale_linetype_manual(name = "", values = c(1,2)) +
  #Separate data sets
  facet_grid_sc(rows=vars(group), scales=list(y=scales_y)) +
  #Beautify
  theme_classic(base_size = 12) +
  labs(y="-log10( FDR )", x="Chromosome", title="") +
  theme(legend.position = "bottom",
        panel.border=element_rect(color="black", size=1, fill="transparent"),
        axis.text.x = element_text(angle=45, size=7, vjust=1, hjust=1)) #

man.plot

 #### Save ####
ggsave("publication/Fig1.RSTR.LTBI.Manhattan.pdf", man.plot, height=6, width=10)
ggsave("publication/Fig1.RSTR.LTBI.Manhattan.png", man.plot, height=6, width=10)
