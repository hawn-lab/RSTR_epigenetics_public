library(tidyverse)
library(ggpubr)

load("RNAseq/RSTR_RNAseq_data_for_eQTM.RData")
stat.dat <- data.frame(
  group1 = c("MEDIA\nLTBI","MEDIA\nRSTR","TB\nLTBI"),
  group2 = c("TB\nLTBI","TB\nRSTR","TB\nRSTR"),
  y.position  = c(10.6,11,11.4),
  symbol = c("**","**","*")
)


p1 <- as.data.frame(dat.combined.voom$E) %>% 
  rownames_to_column("geneName") %>% 
  filter(geneName == "ABCA1") %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets) %>% 
  mutate(x=paste(condition, Sample_Group, sep="\n")) %>% 
  
  ggplot(aes(x=x, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2) +
  stat_pvalue_manual(stat.dat, label="symbol") +
  theme_classic() +
  labs(x="", y="Log2 normalized expression")

ggsave(p1, filename = "publication/FigX.ABCA1.png", width=2.5, height=4)
