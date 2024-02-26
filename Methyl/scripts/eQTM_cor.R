# Correlation

```{r eval=FALSE, echo=FALSE}
library(Hmisc)
DMR.M.overlap.mat <- DMR.M.overlap %>% 
  dplyr::select(-Sample_Name) %>% 
  pivot_wider(names_from = FULLIDNO, values_from = methyl) %>% 
  column_to_rownames("DMR") %>% 
  t()
identical(rownames(DMR.M.overlap.mat), colnames(rna.media))
identical(rownames(DMR.M.overlap.mat), colnames(rna.tb))

corr.media <- rcorr(DMR.M.overlap.mat, t(rna.media), type="pearson")
corr.tb <- rcorr(DMR.M.overlap.mat, t(rna.tb), type="pearson")

#Calculate FDR
P.media <- as.data.frame(corr.media$P) %>% 
  rownames_to_column("gene") %>% 
  #Remove duplicate DMR-gene pairs
  filter(!grepl("^DMR_", gene)) %>% 
  select(gene, starts_with("DMR_")) %>% 
  pivot_longer(-gene, names_to = "DMR", values_to = "pval") %>% 
  group_by(DMR) %>% 
  mutate(FDR = p.adjust(pval, method="BH")) %>% 
  ungroup()
P.tb <- as.data.frame(corr.tb$P) %>% 
  rownames_to_column("gene") %>% 
  #Remove duplicate DMR-gene pairs
  filter(!grepl("^DMR_", gene)) %>% 
  select(gene, starts_with("DMR_")) %>% 
  pivot_longer(-gene, names_to = "DMR", values_to = "pval") %>% 
  group_by(DMR) %>% 
  mutate(FDR = p.adjust(pval, method="BH")) %>% 
  ungroup()

# Get corr coefficient
R.media <- as.data.frame(corr.media$r) %>% 
  rownames_to_column("gene") %>% 
  filter(!grepl("^DMR_", gene)) %>% 
  select(gene, starts_with("DMR_")) %>% 
  pivot_longer(-gene, names_to = "DMR", values_to = "R")
R.tb <- as.data.frame(corr.tb$r) %>% 
  rownames_to_column("gene") %>% 
  filter(!grepl("^DMR_", gene)) %>% 
  select(gene, starts_with("DMR_")) %>% 
  pivot_longer(-gene, names_to = "DMR", values_to = "R")

save(P.media, P.tb, R.media, R.tb, file = "results/eQTM_correlation.RData")
```

```{r eval=FALSE}
#Significant corr
##By FDR
P.signif1 <- P.media %>% 
  filter(FDR < 0.2) %>% 
  distinct(DMR, gene)
P.signif2 <- P.tb %>% 
  filter(FDR < 0.2) %>% 
  distinct(DMR, gene) %>% 
  bind_rows(P.signif1) %>% 
  distinct()

corr.media.signif <- R.media %>% 
  #join to get only signif DMR-gene pair
  inner_join(P.signif2) %>% 
  #Add significance data
  inner_join(P.media) %>% 
  mutate(condition = "MEDIA")
#Repeat for tb-infected correlations
corr.tb.signif <- R.tb %>% 
  inner_join(P.signif2) %>% 
  inner_join(P.tb) %>% 
  mutate(condition = "TB")

corr.signif <- full_join(corr.media.signif, corr.tb.signif)

write_csv(corr.signif, file = "results/eQTM_correlation_FDR0.2.csv")
save(corr.signif, P.signif1, P.signif2,
     file="results/eQTM_correlation_FDR0.2.RData")
```

```{r include=FALSE}
load("results/eQTM_correlation_FDR0.2.RData")
```

Significant in MEDIA and/or TB

```{r echo=FALSE}
corr.signif %>% 
  inner_join(P.signif1) %>% 
  mutate_if(is.numeric, ~round(.,3)) %>% 
  select(DMR, gene, condition, everything()) %>% 
  arrange(DMR, gene, condition) %>% 
  DT::datatable(options = list(pageLength=10))
```

Significant only in TB

```{r echo=FALSE}
corr.signif %>% 
  anti_join(P.signif1) %>% 
  mutate_if(is.numeric, ~round(.,3)) %>% 
  select(DMR, gene, condition, everything()) %>% 
  arrange(DMR, gene, condition) %>% 
  DT::datatable(options = list(pageLength=10))
```
