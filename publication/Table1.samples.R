library(tidyverse)
library(broom)

#### Epigenetic samples with kinship data ####
meta <- read_csv("data_refs/2020.11.20RSTR_Hawn_metadata.csv") 

kin <- read_csv("data_refs/kinship_Hawn_all.csv")

atac <- meta %>% 
  filter(ATACseq == "PASS" & FULLIDNO %in% kin$rowname) %>% 
  distinct(FULLIDNO, Sample_Group, M0_KCVSEX, M0_KCVAGE, KCHCA_AGE_YR_CURRENT,
           KCB_BCGSCAR, RISK_SCORE, avgBMI) %>% 
  mutate(hiv = "N") %>% 
  mutate(group="atac")

methyl <- meta %>% 
  filter(methylation == "PASS" & FULLIDNO %in% kin$rowname) %>% 
  distinct(FULLIDNO, Sample_Group, M0_KCVSEX, M0_KCVAGE, KCHCA_AGE_YR_CURRENT,
           KCB_BCGSCAR, RISK_SCORE, avgBMI) %>% 
  mutate(hiv = "N") %>% 
  mutate(group="methyl")

#### Summarise ####
dat <- bind_rows(atac,methyl)

tab <- dat %>% 
  group_by(group) %>% 
  summarise(n.mean = n(),
            RSTR.mean = sum(Sample_Group == "RSTR"),
            LTBI.mean = sum(Sample_Group == "LTBI"),
            M0_KCVSEX.mean = sum(M0_KCVSEX == "F")/n.mean*100,
            KCB_BCGSCAR.mean = sum(KCB_BCGSCAR == "Y", na.rm=TRUE)/n.mean*100,
            HIV.mean = sum(hiv == "Y")/n.mean*100,
            
            M0_KCVAGE.mean = mean(M0_KCVAGE), M0_KCVAGE.sd = sd(M0_KCVAGE),
            KCHCA_AGE_YR_CURRENT.mean = mean(KCHCA_AGE_YR_CURRENT),
            KCHCA_AGE_YR_CURRENT.sd = sd(KCHCA_AGE_YR_CURRENT),
            RISK_SCORE.mean = mean(RISK_SCORE), RISK_SCORE.sd = sd(RISK_SCORE),
            avgBMI.mean = mean(avgBMI, na.rm=TRUE), avgBMI.sd = sd(avgBMI, na.rm=TRUE)
  )

#### Statistics ####
stats <- data.frame()

for(dataset in unique(dat$group)){
  print(dataset)
  dat.temp <- dat %>% filter(group == dataset)
  
  resultS <- table(dat.temp$Sample_Group, dat.temp$M0_KCVSEX) %>% 
    chisq.test() %>% 
    tidy()
  
  resultA1 <- t.test(dat.temp[dat.temp$Sample_Group=="RSTR",]$M0_KCVAGE,
                     dat.temp[dat.temp$Sample_Group=="LTBI",]$M0_KCVAGE) %>% 
    tidy()
  
  resultA2 <- t.test(dat.temp[dat.temp$Sample_Group=="RSTR",]$KCHCA_AGE_YR_CURRENT,
                     dat.temp[dat.temp$Sample_Group=="LTBI",]$KCHCA_AGE_YR_CURRENT) %>% 
    tidy()
  
  resultR <- t.test(dat.temp[dat.temp$Sample_Group=="RSTR",]$RISK_SCORE,
                    dat.temp[dat.temp$Sample_Group=="LTBI",]$RISK_SCORE) %>% 
    tidy()
  
  resultB <- table(dat.temp$Sample_Group, dat.temp$KCB_BCGSCAR) %>% 
    chisq.test() %>% 
    tidy()
  
  resultI <- t.test(dat.temp[dat.temp$Sample_Group=="RSTR",]$avgBMI,
                    dat.temp[dat.temp$Sample_Group=="LTBI",]$avgBMI) %>% 
    tidy()
  
  stats <- data.frame(group = dataset,
                      variable = c("M0_KCVSEX", "M0_KCVAGE", "KCHCA_AGE_YR_CURRENT",
                                   "RISK_SCORE", "KCB_BCGSCAR", "avgBMI"),
                      p = c(resultS$p.value, resultA1$p.value, resultA2$p.value,
                            resultR$p.value, resultB$p.value, resultI$p.value)) %>% 
    bind_rows(stats)
}

tab.stats <- tab %>% 
  pivot_longer(-group) %>% 
  separate(name, into=c("variable","name"), sep="[.]") %>% 
  pivot_wider() %>% 
  full_join(stats)

#### Kinship ####
kin.raw <- read_csv("data_refs/kinship_Hawn_all.csv")

kinA <- kin.raw %>% 
  filter(rowname %in% atac$FULLIDNO) %>% 
  select(rowname, all_of(atac$FULLIDNO)) %>% 
  column_to_rownames()
diag(kinA) <- NA
# kinA[lower.tri(kinA)] <- NA
kinA<- as.data.frame(kinA) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  drop_na() %>% 
  mutate(group="atac")

kinM <- kin.raw %>% 
  filter(rowname %in% methyl$FULLIDNO) %>% 
  select(rowname, all_of(methyl$FULLIDNO)) %>% 
  column_to_rownames()
diag(kinM) <- NA
# kinM[lower.tri(kinM)] <- NA
kinM<- as.data.frame(kinM) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  drop_na() %>% 
  mutate(group="methyl")

kin <- bind_rows(kinA,kinM) %>% 
  mutate(first = ifelse(value > 0.5, "Y","N"),
         third = ifelse(value > 0.125, "Y","N")) %>% 
  group_by(group, rowname) %>% 
  summarise(first = sum(first == "Y"),
            third = sum(third == "Y")) %>% 
  left_join(dat %>% distinct(FULLIDNO, Sample_Group), by=c("rowname"="FULLIDNO"))

kin.summ <- kin %>% 
  group_by(group) %>% 
  summarise(first.mean = mean(first), first.sd = sd(first),
            third.mean = mean(third), third.sd = sd(third)) %>% 
  pivot_longer(-group) %>% 
  separate(name, into=c("variable","name"), sep="[.]") %>% 
  pivot_wider()

first.statA <- t.test(kin[kin$Sample_Group=="RSTR" & kin$group=="atac",]$first,
                      kin[kin$Sample_Group=="LTBI" & kin$group=="atac",]$first) %>% 
  tidy()
third.statA <- t.test(kin[kin$Sample_Group=="RSTR" & kin$group=="atac",]$third,
                      kin[kin$Sample_Group=="LTBI" & kin$group=="atac",]$third) %>% 
  tidy()
first.statM <- t.test(kin[kin$Sample_Group=="RSTR" & kin$group=="methyl",]$first,
                      kin[kin$Sample_Group=="LTBI" & kin$group=="methyl",]$first) %>% 
  tidy()
third.statM <- t.test(kin[kin$Sample_Group=="RSTR" & kin$group=="methyl",]$third,
                      kin[kin$Sample_Group=="LTBI" & kin$group=="methyl",]$third) %>% 
  tidy()

tab.kin <- data.frame(
  group = c("atac","atac","methyl","methyl"),
  variable = c("first","third","first","third"),
  p = c(first.statA$p.value, third.statA$p.value, 
        first.statM$p.value, third.statM$p.value)
  
) %>% 
  full_join(kin.summ)

#### Save ####
bind_rows(tab.stats, tab.kin) %>% 
  write_csv("publication/Table1.csv")
