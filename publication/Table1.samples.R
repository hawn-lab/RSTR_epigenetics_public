library(tidyverse)
library(broom)

#### Epigenetic samples with kinship data ####
meta <- read_csv("data_refs/2020.11.20RSTR_Hawn_metadata.csv") 

kin <- read_csv("data_refs/kinship_Hawn_all.csv")

atac <- meta %>% 
  filter(ATACseq == "PASS" & FULLIDNO %in% kin$rowname) %>% 
  distinct(FULLIDNO, Sample_Group, M0_KCVSEX, M0_KCVAGE, KCHCA_AGE_YR_CURRENT) %>% 
  mutate(group="atac")

methyl <- meta %>% 
  filter(methylation == "PASS" & FULLIDNO %in% kin$rowname) %>% 
  distinct(FULLIDNO, Sample_Group, M0_KCVSEX, M0_KCVAGE, KCHCA_AGE_YR_CURRENT) %>% 
  mutate(group="methyl")

#### Summarise ####
dat <- bind_rows(atac,methyl)

dat %>% 
  group_by(group) %>% 
  summarise(n = n(),
            RSTR = sum(Sample_Group == "RSTR"),
            LTBI = sum(Sample_Group == "LTBI"),
            female = sum(M0_KCVSEX == "F")/n*100,
            age0 = mean(M0_KCVAGE), age0.sd = sd(M0_KCVAGE),
            age = mean(KCHCA_AGE_YR_CURRENT), age.sd = sd(KCHCA_AGE_YR_CURRENT))

#### Statistics ####
for(dataset in unique(dat$group)){
  print(dataset)
  dat.temp <- dat %>% filter(group == dataset)
  
  
  result <- table(dat.temp$Sample_Group, dat.temp$M0_KCVSEX) %>% 
    chisq.test() %>% 
    tidy()
  print("M0_KCVSEX")
  print(result)
  
  result <- t.test(dat.temp[dat.temp$Sample_Group=="RSTR",]$M0_KCVAGE,
         dat.temp[dat.temp$Sample_Group=="LTBI",]$M0_KCVAGE) %>% 
    tidy()
  print("M0_KCVAGE")
  print(result)
  
  result <- t.test(dat.temp[dat.temp$Sample_Group=="RSTR",]$KCHCA_AGE_YR_CURRENT,
                     dat.temp[dat.temp$Sample_Group=="LTBI",]$KCHCA_AGE_YR_CURRENT) %>% 
    tidy()
  print("KCHCA_AGE_YR_CURRENT")
  print(result)
  
}
