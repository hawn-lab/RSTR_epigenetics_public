library(tidyverse)
library(patchwork)
library(ggvenn)

#### Data ####
#Load all model results
attach("Methyl/results/model_selection/base_model.RData")
attach("Methyl/results/model_selection/kinship_model.RData")
attach("Methyl/results/model_selection/kinship_final_model.RData")

M <- kinship_final$lmerel %>% 
  rename(epigen=gene)

M.fit <- kinship_final$lmerel.fit %>% 
  mutate(model = "~ RSTR + age + sex + kinship") %>% 
  bind_rows(mutate(base$lm.fit, model = "~ RSTR")) %>% 
  bind_rows(mutate(kinship$lmerel.fit, model = "~ RSTR + kinship")) %>% 
  rename(epigen=gene)

attach("ATACseq/results/model_selection/base_model.RData")
attach("ATACseq/results/model_selection/kinship_model.RData")
attach("ATACseq/results/model_selection/kinship_final_model.RData")

A <- kinship_final$lmerel %>% 
  rename(epigen=gene)

A.fit <- kinship_final$lmerel.fit %>% 
  mutate(model = "~ RSTR + age + sex + kinship") %>% 
  bind_rows(mutate(base$lm.fit, model = "~ RSTR")) %>% 
  bind_rows(mutate(kinship$lmerel.fit, model = "~ RSTR + kinship")) %>% 
  rename(epigen=gene)

#### Kinship ####
M.sigma <- M.fit %>% 
  distinct(epigen, model, sigma) %>% 
  pivot_wider(names_from = model, values_from = sigma) %>% 
  mutate(kin = `~ RSTR + kinship`-`~ RSTR`, 
         cov = `~ RSTR + age + sex + kinship`-`~ RSTR + kinship`) %>% 
  mutate(group = "Methylation")

A.sigma <- A.fit %>% 
  distinct(epigen, model, sigma) %>% 
  pivot_wider(names_from = model, values_from = sigma) %>% 
  mutate(kin = `~ RSTR + kinship`-`~ RSTR`, 
         cov = `~ RSTR + age + sex + kinship`-`~ RSTR + kinship`) %>% 
  mutate(group = "Accessibility")

plot1 <- bind_rows(M.sigma,A.sigma) %>% 
  ggplot(aes(y=kin, x=group)) +
  geom_violin(fill = "grey70", color = NA, adjust=1000) +
  geom_hline(yintercept = 0, lty="dashed") +
  ylim(-3.5,3.5) +
  theme_classic() +
  labs(y="Change in sigma\n<--- Better fit with kinship   Better fit without kinship --->", 
       x="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot1

#### Age +sex ####
plot2 <-  bind_rows(M.sigma,A.sigma) %>% 
  ggplot(aes(y=cov, x=group)) +
  geom_violin(fill = "grey70", color = NA, adjust=1000) +
  geom_hline(yintercept = 0, lty="dashed") +
  ylim(-3.5,3.5) +
  theme_classic() +
  labs(y="Change in sigma\n<--- Better fit with age + sex   Better fit without age + sex --->", 
       x="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot2

#### Venn data ####
#Load all model results
M.signif <- M %>% 
  filter(FDR < 0.2 & variable != "(1 | FULLIDNO)") %>% 
  mutate(variable = recode(variable, "Sample_Group"=" RSTR vs LTBI",
                           "KCHCA_AGE_YR_CURRENT"= "Age",
                           "M0_KCVSEX"="Sex"))

A.signif <- A %>% 
  filter(FDR < 0.2 & variable != "(1 | FULLIDNO)") %>% 
  mutate(variable = recode(variable, "Sample_Group"=" RSTR vs LTBI",
                           "KCHCA_AGE_YR_CURRENT"= "Age",
                           "M0_KCVSEX"="Sex"))

#### Venns ####
venn.ls <- list()

for (dat in c("A.signif","M.signif")){
  #list to hold gene vectors
  venn_dat <- list()
  dat_filter <- get(dat)
  
  if(dat == "A.signif"){lab <- "Accessibility"} else {lab <- "Methylation"}
  
  #Each variable of interest
  for (var in unique(get(dat)$variable)){
    venn_dat[[var]] <- dat_filter %>%
      dplyr::filter(variable == var) %>%
      dplyr::distinct(epigen) %>% unlist(use.names = FALSE) 
  }
  
  #Plot all venns
  venn.ls[[dat]] <- ggvenn(venn_dat, show_percentage = FALSE, 
                           fill_color = c("white","white","white"),
                           stroke_size = 0.5, text_size = 4, set_name_size = 4) +
    ggtitle(paste(lab, "FDR < 0.2")) +
    theme(plot.title = element_text(size=12))
}

#### Save plot ####
plot_AB <- plot1 + plot2 + 
  plot_annotation(tag_levels = "A")
plot_CD <- venn.ls$A.signif + venn.ls$M.signif + 
  plot_annotation(tag_levels = list(c('C', 'D'), '1'))

ggsave(plot_AB, filename="publication/FigS1AB.model.fit.png",
       height=5, width=3)
ggsave(plot_CD, filename="publication/FigS1CD.model.fit.png",
       height=5, width=5)

#### best fit numbers ####
bind_rows(M.sigma,A.sigma) %>% 
  mutate(best = ifelse(`~ RSTR + kinship`<`~ RSTR`, "~ RSTR + kinship",
                       "~ RSTR")) %>% 
  count(group, best)

bind_rows(M.sigma,A.sigma) %>% 
  mutate(best = ifelse(`~ RSTR + age + sex + kinship`<`~ RSTR + kinship`, 
                       "~ RSTR + age + sex + kinship",
                       "~ RSTR + kinship")) %>% 
  count(group, best)
