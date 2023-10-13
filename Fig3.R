# UMAP for figures 3
rm(list = ls())

library('uwot')
library('dplyr')
select <- dplyr::select
library('ggplot2')
library('survival')
library("survminer")
library('openxlsx')
library('caret')
library('tidyverse')
library('rpart')
library('rpart.plot')
library('randomForest')
library('e1071')
library('openxlsx')
library('pROC')

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}


seed = 42

library("umap") # umap_0.2.9.0
COLOR <- c("#ED7D31", "#4472C4", "#A5A5A5")


load('~/Automated_Immunophenotyping/JPATH/pop_trn_model_pip4_tuneRF.RData')
fea = names(rf_training$trainingData)[-which(names(rf_training$trainingData) == ".outcome")]


## POPLAR 
load("~/Automated_Immunophenotyping/HK_POPLAR_alltiles_100_cutoff=0.25.RData")
comb_res = comb_res[, -which(colnames(comb_res) == "ID")[1]]

POP_manual = read.csv('POPLAR_with_BCOR.csv')[, c("ID", "immunophenotype")]

Xiao_dat = merge(comb_res, POP_manual, by = "ID")

Xiao_dat = Xiao_dat[!is.na(Xiao_dat$immunophenotype),]


Jeff_POP = read.xlsx('20220211_poplar_bins.xlsx', sheet = 2)
Jeff_POP$ID = paste0(Jeff_POP$Sample.ID, '_', Jeff_POP$StarLiMS.ID, '_ckpos_tiles.csv')
Jeff_POP = Jeff_POP[, c(35, 8, 16:24, 26:34, 12)]

combine_POP = Jeff_POP %>% inner_join(Xiao_dat %>% select(-"immunophenotype"), by = "ID") %>% select(-"ID") %>% 
  rename("CD8_density" = "CD8.+.CK.area./.CK.area.*.100")

set.seed(seed)
kmeans_ = kmeans(x = apply(as.matrix(combine_POP[, fea]), 2, scale), centers = 3, iter.max = 50)


TIL_pos_mean = c(mean(combine_POP$TILs_pos[kmeans_$cluster == 1]),
                 mean(combine_POP$TILs_pos[kmeans_$cluster == 2]),
                 mean(combine_POP$TILs_pos[kmeans_$cluster == 3]))

clusters = kmeans_$cluster

clusters[kmeans_$cluster == which.max(TIL_pos_mean)] = "Inflamed"
clusters[kmeans_$cluster == which.min(TIL_pos_mean)] = "Desert"
clusters[!kmeans_$cluster %in% c(which.min(TIL_pos_mean), which.max(TIL_pos_mean))] = "Excluded"


## rename the cluster name to the 3 IP category
kmeans_pred = factor(clusters, levels = c("Inflamed", "Excluded", "Desert"))


combine_POP$Kmeans = kmeans_pred





## OAK 
load("~/Automated_Immunophenotyping/HK_OAK_alltiles_100_cutoff=0.25.RData")
Xiao_OAK = comb_res[, -which(colnames(comb_res) == "ID")[1]]

Jeff_OAK = readxl::read_excel('COLLATE_OAK_HK_final_20220606153656.xlsx', sheet = 1, col_types = c(rep("text", 10), rep("numeric", 27))) %>% mutate(Final.call = `Final consensus manual phenotype`)

Jeff_OAK$ID = paste0(Jeff_OAK$Sample.ID, '_', Jeff_OAK$StarLiMS.ID, '_ckpos_tiles.csv')

Hauke_dat = read.csv('Patient-MetaData_GO28915_OAK_v070_2112120_for-Darya - PatientDataTable_GO28915_OAK_new_210828.csv')

combine_OAK = Xiao_OAK %>% inner_join(Jeff_OAK[, c(4, 14, 2,3, 18:39)], by = "ID") %>%
  inner_join(Hauke_dat %>% select(USUBJID, ACTARM, OS_MONTHS, OS_CENSOR, PFS_MONTHS, PFS_CENSOR), by = "USUBJID")  %>% 
  mutate(Final.call = factor(Final.call, levels = c("Inflamed", "Excluded", "Desert"))) %>%
  rename("CD8_density" = "CD8 + CK area / CK area * 100")

names(combine_OAK)[56:75] = gsub(" ", ".", names(combine_OAK)[56:75])

combine_OAK = combine_OAK[complete.cases(combine_OAK), ] %>% 
  select(-c("ID", "USUBJID", "Cutoff prediction","SVM prediction", 
            "ACTARM", "OS_MONTHS", "OS_CENSOR","PFS_MONTHS", "PFS_CENSOR")) %>% 
  rename("Manual_Immunophenotype" = "Final.call" )

combine_OAK = combine_OAK[, names(combine_OAK) %in% names(combine_POP)] %>% unique()


set.seed(seed)
kmeans_ = kmeans(x = apply(as.matrix(combine_OAK[, fea]), 2, scale), centers = 3, iter.max = 50)


TIL_pos_mean = c(mean(combine_OAK$TILs_pos[kmeans_$cluster == 1]),
                 mean(combine_OAK$TILs_pos[kmeans_$cluster == 2]),
                 mean(combine_OAK$TILs_pos[kmeans_$cluster == 3]))

clusters = kmeans_$cluster

clusters[kmeans_$cluster == which.max(TIL_pos_mean)] = "Inflamed"
clusters[kmeans_$cluster == which.min(TIL_pos_mean)] = "Desert"
clusters[!kmeans_$cluster %in% c(which.min(TIL_pos_mean), which.max(TIL_pos_mean))] = "Excluded"


## rename the cluster name to the 3 IP category
kmeans_pred = factor(clusters, levels = c("Inflamed", "Excluded", "Desert"))


combine_OAK$Kmeans = kmeans_pred





## IMPassion130
Jeff_IMP130 = readxl::read_xlsx('~/Automated_Immunophenotyping/20221002_IMP130_tracker.xlsx', sheet = 1, col_types = c(rep("text", 12), rep("numeric", 27))) %>%     mutate(Final.call = factor(Immunophenotype_HK, levels = c("Inflamed", "Excluded", "Desert"))) %>% 
  rename("CD8_density" = "CD8 + CK area / CK area * 100")
Jeff_IMP130 = Jeff_IMP130[-c(1:8, 10, 11, 14, 16, 17), ]
Jeff_IMP130$ID = paste0(Jeff_IMP130$Block, '_', as.numeric(Jeff_IMP130$`Patient Number`), '_ckpos_tiles.csv')
Jeff_IMP130 = Jeff_IMP130[, c(21:29, 31:39, 40, 41, 16)]; names(Jeff_IMP130)[1:18] = gsub(" ", ".", names(Jeff_IMP130)[1:18])

Xiao_IMP130 = read.csv('IMPassion130_res.csv'); Xiao_IMP130 = Xiao_IMP130[, -which(colnames(Xiao_IMP130) == "ID.1")[1]]

combine_IMP130 = merge(Jeff_IMP130, Xiao_IMP130, by = "ID") %>% 
  rename("Manual_Immunophenotype" = "Final.call" )
combine_IMP130 = combine_IMP130[complete.cases(combine_IMP130),]

tmp1 = read.csv('~/Automated_Immunophenotyping/20221002_IMP130_tracker_ANON.csv') %>% select(uni_id = hashed.uniId, HK_block = 'Specimen.._HK')
tmp2 = read.csv('~/Automated_Immunophenotyping/wo29522_clin_jea_hk_10022022.csv')

tmp3 = merge(tmp1, tmp2, by = "uni_id")

dat_IMP130 = merge(combine_IMP130 %>% mutate(HK_block = substr(ID, 1, 9)), 
                   tmp3, 
                   by = "HK_block") %>% unique()



set.seed(seed)
kmeans_ = kmeans(x = apply(as.matrix(dat_IMP130[, fea]), 2, scale), centers = 3, iter.max = 50)

TIL_pos_mean = c(mean(combine_IMP130$TILs_pos[kmeans_$cluster == 1]),
                 mean(combine_IMP130$TILs_pos[kmeans_$cluster == 2]),
                 mean(combine_IMP130$TILs_pos[kmeans_$cluster == 3]))

clusters = kmeans_$cluster

clusters[kmeans_$cluster == which.max(TIL_pos_mean)] = "Inflamed"
clusters[kmeans_$cluster == which.min(TIL_pos_mean)] = "Desert"
clusters[!kmeans_$cluster %in% c(which.min(TIL_pos_mean), which.max(TIL_pos_mean))] = "Excluded"


## rename the cluster name to the 3 IP category
kmeans_pred = factor(clusters, levels = c("Inflamed", "Excluded", "Desert"))


dat_IMP130$Kmeans = kmeans_pred
combine_IMP130 = dat_IMP130



binned_density_fea <- names(combine_POP)[2:19]
spat_fea <- names(combine_POP)[21:70]
comb_fea <- c(binned_density_fea, spat_fea)

umap_plot <- function(this_dat = combine_POP, fea = binned_density_fea){
  
  this_dat$Immunophenotype = factor(this_dat$Manual_Immunophenotype, levels = c("Inflamed", "Excluded", "Desert"))
  
  
  umap_fit <- this_dat %>% 
    select(fea) %>% 
    scale() %>%
    umap()
  
  
  umap_df <- umap_fit$layout %>%
    as.data.frame() %>%
    rename(UMAP1 = "V1",
           UMAP2 = "V2") %>%
    mutate(ID = row_number())%>%
    inner_join(this_dat %>% mutate(ID = row_number()), by="ID")
  
  
  this_fig <- ggplot(umap_df) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = Immunophenotype)) +
    scale_color_manual(values = COLOR) +
    labs(x = "UMAP-X",
         y = "UMAP-Y") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme(legend.position="top")
  
  
  
  return(this_fig)
  
}

umap_plot_Kmeans <- function(this_dat = combine_POP, fea = comb_fea){
  
  this_dat$Immunophenotype = this_dat$Kmeans
  
  
  umap_fit <- this_dat %>% 
    select(fea) %>% 
    scale() %>%
    umap()
  
  
  umap_df <- umap_fit$layout %>%
    as.data.frame() %>%
    rename(UMAP1 = "V1",
           UMAP2 = "V2") %>%
    mutate(ID = row_number())%>%
    inner_join(this_dat %>% mutate(ID = row_number()), by="ID")
  
  
  this_fig <- ggplot(umap_df) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = Immunophenotype)) +
    scale_color_manual(values = COLOR) +
    labs(x = "UMAP-X",
         y = "UMAP-Y") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme(legend.position="top")
  
  
  
  return(this_fig)
  
}


hist_plot <-  function(this_dat = combine_POP) {
  
  this_dat$Immunophenotype = factor(this_dat$Manual_Immunophenotype, levels = c("Inflamed", "Excluded", "Desert"))
  this_dat$CD8_density = log(this_dat$CD8_density)
  
  this_fig = ggplot(this_dat) +
    geom_density(aes(x = CD8_density, fill = Immunophenotype),
                 alpha = 0.75) +
    scale_fill_manual(values = c("Inflamed" = COLOR[1],
                                 "Excluded" = COLOR[2],
                                 "Desert" = COLOR[3])) +
    labs(x = "CD8 density(log scale)",
         y = "Smoothed counts") + 
    xlim(-10, 5) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    theme(legend.position="top")
  
  return(this_fig)
}




hist_POP = hist_plot(this_dat = combine_POP)
hist_OAK = hist_plot(this_dat = combine_OAK)
hist_IMP130 = hist_plot(this_dat = combine_IMP130)


## POPLAR
umap_POP_bin = umap_plot(this_dat = combine_POP, fea = binned_density_fea)
umap_POP_spat = umap_plot(this_dat = combine_POP, fea = spat_fea)
umap_POP_comb = umap_plot(this_dat = combine_POP, fea = comb_fea)
umap_POP_Kmeans = umap_plot_Kmeans(this_dat = combine_POP, fea = comb_fea)

## OAK
umap_OAK_bin = umap_plot(this_dat = combine_OAK, fea = binned_density_fea)
umap_OAK_spat = umap_plot(this_dat = combine_OAK, fea = spat_fea)
umap_OAK_comb = umap_plot(this_dat = combine_OAK, fea = comb_fea)
umap_OAK_Kmeans = umap_plot_Kmeans(this_dat = combine_OAK, fea = comb_fea)


## IMPassion130
umap_IMP130_bin = umap_plot(this_dat = combine_IMP130, fea = binned_density_fea)
umap_IMP130_spat = umap_plot(this_dat = combine_IMP130, fea = spat_fea)
umap_IMP130_comb = umap_plot(this_dat = combine_IMP130, fea = comb_fea)
umap_IMP130_Kmeans = umap_plot_Kmeans(this_dat = combine_IMP130, fea = comb_fea)








library("gridExtra")

png("~/Automated_Immunophenotyping/JPATH/Fig3_new.png", width = 5000, height = 6000, res = 300)
# Arrange multiple ggplots and print the output
fig3 = grid.arrange(hist_POP, umap_POP_bin, umap_POP_spat, umap_POP_comb,
                    hist_OAK, umap_OAK_bin, umap_OAK_spat, umap_OAK_comb,
                    hist_IMP130, umap_IMP130_bin, umap_IMP130_spat, umap_IMP130_comb,
                    umap_POP_Kmeans, umap_OAK_Kmeans, umap_IMP130_Kmeans,
                    ncol = 3, nrow = 5, 
                    layout_matrix = cbind(c(1,2,3,4,13), c(5,6,7,8,14), c(9,10,11,12,15)))
print(fig3)

dev.off()


