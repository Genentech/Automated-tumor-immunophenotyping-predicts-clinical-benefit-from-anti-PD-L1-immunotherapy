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
set.seed(seed)

load("~/Automated_Immunophenotyping/HK_POPLAR_alltiles_100_cutoff=0.25.RData")
comb_res = comb_res[, -which(colnames(comb_res) == "ID")[1]]

POP_manual = read.csv('~/Automated_Immunophenotyping/POPLAR_with_BCOR.csv')[, c("ID", "immunophenotype")]

Xiao_dat = merge(comb_res, POP_manual, by = "ID")

Xiao_dat = Xiao_dat[!is.na(Xiao_dat$immunophenotype),]


Jeff_POP = read.xlsx('~/Automated_Immunophenotyping/20220211_poplar_bins.xlsx', sheet = 2)
Jeff_POP$ID = paste0(Jeff_POP$Sample.ID, '_', Jeff_POP$StarLiMS.ID, '_ckpos_tiles.csv')
Jeff_POP = Jeff_POP[, c(35, 16:24, 26:34, 12)]

combine_POP = Jeff_POP %>% inner_join(Xiao_dat %>% select(-"immunophenotype"), by = "ID") %>% select(-"ID")


#5 folds repeat 10 times
control <- trainControl(method = 'repeatedcv', 
                        number = 5, 
                        repeats = 3,
                        search = "grid")


#Number randomely variable selected is mtry
mtry <- round(sqrt(ncol(combine_POP)-1))
tunegrid <- expand.grid(.mtry = c((mtry-1):(mtry+5)),
                        .min.node.size = c(1:5),
                        .splitrule = "gini")
modellist <- list()

ntree_candidates <- c(100, 150, 200, 250, 300, 350, 400, 450, 500)
for (ntree in ntree_candidates){
  set.seed(seed)
  rf_training <- train(Manual_Immunophenotype ~ ., 
                       data = combine_POP, 
                       method = 'ranger', 
                       preProcess = NULL,
                       num.trees = ntree,
                       tuneGrid = tunegrid, 
                       trControl = control)
  
  key <- toString(ntree)
  modellist[[key]] <- rf_training
}



results <- resamples(modellist)
summary(results)

best_id = which.max(as.data.frame(summary(results)$statistics$Accuracy)$Mean)

best_ntree = ntree_candidates[best_id]

best_res = modellist[[best_id]]$results

set.seed(seed)
rf_training <- train(Manual_Immunophenotype ~ ., 
                     data = combine_POP, 
                     method = 'ranger', 
                     importance = "permutation",
                     preProcess = NULL,
                     num.trees = best_ntree,
                     tuneGrid = expand.grid(mtry = best_res$mtry[which.max(best_res$Accuracy)],
                                            min.node.size = best_res$min.node.size[which.max(best_res$Accuracy)], 
                                            splitrule = best_res$splitrule[which.max(best_res$Accuracy)]),
                     trControl = trainControl(classProbs = TRUE)
)

print(rf_training$results)

rfImp = varImp(rf_training, scale = F)
plot(rfImp, top = 5, main = "Feature importance ranking via permutation method") 


fea_names = rfImp$importance %>% rownames()

fea_names[fea_names == "`0.04.ckpos.bin./.tot.ckpos.tiles`"] = "proportion of tiles in CK+ regions with 0.02-0.04 CD8 density"
fea_names[fea_names == "`0.08.ckpos.bin./.tot.ckpos.tiles`"] = "proportion of tiles in CK+ regions with 0.06-0.08 CD8 density"
fea_names[fea_names == "BC_pos"] = "Bhattacharvva coefficient in CK+ regions"
fea_names[fea_names == "TIL_pos"] = "CD8-CK ratio in CK+ regions"
fea_names[fea_names == "SL_pos"] = "Sørensen index CK+ regions"

rownames(rfImp$importance) = fea_names


png('JPATH/suppl_fig1.png', width = 1800, height = 1800, res = 300)
plot(rfImp, top = 5, main = "Feature importance ranking via permutation method") 
dev.off()





rf_training$variable.importance








