# run with no expression features

# import data for run with no expression features
T2D_train   <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/T2D_like_training_ML_table.txt")
T2D_test    <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/T2D_like_testing_ML_table.txt")

lipid_train <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/lipid_training_ML_table.txt")
lipid_test  <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/lipid_testing_ML_table.txt")

# make datatables
T2D_train_data   <- T2D_train[,c(1,2)]  
T2D_test_data    <- T2D_test[,c(1,2)]   
lipid_train_data <- lipid_train[,c(1,2)]
lipid_test_data  <- lipid_test[,c(1,2)]  

#lipid model!
T2D_train_data$lipid   <- (-3.11469695 + T2D_train$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + T2D_train$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + T2D_train$GSE64233_MED1_V_final.bed*0.22376704 + T2D_train$FAT.ADIP.NUC.EnhA*0.082257626 + T2D_train$E118.H3K27ac.gappedPeak*0.021026009)
T2D_test_data$lipid    <- (-3.11469695 + T2D_test$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + T2D_test$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + T2D_test$GSE64233_MED1_V_final.bed*0.22376704 + T2D_test$FAT.ADIP.NUC.EnhA*0.082257626 + T2D_test$E118.H3K27ac.gappedPeak*0.021026009)
lipid_train_data$lipid <- (-3.11469695 + lipid_train$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + lipid_train$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + lipid_train$GSE64233_MED1_V_final.bed*0.22376704 + lipid_train$FAT.ADIP.NUC.EnhA*0.082257626 + lipid_train$E118.H3K27ac.gappedPeak*0.021026009)
lipid_test_data$lipid  <- (-3.11469695 + lipid_test$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + lipid_test$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + lipid_test$GSE64233_MED1_V_final.bed*0.22376704 + lipid_test$FAT.ADIP.NUC.EnhA*0.082257626 + lipid_test$E118.H3K27ac.gappedPeak*0.021026009)

#T2D model!
T2D_train_data$T2D   <- (-2.913749265 + T2D_train$FOXA2*0.668153897 + T2D_train$islet_diff_CpGregions_hg19.bed*0.448044729 + T2D_train$NKX2.2*0.329245952 + T2D_train$PancIslt.EnhA*0.118296882 + T2D_train$PDX1*0.106744796 + T2D_train$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
T2D_test_data$T2D    <- (-2.913749266 + T2D_test$FOXA2*0.668153897 + T2D_test$islet_diff_CpGregions_hg19.bed*0.448044731 + T2D_test$NKX2.2*0.329245952 + T2D_test$PancIslt.EnhA*0.118296882 + T2D_test$PDX1*0.106744797 + T2D_test$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
lipid_train_data$T2D <- (-2.913749266 + lipid_train$FOXA2*0.668153897 + lipid_train$islet_diff_CpGregions_hg19.bed*0.448044731 + lipid_train$NKX2.2*0.329245952 + lipid_train$PancIslt.EnhA*0.118296882 + lipid_train$PDX1*0.106744797 + lipid_train$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
lipid_test_data$T2D  <- (-2.913749266 + lipid_test$FOXA2*0.668153897 + lipid_test$islet_diff_CpGregions_hg19.bed*0.448044731 + lipid_test$NKX2.2*0.329245952 + lipid_test$PancIslt.EnhA*0.118296882 + lipid_test$PDX1*0.106744797 + lipid_test$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)

library(dplyr)
library(tidyr)
library(ROCR)

####################################
#ROC plots!

#T2D_training
T2D_training.pred <- prediction(T2D_train_data$T2D,T2D_train_data$type)
T2D_training.perf  = performance(T2D_training.pred, measure = "tpr", x.measure = "fpr")
T2D_training.auc   = performance(T2D_training.pred, measure = "auc")

#lipid_training
lipid_training.pred <- prediction(lipid_train_data$lipid,lipid_train_data$type)
lipid_training.perf  = performance(lipid_training.pred, measure = "tpr", x.measure = "fpr")
lipid_training.auc   = performance(lipid_training.pred, measure = "auc")

#training AUCs
T2D_training.auc@y.values
lipid_training.auc@y.values

#training plots
pdf('/Users/rajadurty/Git/gene_expression_lookup/training_ROC_noExpression.pdf', width = 5, height = 5)

plot(lipid_training.perf, main = toupper("ROC for training"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(T2D_training.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
legend('bottomright', legend=c("lipid", "T2D_like", "random"), fill = c("purple1", "forestgreen", "black"))

dev.off()

#T2D models on holdouts
T2D_model_T2D_testing.pred <-  prediction(T2D_test_data$T2D,T2D_test_data$type)
T2D_model_T2D_testing.perf  = performance(T2D_model_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
T2D_model_T2D_testing.auc   = performance(T2D_model_T2D_testing.pred, measure = "auc")

T2D_model_lipid_testing.pred <-  prediction(lipid_test_data$T2D,lipid_test_data$type)
T2D_model_lipid_testing.perf  = performance(T2D_model_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
T2D_model_lipid_testing.auc   = performance(T2D_model_lipid_testing.pred, measure = "auc")

T2D_model_T2D_testing.auc@y.values
T2D_model_lipid_testing.auc@y.values

#lipid models on holdouts
lipid_model_T2D_testing.pred <-  prediction(T2D_test_data$lipid,T2D_test_data$type)
lipid_model_T2D_testing.perf  = performance(lipid_model_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
lipid_model_T2D_testing.auc   = performance(lipid_model_T2D_testing.pred, measure = "auc")

lipid_model_lipid_testing.pred <-  prediction(lipid_test_data$lipid,lipid_test_data$type)
lipid_model_lipid_testing.perf  = performance(lipid_model_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
lipid_model_lipid_testing.auc   = performance(lipid_model_lipid_testing.pred, measure = "auc")

lipid_model_T2D_testing.auc@y.values
lipid_model_lipid_testing.auc@y.values

#plot it!
pdf('/Users/rajadurty/Git/gene_expression_lookup/testing_ROCs_noExpression.pdf', width = 15, height = 5)

par(mfrow=c(1,2))
plot(T2D_model_lipid_testing.perf, main=toupper("T2D_like model on holdouts"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(T2D_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid holdout", "T2D_like holdout", "random"), fill = c("purple1", "forestgreen", "black"))

text(1, .5, "AUC", adj = 1)
text(1, .4, paste("T2D",":",round(T2D_model_T2D_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste("lipid",":",round(T2D_model_lipid_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

plot(lipid_model_lipid_testing.perf, main=toupper("lipid model on holdouts"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(lipid_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid holdout", "T2D_like holdout", "random"), fill = c("purple1", "forestgreen", "black"))

text(1, .5, "AUC", adj = 1)
text(1, .4, paste("T2D",":",round(lipid_model_T2D_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste("lipid",":",round(lipid_model_lipid_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

dev.off()