# run with no expression features

# import data for run with no expression features
T2D_train   <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/T2D_like_training_ML_table.txt")
T2D_test    <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/T2D_like_testing_ML_table.txt")

lipid_train <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/lipid_training_ML_table.txt")
lipid_test  <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/lipid_testing_ML_table.txt")

# alz_train   <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer\\alzheimers_training_ML_table.txt")
# alz_test    <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer\\alzheimers_testing_ML_table.txt")

# make datatables
T2D_train_data   <- T2D_train[,c(1,2)]  
T2D_test_data    <- T2D_test[,c(1,2)]   
lipid_train_data <- lipid_train[,c(1,2)]
lipid_test_data  <- lipid_test[,c(1,2)] 
# alz_train_data   <- alz_train[,c(1,2)]  
# alz_test_data    <- alz_test[,c(1,2)]   

#lipid model!
T2D_train_data$lipid   <- (-3.11469695 + T2D_train$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + T2D_train$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + T2D_train$GSE64233_MED1_V_final.bed*0.22376704 + T2D_train$FAT.ADIP.NUC.EnhA*0.082257626 + T2D_train$E118.H3K27ac.gappedPeak*0.021026009)
T2D_test_data$lipid    <- (-3.11469695 + T2D_test$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + T2D_test$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + T2D_test$GSE64233_MED1_V_final.bed*0.22376704 + T2D_test$FAT.ADIP.NUC.EnhA*0.082257626 + T2D_test$E118.H3K27ac.gappedPeak*0.021026009)
lipid_train_data$lipid <- (-3.11469695 + lipid_train$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + lipid_train$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + lipid_train$GSE64233_MED1_V_final.bed*0.22376704 + lipid_train$FAT.ADIP.NUC.EnhA*0.082257626 + lipid_train$E118.H3K27ac.gappedPeak*0.021026009)
lipid_test_data$lipid  <- (-3.11469695 + lipid_test$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + lipid_test$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292066 + lipid_test$GSE64233_MED1_V_final.bed*0.22376704 + lipid_test$FAT.ADIP.NUC.EnhA*0.082257626 + lipid_test$E118.H3K27ac.gappedPeak*0.021026009)

# alz_train_data$lipid   <- (-3.11469695 + alz_train$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + alz_train$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292069 + alz_train$GSE64233_MED1_V_final.bed*0.22376704 + alz_train$FAT.ADIP.NUC.EnhA*0.082257626 + alz_train$E118.H3K27ac.gappedPeak*0.021026009)
# alz_test_data$lipid    <- (-3.11469695 + alz_test$hepatocyte_H3K27Ac_DMSO.bed*0.669428047 + alz_test$HepG2_ChIP.seq_SREBF1_ENCFF001VLA.bed*0.242292069 + alz_test$GSE64233_MED1_V_final.bed*0.22376704 + alz_test$FAT.ADIP.NUC.EnhA*0.082257626 + alz_test$E118.H3K27ac.gappedPeak*0.021026009)

#T2D model!
T2D_train_data$T2D   <- (-2.913749265 + T2D_train$FOXA2*0.668153897 + T2D_train$islet_diff_CpGregions_hg19.bed*0.448044729 + T2D_train$NKX2.2*0.329245952 + T2D_train$PancIslt.EnhA*0.118296882 + T2D_train$PDX1*0.106744796 + T2D_train$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
T2D_test_data$T2D    <- (-2.913749266 + T2D_test$FOXA2*0.668153897 + T2D_test$islet_diff_CpGregions_hg19.bed*0.448044731 + T2D_test$NKX2.2*0.329245952 + T2D_test$PancIslt.EnhA*0.118296882 + T2D_test$PDX1*0.106744797 + T2D_test$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
lipid_train_data$T2D <- (-2.913749266 + lipid_train$FOXA2*0.668153897 + lipid_train$islet_diff_CpGregions_hg19.bed*0.448044731 + lipid_train$NKX2.2*0.329245952 + lipid_train$PancIslt.EnhA*0.118296882 + lipid_train$PDX1*0.106744797 + lipid_train$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
lipid_test_data$T2D  <- (-2.913749266 + lipid_test$FOXA2*0.668153897 + lipid_test$islet_diff_CpGregions_hg19.bed*0.448044731 + lipid_test$NKX2.2*0.329245952 + lipid_test$PancIslt.EnhA*0.118296882 + lipid_test$PDX1*0.106744797 + lipid_test$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
# alz_train_data$T2D   <- (-2.913749266 + alz_train$FOXA2*0.668153897 + alz_train$islet_diff_CpGregions_hg19.bed*0.448044731 + alz_train$NKX2.2*0.329245952 + alz_train$PancIslt.EnhA*0.118296882 + alz_train$PDX1*0.106744797 + alz_train$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)
# alz_test_data$T2D    <- (-2.913749266 + alz_test$FOXA2*0.668153897 + alz_test$islet_diff_CpGregions_hg19.bed*0.448044731 + alz_test$NKX2.2*0.329245952 + alz_test$PancIslt.EnhA*0.118296882 + alz_test$PDX1*0.106744797 + alz_test$prostate_gland_transcription_profiling_by_array_assay_ENCFF001SOE.bed*0.053514704)

# #alzheimers model!
# T2D_train_data$alzheimers   <- (-2.646907513 + T2D_train$DLD1_FOXO3_endogenous.bed*2.127621036 + T2D_train$E049.H3K4me1.gappedPeak*-0.025403905 + T2D_train$keratinocyte_transcription_profiling_by_array_assay_ENCFF001SOB.bed*-0.04913624 + T2D_train$E076.H3K4me1.gappedPeak*-0.067596955)
# T2D_test_data$alzheimers    <- (-2.646907513 + T2D_test$DLD1_FOXO3_endogenous.bed*2.127621036 + T2D_test$E049.H3K4me1.gappedPeak*-0.025403905 + T2D_test$keratinocyte_transcription_profiling_by_array_assay_ENCFF001SOB.bed*-0.04913624 + T2D_test$E076.H3K4me1.gappedPeak*-0.067596955)
# lipid_train_data$alzheimers <- (-2.646907513 + lipid_train$DLD1_FOXO3_endogenous.bed*2.127621036 + lipid_train$E049.H3K4me1.gappedPeak*-0.025403905 + lipid_train$keratinocyte_transcription_profiling_by_array_assay_ENCFF001SOB.bed*-0.04913624 + lipid_train$E076.H3K4me1.gappedPeak*-0.067596955)
# lipid_test_data$alzheimers  <- (-2.646907513 + lipid_test$DLD1_FOXO3_endogenous.bed*2.127621036 + lipid_test$E049.H3K4me1.gappedPeak*-0.025403905 + lipid_test$keratinocyte_transcription_profiling_by_array_assay_ENCFF001SOB.bed*-0.04913624 + lipid_test$E076.H3K4me1.gappedPeak*-0.067596955)
# alz_train_data$alzheimers   <- (-2.646907513 + alz_train$DLD1_FOXO3_endogenous.bed*2.127621036 + alz_train$E049.H3K4me1.gappedPeak*-0.025403905 + alz_train$keratinocyte_transcription_profiling_by_array_assay_ENCFF001SOB.bed*-0.04913624 + alz_train$E076.H3K4me1.gappedPeak*-0.067596955)
# alz_test_data$alzheimers    <- (-2.646907513 + alz_test$DLD1_FOXO3_endogenous.bed*2.127621036 + alz_test$E049.H3K4me1.gappedPeak*-0.025403905 + alz_test$keratinocyte_transcription_profiling_by_array_assay_ENCFF001SOB.bed*-0.04913624 + alz_test$E076.H3K4me1.gappedPeak*-0.067596955)
library(dplyr)
library(tidyr)

#loads ROCR package to calculate AUCs & ROC curves
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

# #alzheimers_training
# alz_training.pred <- prediction(alz_train_data$alzheimers,alz_train_data$type)
# alz_training.perf  = performance(alz_training.pred, measure = "tpr", x.measure = "fpr")
# alz_training.auc   = performance(alz_training.pred, measure = "auc")

#training AUCs
T2D_training.auc@y.values
lipid_training.auc@y.values
# alz_training.auc@y.values

#training plots
pdf('/Users/rajadurty/Git/gene_expression_lookup/training_ROC_noExpression.pdf', width = 5, height = 5)

# plot(alz_training.perf, main=toupper("ROC for training"), col = "orangered1", cex.axis = 3, lwd = 2)
# par(new=TRUE)
plot(lipid_training.perf, main = toupper("ROC for training"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(T2D_training.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
legend('bottomright', legend=c("lipid", "T2D_like", "random"), fill = c("purple1", "forestgreen", "black"))
# legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "random"), fill = c("orangered1", "purple1", "forestgreen", "black"))

dev.off()

#T2D models on holdouts
T2D_model_T2D_testing.pred <-  prediction(T2D_test_data$T2D,T2D_test_data$type)
T2D_model_T2D_testing.perf  = performance(T2D_model_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
T2D_model_T2D_testing.auc   = performance(T2D_model_T2D_testing.pred, measure = "auc")

T2D_model_lipid_testing.pred <-  prediction(lipid_test_data$T2D,lipid_test_data$type)
T2D_model_lipid_testing.perf  = performance(T2D_model_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
T2D_model_lipid_testing.auc   = performance(T2D_model_lipid_testing.pred, measure = "auc")

# T2D_model_alz_testing.pred <-  prediction(alz_test_data$T2D,alz_test_data$type)
# T2D_model_alz_testing.perf  = performance(T2D_model_alz_testing.pred, measure = "tpr", x.measure = "fpr")
# T2D_model_alz_testing.auc   = performance(T2D_model_alz_testing.pred, measure = "auc")

T2D_model_T2D_testing.auc@y.values
T2D_model_lipid_testing.auc@y.values
# T2D_model_alz_testing.auc@y.values

#lipid models on holdouts
lipid_model_T2D_testing.pred <-  prediction(T2D_test_data$lipid,T2D_test_data$type)
lipid_model_T2D_testing.perf  = performance(lipid_model_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
lipid_model_T2D_testing.auc   = performance(lipid_model_T2D_testing.pred, measure = "auc")

lipid_model_lipid_testing.pred <-  prediction(lipid_test_data$lipid,lipid_test_data$type)
lipid_model_lipid_testing.perf  = performance(lipid_model_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
lipid_model_lipid_testing.auc   = performance(lipid_model_lipid_testing.pred, measure = "auc")

# lipid_model_alz_testing.pred <-  prediction(alz_test_data$lipid,alz_test_data$type)
# lipid_model_alz_testing.perf  = performance(lipid_model_alz_testing.pred, measure = "tpr", x.measure = "fpr")
# lipid_model_alz_testing.auc   = performance(lipid_model_alz_testing.pred, measure = "auc")

lipid_model_T2D_testing.auc@y.values
lipid_model_lipid_testing.auc@y.values
# lipid_model_alz_testing.auc@y.values


#alzheimers models on holdouts
# alz_model_T2D_testing.pred <-  prediction(T2D_test_data$alzheimers,T2D_test_data$type)
# alz_model_T2D_testing.perf  = performance(alz_model_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
# alz_model_T2D_testing.auc   = performance(alz_model_T2D_testing.pred, measure = "auc")

# alz_model_lipid_testing.pred <-  prediction(lipid_test_data$alzheimers,lipid_test_data$type)
# alz_model_lipid_testing.perf  = performance(alz_model_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
# alz_model_lipid_testing.auc   = performance(alz_model_lipid_testing.pred, measure = "auc")

# alz_model_alz_testing.pred <-  prediction(alz_test_data$alzheimers,alz_test_data$type)
# alz_model_alz_testing.perf  = performance(alz_model_alz_testing.pred, measure = "tpr", x.measure = "fpr")
# alz_model_alz_testing.auc   = performance(alz_model_alz_testing.pred, measure = "auc")

# alz_model_T2D_testing.auc@y.values
# alz_model_lipid_testing.auc@y.values
# alz_model_alz_testing.auc@y.values

#plot it!
pdf('/Users/rajadurty/Git/gene_expression_lookup/testing_ROCs_noExpression.pdf', width = 15, height = 5)

par(mfrow=c(1,2)) # TODO: SHOULD THIS BE 2?? (was 3 originally, in Kim's version)
# plot(T2D_model_alz_testing.perf, main=toupper("T2D_like model on holdouts"), col = "orangered1", cex.axis = 3, lwd = 2)
# par(new=TRUE)
plot(T2D_model_lipid_testing.perf, main=toupper("T2D_like model on holdouts"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(T2D_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
legend('bottomright', legend=c("lipid holdout", "T2D_like holdout", "random"), fill = c("purple1", "forestgreen", "black"))
# legend('bottomright', legend=c("alzheimers holdout", "lipid holdout", "T2D_like holdout", "random"), fill = c("orangered1", "purple1", "forestgreen", "black"))

# plot(lipid_model_alz_testing.perf, main=toupper("lipid model on holdouts"), col = "orangered1", cex.axis = 3, lwd = 2)
# par(new=TRUE)
plot(lipid_model_lipid_testing.perf, main=toupper("lipid model on holdouts"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(lipid_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)
legend('bottomright', legend=c("lipid holdout", "T2D_like holdout", "random"), fill = c("purple1", "forestgreen", "black"))
# legend('bottomright', legend=c("alzheimers holdout", "lipid holdout", "T2D_like holdout", "random"), fill = c("orangered1", "purple1", "forestgreen", "black"))

# plot(alz_model_alz_testing.perf, main=toupper("Alzheimers model on holdouts"), col = "orangered1", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(alz_model_lipid_testing.perf, col = "purple1", lwd = 2)
# par(new=TRUE)
# plot(alz_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("alzheimers holdout", "lipid holdout", "T2D_like holdout", "random"), fill = c("orangered1", "purple1", "forestgreen", "black"))

dev.off()


#add GWAVA, CADD, & DeepSEA to testing roc plots --> SKIPPING THIS FOR NOW (future directions)

# T2D_GCD    <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/T2D_like_testing_scores_table_noExpression.txt")

# lipid_GCD  <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/lipid_testing_scores_table_noExpression.txt")

# alz_GCD    <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer\\alzheimers_testing_scores_table.txt")


#ROC plots!

#CADD on holdouts
# CADD_T2D_testing.pred <-  prediction(T2D_GCD$CADD,T2D_GCD$type)
# CADD_T2D_testing.perf  = performance(CADD_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
# CADD_T2D_testing.auc   = performance(CADD_T2D_testing.pred, measure = "auc")

# CADD_lipid_testing.pred <-  prediction(lipid_GCD$CADD,lipid_GCD$type)
# CADD_lipid_testing.perf  = performance(CADD_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
# CADD_lipid_testing.auc   = performance(CADD_lipid_testing.pred, measure = "auc")

# CADD_alz_testing.pred <-  prediction(alz_GCD$CADD,alz_GCD$type)
# CADD_alz_testing.perf  = performance(CADD_alz_testing.pred, measure = "tpr", x.measure = "fpr")
# CADD_alz_testing.auc   = performance(CADD_alz_testing.pred, measure = "auc")

# CADD_T2D_testing.auc@y.values
# CADD_lipid_testing.auc@y.values
# CADD_alz_testing.auc@y.values


#GWAVA on holdouts
# GWAVA_T2D_testing.pred <-  prediction(T2D_GCD$GWAVA,T2D_GCD$type)
# GWAVA_T2D_testing.perf  = performance(GWAVA_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
# GWAVA_T2D_testing.auc   = performance(GWAVA_T2D_testing.pred, measure = "auc")

# GWAVA_lipid_testing.pred <-  prediction(lipid_GCD$GWAVA,lipid_GCD$type)
# GWAVA_lipid_testing.perf  = performance(GWAVA_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
# GWAVA_lipid_testing.auc   = performance(GWAVA_lipid_testing.pred, measure = "auc")

# GWAVA_alz_testing.pred <-  prediction(alz_GCD$GWAVA,alz_GCD$type)
# GWAVA_alz_testing.perf  = performance(GWAVA_alz_testing.pred, measure = "tpr", x.measure = "fpr")
# GWAVA_alz_testing.auc   = performance(GWAVA_alz_testing.pred, measure = "auc")

# GWAVA_T2D_testing.auc@y.values
# GWAVA_lipid_testing.auc@y.values
# GWAVA_alz_testing.auc@y.values

#DeepSEA on holdouts
# DeepSEA_T2D_testing.pred <-  prediction(T2D_GCD$DeepSEA,T2D_GCD$type)
# DeepSEA_T2D_testing.perf  = performance(DeepSEA_T2D_testing.pred, measure = "tpr", x.measure = "fpr")
# DeepSEA_T2D_testing.auc   = performance(DeepSEA_T2D_testing.pred, measure = "auc")

# DeepSEA_lipid_testing.pred <-  prediction(lipid_GCD$DeepSEA,lipid_GCD$type)
# DeepSEA_lipid_testing.perf  = performance(DeepSEA_lipid_testing.pred, measure = "tpr", x.measure = "fpr")
# DeepSEA_lipid_testing.auc   = performance(DeepSEA_lipid_testing.pred, measure = "auc")

# DeepSEA_alz_testing.pred <-  prediction(alz_GCD$DeepSEA,alz_GCD$type)
# DeepSEA_alz_testing.perf  = performance(DeepSEA_alz_testing.pred, measure = "tpr", x.measure = "fpr")
# DeepSEA_alz_testing.auc   = performance(DeepSEA_alz_testing.pred, measure = "auc")

# DeepSEA_T2D_testing.auc@y.values
# DeepSEA_lipid_testing.auc@y.values
# DeepSEA_alz_testing.auc@y.values


#and add to plots. 
# pdf('/Users/rajadurty/Git/gene_expression_lookup/testing_ROCs_GCD_noExpression.pdf', width = 15, height = 5)

# par(mfrow=c(1,2))
# plot(alz_model_T2D_testing.perf, main=toupper("models on T2D holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(lipid_model_T2D_testing.perf, main=toupper("models on T2D holdouts"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(T2D_model_T2D_testing.perf, col = "#1b9e77", lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_T2D_testing.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_T2D_testing.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_T2D_testing.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid", "T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#7570b3", "#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))
# legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#d95f02", "#7570b3", "#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))

# plot(alz_model_lipid_testing.perf, main=toupper("models on lipid holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(lipid_model_lipid_testing.perf, main=toupper("models on lipid holdouts"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(T2D_model_lipid_testing.perf, col = "#1b9e77", lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_lipid_testing.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_lipid_testing.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_lipid_testing.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid", "T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#7570b3", "#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))
# legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#d95f02", "#7570b3", "#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))

# plot(alz_model_alz_testing.perf, main=toupper("models on Alzheimers holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(lipid_model_alz_testing.perf, main=toupper("models on Alzheimers holdouts"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(T2D_model_alz_testing.perf, col = "#1b9e77", lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_alz_testing.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_alz_testing.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_alz_testing.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#d95f02", "#7570b3", "#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))

# dev.off()

#just my models
# pdf('/Users/rajadurty/Git/gene_expression_lookup/testing_ROCs_same_holdouts_noExpression.pdf', width = 15, height = 5)

# par(mfrow=c(1,2))
# # plot(alz_model_T2D_testing.perf, main=toupper("models on T2D holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# # par(new=TRUE)
# plot(lipid_model_T2D_testing.perf, main=toupper("models on T2D holdouts"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(T2D_model_T2D_testing.perf, col = "#1b9e77", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid", "T2D_like", "random"), fill = c("#7570b3", "#1b9e77", "black"))
# # legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "random"), fill = c("#d95f02", "#7570b3", "#1b9e77", "black"))

# # plot(alz_model_lipid_testing.perf, main=toupper("models on lipid holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# # par(new=TRUE)
# plot(lipid_model_lipid_testing.perf, main=toupper("models on lipid holdouts"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(T2D_model_lipid_testing.perf, col = "#1b9e77", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid", "T2D_like", "random"), fill = c("#7570b3", "#1b9e77", "black"))
# legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "random"), fill = c("#d95f02", "#7570b3", "#1b9e77", "black"))

# plot(alz_model_alz_testing.perf, main=toupper("models on Alzheimers holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(lipid_model_alz_testing.perf, col = "#7570b3", lwd = 2)
# par(new=TRUE)
# plot(T2D_model_alz_testing.perf, col = "#1b9e77", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("alzheimers", "lipid", "T2D_like", "random"), fill = c("#d95f02", "#7570b3", "#1b9e77", "black"))

# dev.off()

# #just on target model & scores
# pdf('/Users/rajadurty/Git/gene_expression_lookup/testing_ROCs_GCD_only_noExpression.pdf', width = 15, height = 5)

# par(mfrow=c(1,2))
# plot(T2D_model_T2D_testing.perf, main=toupper("models on T2D holdouts"), col = "#1b9e77", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_T2D_testing.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_T2D_testing.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_T2D_testing.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))

# plot(lipid_model_lipid_testing.perf, main=toupper("models on lipid holdouts"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_lipid_testing.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_lipid_testing.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_lipid_testing.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#7570b3", "#66a61e", "#e7298a", "#e6ab02", "black"))

# plot(alz_model_alz_testing.perf, main=toupper("models on Alzheimers holdouts"), col = "#d95f02", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_alz_testing.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_alz_testing.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_alz_testing.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("alzheimers", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#d95f02", "#66a61e", "#e7298a", "#e6ab02", "black"))

# dev.off()




















# #and for training data:

# T2D_tr_GCD    <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/T2D_like_training_scores_table_noExpression.txt")

# lipid_tr_GCD  <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/lipid_training_scores_table_noExpression.txt")

# # alz_tr_GCD    <-read.delim(file="C:\\Users\\Kim\\Documents\\ML\\V2\\NHGRI\\models\\T2D_like_lipid_alzheimer\\alzheimers_training_scores_table.txt")

# #CADD 
# CADD_T2D_training.pred <-  prediction(T2D_tr_GCD$CADD,T2D_tr_GCD$type)
# CADD_T2D_training.perf  = performance(CADD_T2D_training.pred, measure = "tpr", x.measure = "fpr")
# CADD_T2D_training.auc   = performance(CADD_T2D_training.pred, measure = "auc")

# CADD_lipid_training.pred <-  prediction(lipid_tr_GCD$CADD,lipid_tr_GCD$type)
# CADD_lipid_training.perf  = performance(CADD_lipid_training.pred, measure = "tpr", x.measure = "fpr")
# CADD_lipid_training.auc   = performance(CADD_lipid_training.pred, measure = "auc")

# # CADD_alz_training.pred <-  prediction(alz_tr_GCD$CADD,alz_tr_GCD$type)
# # CADD_alz_training.perf  = performance(CADD_alz_training.pred, measure = "tpr", x.measure = "fpr")
# # CADD_alz_training.auc   = performance(CADD_alz_training.pred, measure = "auc")

# CADD_T2D_training.auc@y.values
# CADD_lipid_training.auc@y.values
# # CADD_alz_training.auc@y.values


# #GWAVA
# GWAVA_T2D_training.pred <-  prediction(T2D_tr_GCD$GWAVA,T2D_tr_GCD$type)
# GWAVA_T2D_training.perf  = performance(GWAVA_T2D_training.pred, measure = "tpr", x.measure = "fpr")
# GWAVA_T2D_training.auc   = performance(GWAVA_T2D_training.pred, measure = "auc")

# GWAVA_lipid_training.pred <-  prediction(lipid_tr_GCD$GWAVA,lipid_tr_GCD$type)
# GWAVA_lipid_training.perf  = performance(GWAVA_lipid_training.pred, measure = "tpr", x.measure = "fpr")
# GWAVA_lipid_training.auc   = performance(GWAVA_lipid_training.pred, measure = "auc")

# # GWAVA_alz_training.pred <-  prediction(alz_tr_GCD$GWAVA,alz_tr_GCD$type)
# # GWAVA_alz_training.perf  = performance(GWAVA_alz_training.pred, measure = "tpr", x.measure = "fpr")
# # GWAVA_alz_training.auc   = performance(GWAVA_alz_training.pred, measure = "auc")

# GWAVA_T2D_training.auc@y.values
# GWAVA_lipid_training.auc@y.values
# # GWAVA_alz_training.auc@y.values

# #DeepSEA
# DeepSEA_T2D_training.pred <-  prediction(T2D_tr_GCD$DeepSEA,T2D_tr_GCD$type)
# DeepSEA_T2D_training.perf  = performance(DeepSEA_T2D_training.pred, measure = "tpr", x.measure = "fpr")
# DeepSEA_T2D_training.auc   = performance(DeepSEA_T2D_training.pred, measure = "auc")

# DeepSEA_lipid_training.pred <-  prediction(lipid_tr_GCD$DeepSEA,lipid_tr_GCD$type)
# DeepSEA_lipid_training.perf  = performance(DeepSEA_lipid_training.pred, measure = "tpr", x.measure = "fpr")
# DeepSEA_lipid_training.auc   = performance(DeepSEA_lipid_training.pred, measure = "auc")

# # DeepSEA_alz_training.pred <-  prediction(alz_tr_GCD$DeepSEA,alz_tr_GCD$type)
# # DeepSEA_alz_training.perf  = performance(DeepSEA_alz_training.pred, measure = "tpr", x.measure = "fpr")
# # DeepSEA_alz_training.auc   = performance(DeepSEA_alz_training.pred, measure = "auc")

# DeepSEA_T2D_training.auc@y.values
# DeepSEA_lipid_training.auc@y.values
# # DeepSEA_alz_training.auc@y.values


# pdf('/Users/rajadurty/Git/gene_expression_lookup/training_ROC_GCD_noExpression.pdf', width = 15, height = 5)
# par(mfrow=c(1,2))

# plot(T2D_training.perf, main=toupper("T2D Training"), col = "#1b9e77", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_T2D_training.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_T2D_training.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_T2D_training.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("T2D_like", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#1b9e77", "#66a61e", "#e7298a", "#e6ab02", "black"))

# plot(lipid_training.perf, main=toupper("Lipid Training"), col = "#7570b3", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_lipid_training.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_lipid_training.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_lipid_training.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("lipid", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#7570b3", "#66a61e", "#e7298a", "#e6ab02", "black"))

# plot(alz_training.perf, main=toupper("Alzheimer's Training"), col = "#d95f02", cex.axis = 3, lwd = 2)
# par(new=TRUE)
# plot(DeepSEA_alz_training.perf, col = "#e6ab02", lwd = 2)
# par(new=TRUE)
# plot(CADD_alz_training.perf, col = "#66a61e", lwd = 2)
# par(new=TRUE)
# plot(GWAVA_alz_training.perf, col = "#e7298a", lwd = 2)
# par(new=TRUE)
# abline(a=0, b= 1, lty = 2, lwd=2)
# legend('bottomright', legend=c("alzheimers", "CADD", "GWAVA", "DeepSEA", "random"), fill = c("#d95f02", "#66a61e", "#e7298a", "#e6ab02", "black"))

# dev.off()