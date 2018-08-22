# run with only expression features

# import data for run with no expression features
T2D_train   <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/T2D_like_training_expression_ML_table.txt")
T2D_test    <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/T2D_like_testing_expression_ML_table.txt")

lipid_train <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/lipid_training_expression_ML_table.txt")
lipid_test  <-read.delim(file="/Users/rajadurty/Git/gene_expression_lookup/ML_tables/lipid_testing_expression_ML_table.txt")

# make datatables
T2D_train_data   <- T2D_train[,c(1,2)]  
T2D_test_data    <- T2D_test[,c(1,2)]   
lipid_train_data <- lipid_train[,c(1,2)]
lipid_test_data  <- lipid_test[,c(1,2)]

# lipid model
T2D_train_data$lipid <- (-2.911810303 + T2D_train$Liver*1.08214975 + T2D_train$Adrenal_Gland*0.467348588 + T2D_train$Adipose_Visceral_.Omentum.*0.373542229 + T2D_train$Cells_Transformed_fibroblasts*0.149519192 + T2D_train$Nerve_Tibial*0.13125702 + T2D_train$Uterus*0.050740519 + T2D_train$Bladder*0.049416367 + T2D_train$Muscle_Skeletal*-0.002444573 + T2D_train$Pancreas*-0.014148842 + T2D_train$Testis*-0.017890856 + T2D_train$Vagina*-0.046478713 + T2D_train$Stomach*-0.049542952 + T2D_train$Brain_Hypothalamus*-0.077216256 + T2D_train$Lung*-0.079275173 + T2D_train$Prostate*-0.109969625 + T2D_train$Kidney_Cortex*-0.145067525 + T2D_train$Artery_Tibial*-0.179748111 + T2D_train$Minor_Salivary_Gland*-0.184420825 + T2D_train$Whole_Blood*-0.226388792 + T2D_train$Brain_Nucleus_accumbens_.basal_ganglia.*-0.335715303 + T2D_train$Heart_Atrial_Appendage*-0.365555333)
T2D_test_data$lipid <- (-2.911810303 + T2D_test$Liver*1.08214975 + T2D_test$Adrenal_Gland*0.467348588 + T2D_test$Adipose_Visceral_.Omentum.*0.373542229 + T2D_test$Cells_Transformed_fibroblasts*0.149519192 + T2D_test$Nerve_Tibial*0.13125702 + T2D_test$Uterus*0.050740519 + T2D_test$Bladder*0.049416367 + T2D_test$Muscle_Skeletal*-0.002444573 + T2D_test$Pancreas*-0.014148842 + T2D_test$Testis*-0.017890856 + T2D_test$Vagina*-0.046478713 + T2D_test$Stomach*-0.049542952 + T2D_test$Brain_Hypothalamus*-0.077216256 + T2D_test$Lung*-0.079275173 + T2D_test$Prostate*-0.109969625 + T2D_test$Kidney_Cortex*-0.145067525 + T2D_test$Artery_Tibial*-0.179748111 + T2D_test$Minor_Salivary_Gland*-0.184420825 + T2D_test$Whole_Blood*-0.226388792 + T2D_test$Brain_Nucleus_accumbens_.basal_ganglia.*-0.335715303 + T2D_test$Heart_Atrial_Appendage*-0.365555333)
lipid_train_data$lipid <- (-2.911810303 + lipid_train$Liver*1.08214975 + lipid_train$Adrenal_Gland*0.467348588 + lipid_train$Adipose_Visceral_.Omentum.*0.373542229 + lipid_train$Cells_Transformed_fibroblasts*0.149519192 + lipid_train$Nerve_Tibial*0.13125702 + lipid_train$Uterus*0.050740519 + lipid_train$Bladder*0.049416367 + lipid_train$Muscle_Skeletal*-0.002444573 + lipid_train$Pancreas*-0.014148842 + lipid_train$Testis*-0.017890856 + lipid_train$Vagina*-0.046478713 + lipid_train$Stomach*-0.049542952 + lipid_train$Brain_Hypothalamus*-0.077216256 + lipid_train$Lung*-0.079275173 + lipid_train$Prostate*-0.109969625 + lipid_train$Kidney_Cortex*-0.145067525 + lipid_train$Artery_Tibial*-0.179748111 + lipid_train$Minor_Salivary_Gland*-0.184420825 + lipid_train$Whole_Blood*-0.226388792 + lipid_train$Brain_Nucleus_accumbens_.basal_ganglia.*-0.335715303 + lipid_train$Heart_Atrial_Appendage*-0.365555333)
lipid_test_data$lipid <- (-2.911810303 + lipid_test$Liver*1.08214975 + lipid_test$Adrenal_Gland*0.467348588 + lipid_test$Adipose_Visceral_.Omentum.*0.373542229 + lipid_test$Cells_Transformed_fibroblasts*0.149519192 + lipid_test$Nerve_Tibial*0.13125702 + lipid_test$Uterus*0.050740519 + lipid_test$Bladder*0.049416367 + lipid_test$Muscle_Skeletal*-0.002444573 + lipid_test$Pancreas*-0.014148842 + lipid_test$Testis*-0.017890856 + lipid_test$Vagina*-0.046478713 + lipid_test$Stomach*-0.049542952 + lipid_test$Brain_Hypothalamus*-0.077216256 + lipid_test$Lung*-0.079275173 + lipid_test$Prostate*-0.109969625 + lipid_test$Kidney_Cortex*-0.145067525 + lipid_test$Artery_Tibial*-0.179748111 + lipid_test$Minor_Salivary_Gland*-0.184420825 + lipid_test$Whole_Blood*-0.226388792 + lipid_test$Brain_Nucleus_accumbens_.basal_ganglia.*-0.335715303 + lipid_test$Heart_Atrial_Appendage*-0.365555333)

# T2D model
T2D_train_data$T2D <- (-2.803204557 + T2D_train$Pancreas*0.602849064 + T2D_train$Thyroid*0.447174811 + T2D_train$Adrenal_Gland*0.436789032 + T2D_train$Brain_Frontal_Cortex_.BA9.*0.269929025 + T2D_train$Small_Intestine_Terminal_Ileum*0.142602349 + T2D_train$Cells_Transformed_fibroblasts*0.060990289 + T2D_train$Testis*-0.011485718 + T2D_train$Nerve_Tibial*-0.032595738 + T2D_train$Whole_Blood*-0.205155212 + T2D_train$Minor_Salivary_Gland*-0.303130856 + T2D_train$Brain_Putamen_.basal_ganglia.*-0.343324335 + T2D_train$Lung*-0.491940496)
T2D_test_data$T2D <- (-2.803204557 + T2D_test$Pancreas*0.602849064 + T2D_test$Thyroid*0.447174811 + T2D_test$Adrenal_Gland*0.436789032 + T2D_test$Brain_Frontal_Cortex_.BA9.*0.269929025 + T2D_test$Small_Intestine_Terminal_Ileum*0.142602349 + T2D_test$Cells_Transformed_fibroblasts*0.060990289 + T2D_test$Testis*-0.011485718 + T2D_test$Nerve_Tibial*-0.032595738 + T2D_test$Whole_Blood*-0.205155212 + T2D_test$Minor_Salivary_Gland*-0.303130856 + T2D_test$Brain_Putamen_.basal_ganglia.*-0.343324335 + T2D_test$Lung*-0.491940496)
lipid_train_data$T2D <- (-2.803204557 + lipid_train$Pancreas*0.602849064 + lipid_train$Thyroid*0.447174811 + lipid_train$Adrenal_Gland*0.436789032 + lipid_train$Brain_Frontal_Cortex_.BA9.*0.269929025 + lipid_train$Small_Intestine_Terminal_Ileum*0.142602349 + lipid_train$Cells_Transformed_fibroblasts*0.060990289 + lipid_train$Testis*-0.011485718 + lipid_train$Nerve_Tibial*-0.032595738 + lipid_train$Whole_Blood*-0.205155212 + lipid_train$Minor_Salivary_Gland*-0.303130856 + lipid_train$Brain_Putamen_.basal_ganglia.*-0.343324335 + lipid_train$Lung*-0.491940496)
lipid_test_data$T2D <- (-2.803204557 + lipid_test$Pancreas*0.602849064 + lipid_test$Thyroid*0.447174811 + lipid_test$Adrenal_Gland*0.436789032 + lipid_test$Brain_Frontal_Cortex_.BA9.*0.269929025 + lipid_test$Small_Intestine_Terminal_Ileum*0.142602349 + lipid_test$Cells_Transformed_fibroblasts*0.060990289 + lipid_test$Testis*-0.011485718 + lipid_test$Nerve_Tibial*-0.032595738 + lipid_test$Whole_Blood*-0.205155212 + lipid_test$Minor_Salivary_Gland*-0.303130856 + lipid_test$Brain_Putamen_.basal_ganglia.*-0.343324335 + lipid_test$Lung*-0.491940496)

library(dplyr)
library(tidyr)
library(ROCR)

# T2D_training
T2D_training.pred <- prediction(T2D_train_data$T2D,T2D_train_data$type)
T2D_training.perf  = performance(T2D_training.pred, measure = "tpr", x.measure = "fpr")
T2D_training.auc   = performance(T2D_training.pred, measure = "auc")

# lipid_training
lipid_training.pred <- prediction(lipid_train_data$lipid,lipid_train_data$type)
lipid_training.perf  = performance(lipid_training.pred, measure = "tpr", x.measure = "fpr")
lipid_training.auc   = performance(lipid_training.pred, measure = "auc")

#training AUCs
T2D_training.auc@y.values
lipid_training.auc@y.values

#training plots
pdf('/Users/rajadurty/Git/gene_expression_lookup/training_ROC_onlyExpression.pdf', width = 5, height = 5)

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
pdf('/Users/rajadurty/Git/gene_expression_lookup/testing_ROCs_onlyExpression.pdf', width = 15, height = 5)

par(mfrow=c(1,2))
plot(T2D_model_lipid_testing.perf, main=toupper("T2D_like model on holdouts"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(T2D_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)

text(1, .5, "AUC", adj = 1)
text(1, .4, paste("T2D",":",round(T2D_model_T2D_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste("lipid",":",round(T2D_model_lipid_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

plot(lipid_model_lipid_testing.perf, main=toupper("lipid model on holdouts"), col = "purple1", cex.axis = 3, lwd = 2)
par(new=TRUE)
plot(lipid_model_T2D_testing.perf, col = "forestgreen", lwd = 2)
par(new=TRUE)
abline(a=0, b= 1, lty = 2, lwd=2)

text(1, .5, "AUC", adj = 1)
text(1, .4, paste("T2D",":",round(lipid_model_T2D_testing.auc@y.values[[1]],2),sep=""), col = "forestgreen", adj = 1)
text(1, .3, paste("lipid",":",round(lipid_model_lipid_testing.auc@y.values[[1]],2),sep=""), col = "purple1", adj = 1)

dev.off()