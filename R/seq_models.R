
oneHot<-function (sequence_string, nucleotide = c("A", "C", "G", "T")){
     sequence_string <- toupper(sequence_string)
     seq_len <- nchar(sequence_string)
     seq_split <- unlist(strsplit(x = sequence_string, split = ""))
     seq_mat <- matrix(data = rep(0, seq_len*(length(nucleotide))), nrow = length(nucleotide))
     rownames(seq_mat) <- nucleotide
     colnames(seq_mat) <- seq_split
     for (i in nucleotide) {
          seq_mat[rownames(seq_mat) == i, colnames(seq_mat) == 
                       i] <- 1
     }
     return(seq_mat)
}

create_merged_data_seq <- function(seq_list, TF_labels, nuc, sequence_length, num_cores){
     encode_vec_length <- length(nuc)*(sequence_length)
     datasets_enc <- lapply(1:length(seq_list), function(s){
          seq_data <- as.data.frame(get(seq_list[s]))
          # print(dim(seq_data))
          # remove sequences with N
          n_id <- lapply(1:nrow(seq_data), function(z){
               grep("N", seq_data[z,])
          })
          n_id1<-grep("1", n_id)
          if(length(n_id1) >=1){
               seq_data<- as.data.frame(seq_data[-n_id1,])
          } else {
               seq_data <- seq_data
          }
          # print(dim(seq_data))
          encoded <- mclapply(1:nrow(seq_data), function(y){
               nuc_coded <- oneHot(seq_data[y,], nucleotide = nuc)
               c(nuc_coded)
          }, mc.cores = num_cores)
          len_enc <- lapply(1:length(encoded),function(z){
               length(encoded[[z]])
          })
          encoded <- encoded[grep(paste(encode_vec_length),len_enc)]
          encoded <- do.call(rbind, encoded)
          encoded <- cbind(encoded, rep(TF_labels[s], nrow(encoded)))
     })
     datasets_enc <- do.call(rbind, datasets_enc)
     
     coln<-c()
     for(one_pos in 1:(sequence_length)){
          for(n in 1:length(nuc)){
               coln<-c(coln,paste("pos",one_pos,nuc[n],sep="_"))
          }
     }
     colnames(datasets_enc) <- c(coln, "TFlabel")
     
     datasets_enc <- as.data.frame(datasets_enc)
     
     datasets_enc
     
}

seq_model <- function(tf_n, merged_data, n_threads){
     set.seed(12256803)  # Setting seed for reproducibility
     trainIndex <- createDataPartition(merged_data$TFlabel, p = 0.7, list = FALSE)
     
     # split data into test and train
     trainData <- merged_data[trainIndex, ]
     testData <- merged_data[-trainIndex, ]
     
     # prepare features and labels
     X_train <- trainData[, -ncol(trainData)]
     y_train <- trainData[, ncol(trainData)]
     X_train_matrix <- as.matrix(X_train)
     X_test <- testData[,-ncol(testData)]
     y_test <- testData[,ncol(testData)]
     
     # Create a LightGBM dataset from the matrix
     trainData <- lgb.Dataset(data = X_train_matrix, label = y_train)
     
     # try multiple learning rates to check which performs better
     learn_rate <- seq(0.05,0.5, 0.05)
     all_res <- lapply(1:length(learn_rate), function(l){
          # Define LightGBM parameters
          params <- list(
               objective = "binary",
               metric = "binary_logloss",
               learning_rate = learn_rate[l],
               num_threads = n_threads
          )
          
          # cross-validation parameters
          cv_result <- lgb.cv(
               params = params,
               data = trainData,
               nrounds = 100,           # Number of boosting rounds
               nfold = 5,               # Number of folds in cross-validation
               stratified = TRUE,       # Stratified sampling in folds
               verbose = 1,             # Print out information
               early_stopping_rounds = 10 # Stop if no improvement in 10 rounds
          )
          
          best_iter <- cv_result$best_iter
          print(paste("Best number of boosting rounds:", best_iter))
          
          # Train a LightGBM model
          model <- lgb.train(
               params = params,
               data = trainData,
               nrounds = best_iter
          )
          
          # generate predictions
          tr_predictions <- predict(model, as.matrix(X_train))
          te_predictions <- predict(model, as.matrix(X_test))
          
          pred_names <- c("tr_predictions", "te_predictions")
          y_vals <- c("y_train", "y_test")
          
          all_pred_res <- c()
          for(p in 1:length(pred_names)){
               predictions <- get(pred_names[p])
               actual <- get(y_vals[p])
               y_true <- get(y_vals[p])
               y_prob <- predictions
               
               # ROC Curve and AUC
               roc_obj <- roc(y_true, y_prob)
               best_t <- unlist(coords(roc_obj, "best", ret = "threshold"))
               
               plot.roc(roc_obj, col = "blue4", lwd = 3, print.auc = T, cex.axis=1.7, 
                        cex.lab=1.7, main = paste(tf_n, pred_names[p], "learning rate", learn_rate[l], sep=" "))
               
               # Calculate and print AUC
               auc_value <- auc(roc_obj)
               print(paste("AUC:", auc_value))
               
               y_pred <- ifelse(predictions >= best_t[1], 1, 0)
               conf_matrix <- confusionMatrix(as.factor(y_pred), as.factor(get(y_vals[p])))
               # print(conf_matrix)
               # Plot Confusion Matrix
               cm <- as.table(conf_matrix$table)
               heatmap(cm, Rowv = NA, Colv = NA, scale = "none", col = cm.colors(2), 
                       margins = c(5,5))
               
               
               
               # Precision-Recall Curve
               pr_curve <- pr.curve(scores.class0 = y_prob, weights.class0 = y_true, 
                                    curve = TRUE)
               plot(pr_curve, col = "blue4", lwd = 3, cex.axis=1.7, 
                    cex.lab=1.7, main = paste(tf_n, pred_names[p], "learning rate", 
                                              learn_rate[l], "Precision-Recall Curve", sep=" "))
               
               # F1 Scores and precision curves for thresholds
               thresholds <- seq(0, 1, by = 0.01)
               f1_scores <- sapply(thresholds, function(t) {
                    y_pred <- ifelse(y_prob >= t, 1, 0)
                    tp <- sum(y_true == 1 & y_pred == 1)
                    fp <- sum(y_true == 0 & y_pred == 1)
                    fn <- sum(y_true == 1 & y_pred == 0)
                    
                    precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
                    recall <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
                    
                    if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
                         return(2 * (precision * recall) / (precision + recall))
                    } else {
                         return(NA)
                    }
               })
               
               precision_scores <- sapply(thresholds, function(t) {
                    y_pred <- ifelse(y_prob >= t, 1, 0)
                    tp <- sum(y_true == 1 & y_pred == 1)
                    fp <- sum(y_true == 0 & y_pred == 1)
                    
                    precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
                    return(precision)
               })
               
               plot_data <- data.frame(
                    Threshold = thresholds,
                    F1_Score = f1_scores,
                    Precision = precision_scores
               )
               
               plot_data <- na.omit(plot_data)
               
               pl <- ggplot(plot_data, aes(x = Threshold)) +
                    geom_line(aes(y = F1_Score, color = "F1 Score"), linewidth = 2) +
                    geom_line(aes(y = Precision, color = "Precision"), linewidth = 2) +
                    labs(title = paste(tf_n, "learning rate", learn_rate[l], 
                                       "F1 Score and Precision vs. Threshold", sep=" "), x = "Threshold", 
                         y = "Score", subtitle = paste(pred_names[p])) +
                    theme_minimal() +
                    theme(text = element_text(size = 18),
                          axis.title = element_text(size = 18),
                          axis.text = element_text(size = 18),
                          plot.title = element_text(size = 14),
                          plot.background = element_rect(fill="white"),
                          legend.position = "bottom",
                          legend.text = element_text(size = 18),
                          panel.background = element_rect(fill = "white", colour = "grey95")) +
                    scale_color_manual(name = "Metric", 
                                       values = c("F1 Score" = "green4", "Precision" = "blue")) 
               print(pl)
               
               # Classification Report
               classification_report <- confusionMatrix(as.factor(y_pred), as.factor(y_true),
                                                        positive = "1")
               
               results <- t(as.data.frame(c(tf_n, paste(pred_names[p]), "nuc", learn_rate[l],
                                            best_t[1], t(as.data.frame(classification_report$byClass)), 
                                            auc_value, pr_curve$auc.integral)))
               rownames(results) <- c()
               colnames(results) <- c("TFname", "Dataset", "Model", "LearningRate", "threshold",
                                      names(classification_report$byClass), "AUCROC", "AUCPR")
               all_pred_res <- rbind(all_pred_res, results)
          }
          all_pred_res
     })
     all_res <- do.call(rbind, all_res)
     all_res
}

setwd("/home/manisha/manisha-extra/shape_states")
# setwd("~/Desktop/shape_clusters/hepg2/raw_bed/")

library(parallel)
library(pROC)
library(ROCR)
library(PRROC)
library(caret)
library(lightgbm)
library(ggplot2)

nuc <- c("A", "C", "G", "T")

# list of TF names
tf_names <- read.table("encode_tf_list_500thresh", header=F, stringsAsFactors = F)
tf_names <- as.data.frame(tf_names[,1])
# prepare alphabet files for TFBS and shuffled datasets using shape and random shape
num_cores <- 2
options(scipen=999)
nsites <- as.numeric(100000)


for(x in 46:nrow(tf_names)){
     
     # read the tf binding site sequences and used them to generate cluster names
     # on the basis of cluster labels
     set.seed(125637)
     tf_fasta <- read.table(paste(tf_names[x,],"_tfbs.fasta",sep=""),
                            header = F, stringsAsFactors = F)
     tf_seq <- as.data.frame(tf_fasta[seq(2,nrow(tf_fasta), 2),])
     if(nrow(tf_seq) > nsites){
          tf_seq <- as.data.frame(tf_seq[sample(1:nrow(tf_seq), nsites),])
     }
     
     r_fasta <- read.table(paste("shuffled_", tf_names[x,], ".fasta", sep=""),
                           header = F, stringsAsFactors = F)
     r_seq <- as.data.frame(r_fasta[seq(2,nrow(r_fasta), 2),])
     if(nrow(r_seq) > nsites){
          r_seq <- as.data.frame(r_seq[sample(1:nrow(r_seq), nsites),])
     }
     
     tf_fasta <- r_fasta <- c()
     
     tf_seq_list <- c("tf_seq", "r_seq")
     TF_labels <- c(1,0)
     
     #KEEP IN MIND THE NUMBER OF CORES
     merged_data <- create_merged_data_seq(tf_seq_list, TF_labels, nuc, 51, num_cores)
     
     pdf(paste("seq_models", paste(tf_names[x,]),  ".pdf", sep="_"))
     
     mod_res <- seq_model(paste(tf_names[x,]), merged_data, num_cores)
     
     dev.off()
     
     tf_data <- r_data <- merged_data <- c()
     # NAME THE DATA ACCORDING TO INPUT FEATURES (SHAPE OR RANDOM SHAPE)
     mod_res <- cbind(paste("sequence"), mod_res)
     write.table(mod_res, paste("seq_models_parametrization_", nsites,"sites", sep=""), 
                 row.names = F, sep="\t", quote=F, append=T, col.names = F)
     
}




