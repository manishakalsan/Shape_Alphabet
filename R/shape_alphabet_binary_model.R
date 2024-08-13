# create TF specific merged data using random and tfbs sequences   
# develop a TF specific LightGBM (binary) model and test it
# this script takes a matrix of the converted alphabet and 
# hot encodes it to create features based on the sequences 
# for clustering methods with k>=4 

library(pROC)
library(lightgbm)
library(parallel)
library(seqinr)

options(scipen = 999)

# relevant functions
oneHot<-function (s, n = c("A", "C", "G", "T")) {
     s <- toupper(s)
     sequence_len <- nchar(s)
     seq_split <- unlist(strsplit(x = s, split = ""))
     seq_mat <- matrix(data = rep(0, sequence_len*(length(n))), nrow = length(n))
     rownames(seq_mat) <- n
     colnames(seq_mat) <- seq_split
     for (i in n) {
          seq_mat[rownames(seq_mat) == i, colnames(seq_mat) == i] <- 1
     }
     return(seq_mat)
}
tf_test_data <- function(tf_id, cluster_number, sequ_len){
     test_tfbs_data<- read.table(paste("tfbs_", tf_id, "_shape_alphabet", sep=""), header = T)
     test_backg_data <- read.table(paste("background_", tf_id, "_shape_alphabet", sep=""), header = T)
     
     # extract clustering method specific columns
     tfbs_cl <- as.data.frame(test_tfbs_data[, cluster_number+2])
     back_cl <- as.data.frame(test_backg_data[, cluster_number+2])
     
     # one hot encoding
     data_names <- c("tfbs_cl", "back_cl")
     data_labels <- c(1,0)
     
     cluster_names<-clusters[1:num_clust[cluster_number]]
     coded_df <- lapply(1:length(data_names), function(dn){
          cl_data <- as.data.frame(get(data_names[dn]))
          # check length of sequences
          enc_len<-lapply(1:nrow(cl_data),function(x){
               nchar(as.character(cl_data[x,]))
          })
          # remove longer than desired ones
          if(length(which(unlist(enc_len) != sequ_len)) > 0){
               cl_data <- as.data.frame(cl_data[-(which(unlist(enc_len) != sequ_len)),])
          } else {
               cl_data <- cl_data
          }
          
          # one hot encoding for shape alphabet in binding sites
          coded <-mclapply(1:nrow(cl_data),function(loci){
               one_seq_code<-oneHot(cl_data[loci,], n=c(cluster_names))
               c(one_seq_code)
          }, mc.cores=20)
          coded_df <- do.call(rbind, coded)
          coded_df <- cbind(coded_df, rep(data_labels[dn], nrow(coded_df)))
          coded_df
     })
     tf_backg_df <- do.call(rbind, coded_df)
     
     coln<-c()
     for(one_pos in 1:as.numeric(nchar(tfbs_cl[1,1]))){
          for(one_letter in 1:length(cluster_names)){
               coln<-c(coln,paste("pos",one_pos,cluster_names[one_letter],sep="_"))
          }
     }
     colnames(tf_backg_df)<-c(coln,"TFlabel")
     tf_backg_df
}
tf_random_encoding <- function(tfbs_data, backg_data, cluster_number, sequ_len){
     cluster_names<-clusters[1:num_clust[cluster_number]]
     
     # extract clustering method specific columns
     tfbs_cl <- tfbs_data[, cluster_number+2]
     back_cl <- backg_data[, cluster_number+2]
     
     # one hot encoding
     data_names <- c("tfbs_cl", "back_cl")
     data_labels <- c(1,0)
     
     coded_df<- lapply(1:length(data_names), function(dn){
          data <- as.data.frame(get(data_names[dn]))
          # check length of sequences
          enc_len<-lapply(1:nrow(data),function(x){
               nchar(as.character(data[x,]))
          })
          # remove longer than desired ones
          if(length(which(unlist(enc_len) != sequ_len)) > 0){
               data <- as.data.frame(data[-(which(unlist(enc_len) != sequ_len)),])
          } else {
               data <- data
          }
          
          # one hot encoding for shape alphabet in binding sites
          coded <- mclapply(1:nrow(data),function(loci){
               one_seq_code<-oneHot(data[loci,], n=c(cluster_names))
               c(one_seq_code)
          }, mc.cores=20)
          coded_df <- do.call(rbind, coded)
          coded_df <- cbind(coded_df, rep(data_labels[dn], nrow(coded_df)))
          coded_df
     })
     tf_backg_df <- do.call(rbind, coded_df)
     
     coln<-c()
     for(one_pos in 1:sequ_len){
          for(one_letter in 1:length(cluster_names)){
               coln<-c(coln,paste("pos",one_pos,cluster_names[one_letter],sep="_"))
          }
     }
     colnames(tf_backg_df)<-c(coln,"TFlabel")
     tf_backg_df
     
     
}
seq_model_lgb<-function(merged_data, prefix, cl_num){
     input_ones <- merged_data[which(merged_data$TFlabel == 1), ] 
     input_zeros <- merged_data[which(merged_data$TFlabel == 0), ]  
     set.seed(100) 
     input_ones_tr_rows <- sample(1:nrow(input_ones), 0.7*nrow(input_ones)) 
     input_zeros_tr_rows <- sample(1:nrow(input_zeros), 0.7*nrow(input_ones)) 
     tr_ones <- input_ones[input_ones_tr_rows, ]  
     tr_zeros <- input_zeros[input_zeros_tr_rows, ]
     trainingData <- rbind(tr_ones, tr_zeros) 
     test_ones <- input_ones[-input_ones_tr_rows, ]
     test_zeros <- input_zeros[-input_zeros_tr_rows, ]
     testData <- rbind(test_ones, test_zeros) 
     
     # Split data into predictors (X) and target (y)
     X_train <- trainingData[, -ncol(trainingData)]
     y_train <- trainingData[, ncol(trainingData)]
     X_train_matrix <- as.matrix(X_train)
     X_test<-testData[,-ncol(testData)]
     y_test<-testData[,ncol(testData)]
     
     # Create a LightGBM dataset from the matrix
     train_data <- lgb.Dataset(data = X_train_matrix, label = y_train)
     
     # Define LightGBM parameters
     params <- list(
          objective = "binary",
          metric = "binary_logloss",
          learning_rate = 0.01
     )
     
     # Train a LightGBM model
     model <- lgb.train(params = params, data = train_data)
     
     # evaluation metrics
     predictions <- predict(model, as.matrix(X_test))
     rmse <- round( sqrt(mean((predictions - y_test)^2)), 2)
     mae <- round( mean(abs(predictions - y_test)), 2)
     
     class_labs<-lapply(1:length(predictions),function(x){
          class<-c()
          if(unlist(predictions[[x]][1]) >= 0.5){
               class<-1
          } else {
               class<-0
          }
          class
     })
     class_labs<-unlist(class_labs)
     accuracy <- round( sum(class_labs == y_test) / length(predictions), 2)
     log_loss <-round( -mean(y_test * log(predictions) + (1 - y_test) * log(1 - predictions)), 2)
     roc <- roc(as.factor(y_test), predictions)
     auc_roc <- round(auc(roc),2)
     
     results<-as.data.frame(c(prefix,c(unlist(strsplit(all[cl_num],"_"))),rmse, mae, log_loss, accuracy, auc_roc))
     t(results)
}

# vector of number of clusters
num_clust<-lapply(2:8, function(i){
     c(rep(i,3))
})
num_clust<-unlist(num_clust)

#generate column names
all<-c()
for(i in 2:8){ 
     d<-c()
     a<-paste("spec_",i, sep="")
     b<-paste("km_",i, sep="")
     c<-paste("hc_",i, sep="")
     d<-c(a,b,c)
     all<-c(all,d)
}

clusters<-c("P","Q","R","S","T","U","V","W")

# List of TFs
list_bed <- system(" ls *peaks.bed", intern = T)
list_TFs <- lapply(1:length(list_bed),function(x){
     strsplit(list_bed[x],"_peaks")[[1]][1]
})
list_TFs <- unlist(list_TFs)

meths_acc <- lapply(7:21, function(y){
     tf_model_res <- lapply(1:length(list_TFs), function(x){
          # create TF specific merged data using random and tfbs sequences   
          merged_data <- tf_test_data(list_TFs[x], y, 16)
          # develop a TF specific LightGBM (binary) model and test it
          model_res <- seq_model_lgb(as.data.frame(merged_data), list_TFs[x], y)
          model_res
     })
     tf_model_res_df <- do.call(rbind, tf_model_res)
     colnames(tf_model_res_df) <- c("TF", "method","clusters", "rmse", "mae", "log_loss",
                                    "accuracy", "auc")
     write.table(tf_model_res_df, "prediction_tfbs5k_binary", row.names = F,
                 quote = F, sep="\t", append = T, col.names = F)
     
     tf_model_res_df
     
})
meths_acc_df <- do.call(rbind, meths_acc)
write.table(meths_acc_df ,"predicion_tf_shape_binary", sep="\t")






