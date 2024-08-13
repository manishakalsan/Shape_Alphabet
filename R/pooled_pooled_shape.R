# here we are using the merged shape alphabet of all TFs and 
# background too and developing a model based on boht pooled 
# datasets for each method with k >= 4

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
shape_binary_model <- function( merged_data, cl_num) {
  # prepare training and test sets
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
    learning_rate = 0.1
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
  roc <- roc(y_test, predictions)
  auc_roc <- round(auc(roc),2)
  
  results<-c("pooled_shape",c(unlist(strsplit(all[cl_num],"_"))),rmse, mae, log_loss, accuracy, auc_roc)
  t(results)
  
}


# shape alphabet of pooled random background sequences
back_shape <- read.delim("pooled_background_id_shape_alphabet", header = F)
tfbs_shape <- read.delim("pooled_tfid_shape_alphabet", header = F)
colnames(back_shape) <- colnames(tfbs_shape) <- c("id",tfbs_shape[1,-ncol(tfbs_shape)])
tfbs_shape <- tfbs_shape[-1,-1]
back_shape <- back_shape[-1,-1]

# names of clusters
clusters<-c("P","Q","R","S","T","U","V","W")

# numbers of clusters
num_clust<-lapply(2:8, function(i){
  c(rep(i,3))
})

num_clust<-unlist(num_clust)

all<-c()
for(i in 2:8){ #generate column names
  d<-c()
  a<-paste("spec_",i, sep="")
  b<-paste("km_",i, sep="")
  c<-paste("hc_",i, sep="")
  d<-c(a,b,c)
  all<-c(all,d)
}
sequ_len<-16
all_method_res <- lapply(7:21, function(k){
  tfbs_cl <- tfbs_shape[,k+2]
  back_cl <- back_shape[,k+2]
  
  shape_list <- c("tfbs_cl", "back_cl")
  tf_labels <- c(1, 0)
  
  cluster_names<-clusters[1:num_clust[k]]
  
  # encode both datasets and prepare pooled set of tfbs and 
  # random background for training
  datasets_enc <- lapply(1:length(shape_list), function(x){
    data <- as.data.frame(get(shape_list[x]))
    enc_len<-lapply(1:nrow(data),function(x){
      nchar(as.character(data[x,]))
    })
    # remove longer than desired ones
    if(length(which(unlist(enc_len) != sequ_len)) > 0){
      data <- as.data.frame(data[-(which(unlist(enc_len) != sequ_len)),])
    } else {
      data <- data
    }
    
    encoded <- mclapply(1:nrow(data), function(y){
      nuc_coded <- oneHot(data[y,], n = c(cluster_names))
      c(nuc_coded)
    }, mc.cores = 20)
    encoded <- do.call(rbind, encoded)
    encoded <- cbind(encoded, rep(tf_labels[x], nrow(encoded)))
  })
  datasets_enc <- do.call(rbind, datasets_enc)
  
  coln<-c()
  for(one_pos in 1:sequ_len){
    for(n in 1:length(cluster_names)){
      coln<-c(coln,paste("pos",one_pos,cluster_names[n],sep="_"))
    }
  }
  colnames(datasets_enc) <- c(coln, "TFlabel")
  merged_data <- as.data.frame(datasets_enc)
  
  model_res <- shape_binary_model(merged_data, k)
  
  model_res
  
})
all_method_res_df <- do.call(rbind, all_method_res)
colnames(all_method_res_df) <- c(c("TF", "method","clusters","rmse", 
                                   "mae", "accuracy", "auc"))
write.table(all_method_res_df, "pooledTr_pooledTe_shape_binary", 
            sep="\t",row.names = F, quote= F)




