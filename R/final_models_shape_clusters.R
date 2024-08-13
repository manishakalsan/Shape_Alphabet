# setwd("~/Desktop/shape_clusters/encode/")
# require(ggplot2)
# require(ggseqlogo)
library(parallel)
library(seqinr)
library(pROC)
library(lightgbm)

options(scipen = 999)

clusters <- c("P", "Q", "R", "S", "T", "U", "V", "W")
nuc <- c("A", "C", "G", "T")

cl_methods <- c("spec", "km", "hc")
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

# one hot encoding of shape alphabet or nucleotide sequences
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

assign_clusters <- function(seq_data, tf_name, file_name){
  # remove sequences with N
  n_grep <- grep("N", seq_data[,1])
  if(length(n_grep) >=1){
    seq_data <- as.data.frame(seq_data[-n_grep,])
  } else {
    seq_data <- as.data.frame(seq_data)
  }
  all_greps<-lapply(1:nrow(seq_data), function(nseq){
    aseq<-toupper(as.character(seq_data[nseq,1])) #sequence
    # creates pentamers for a sequence
    sub_seqs<-lapply(1:(nchar(aseq)-4), function(l){
      substr(aseq,l,as.numeric(l+4))
    })
    sub_seqs <- unlist(sub_seqs)
    k_clusters <- mclapply(2:8, function(k){
      # read clustering info files and create strings of shape alphabet
      # for each method for a value of k
      all_method_clusters <- lapply(1:length(cl_methods), function(cl){
        method_clusters <- read.table(paste("new_pc10",cl_methods[cl],k,sep="_"),header=F)
        one_seq_c<-lapply(1:length(sub_seqs), function(m){
          method_cluster_names <- clusters[method_clusters[grep(toupper(sub_seqs[m]),
                                                                method_clusters[,2]),3]] 
        })
        one_seq_c <- paste(one_seq_c, collapse = "")
      })
      all_method_clusters <- unlist(all_method_clusters)
    }, mc.cores = 20)
    k_clusters <- unlist(k_clusters)
    k_clusters <- t(as.data.frame(c(tf_name, paste("Sequence",nseq, sep=""), k_clusters)))
    colnames(k_clusters) <- c("TF","Sequence",all)
    
    k_clusters
  })
  all_greps <- do.call(rbind, all_greps)
  write.table(all_greps, paste(file_name), col.names = T, row.names = F,
              quote=F)
  print(paste("Results written to ", file_name, sep=""))
  
  all_greps
  
}

create_merged_data <- function(seq_list, TF_labels, cluster_names, sequence_length){
  encode_vec_length <- length(cluster_names)*(sequence_length-4)
  datasets_enc <- lapply(1:length(seq_list), function(x){
    seq_data <- as.data.frame(get(seq_list[x]))
    print(dim(seq_data))
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
    print(dim(seq_data))
    encoded <- mclapply(1:nrow(seq_data), function(y){
      nuc_coded <- oneHot(seq_data[y,], n=cluster_names)
      c(nuc_coded)
    }, mc.cores = 20)
    len_enc <- lapply(1:length(encoded),function(z){
      length(encoded[[z]])
    })
    encoded <- encoded[grep(paste(encode_vec_length),len_enc)]
    encoded <- do.call(rbind, encoded)
    encoded <- cbind(encoded, rep(TF_labels[x], nrow(encoded)))
  })
  datasets_enc <- do.call(rbind, datasets_enc)
  
  coln<-c()
  for(one_pos in 1:(sequence_length-4)){
    for(n in 1:length(cluster_names)){
      coln<-c(coln,paste("pos",one_pos,cluster_names[n],sep="_"))
    }
  }
  colnames(datasets_enc) <- c(coln, "TFlabel")
  
  merged_data <- as.data.frame(datasets_enc)
  
}

create_merged_data_seq <- function(seq_list, TF_labels){
  datasets_enc <- lapply(1:length(seq_list), function(x){
    seq_data <- as.data.frame(get(seq_list[x]))
    print(dim(seq_data))
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
    print(dim(seq_data))
    encoded <- mclapply(1:nrow(seq_data), function(y){
      nuc_coded <- oneHot(seq_data[y,])
      c(nuc_coded)
    }, mc.cores = 20)
    encoded <- do.call(rbind, encoded)
    encoded <- cbind(encoded, rep(TF_labels[x], nrow(encoded)))
  })
  datasets_enc <- do.call(rbind, datasets_enc)
  
  coln<-c()
  for(one_pos in 1:20){
    for(n in 1:length(nuc)){
      coln<-c(coln,paste("pos",one_pos,nuc[n],sep="_"))
    }
  }
  colnames(datasets_enc) <- c(coln, "TFlabel")
  
  merged_data <- as.data.frame(datasets_enc)
  
}

create_merged_data_shape <- function(seq_list, TF_labels, cluster_names, sequence_length){
  encode_vec_length <- length(cluster_names)*(sequence_length-4)
  datasets_enc <- lapply(1:length(seq_list), function(x){
    seq_data <- as.data.frame(get(seq_list[x]))
    print(dim(seq_data))
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
    print(dim(seq_data))
    encoded <- mclapply(1:nrow(seq_data), function(y){
      nuc_coded <- oneHot(seq_data[y,], n=cluster_names)
      c(nuc_coded)
    }, mc.cores = 20)
    len_enc <- lapply(1:length(encoded),function(z){
      length(encoded[[z]])
    })
    encoded <- encoded[grep(paste(encode_vec_length),len_enc)]
    encoded <- do.call(rbind, encoded)
    encoded <- cbind(encoded, rep(TF_labels[x], nrow(encoded)))
  })
  datasets_enc <- do.call(rbind, datasets_enc)

  coln<-c()
  for(one_pos in 1:(sequence_length-4)){
    for(n in 1:length(cluster_names)){
      coln<-c(coln,paste("pos",one_pos,cluster_names[n],sep="_"))
    }
  }
  colnames(datasets_enc) <- c(coln, "TFlabel")

  merged_data <- as.data.frame(datasets_enc)

}

pooled_cv_model <- function(merged_data, tf_data, num_parts){ #shape_pooled_spec_model
  
  data_list <- c("merged_data", "tf_data")
  parts_list <- lapply(1:length(data_list), function(z){
    data <- get(data_list[z])
    input_ones <- data[which(data$TFlabel == 1), ] 
    input_zeros <- data[which(data$TFlabel == 0), ]  
    
    set.seed(100) 
    parts_ones <- split(input_ones, cut(sample(seq(nrow(input_ones))), breaks = num_parts, labels = FALSE))
    parts_zeroes <- split(input_zeros, cut(sample(seq(nrow(input_zeros))), breaks = num_parts, labels = FALSE))
    
    list(parts_ones, parts_zeroes)
  })
  
  # create test and train tests
  all_predictions <- c()
  for(l in 1:num_parts){
    train_data <- rbind(do.call(rbind, parts_list[[1]][[1]][-l]),
                        do.call(rbind, parts_list[[1]][[2]][-l]), 
                        do.call(rbind, parts_list[[2]][[1]][-l]),
                        do.call(rbind, parts_list[[2]][[2]][-l]))
    
    test_data <- rbind(as.data.frame(parts_list[[1]][[1]][[l]]),
                       as.data.frame(parts_list[[1]][[2]][[l]]), 
                       as.data.frame(parts_list[[2]][[1]][[l]]),
                       as.data.frame(parts_list[[2]][[2]][[l]]))
    
    
    # prepare features and labels
    X_train <- train_data[, -ncol(train_data)]
    y_train <- train_data[, ncol(train_data)]
    X_train_matrix <- as.matrix(X_train)
    X_test <- test_data[,-ncol(test_data)]
    y_test <- test_data[,ncol(test_data)]
    
    # Create a LightGBM dataset from the matrix
    train_data <- lgb.Dataset(data = X_train_matrix, label = y_train)
    
    # Define LightGBM parameters
    params <- list(
      objective = "binary",
      metric = "binary_logloss",
      learning_rate = 0.1,
      num_threads = 20
    )
    
    # Train a LightGBM model
    model <- lgb.train(params = params, data = train_data)
    
    # generate predictions
    tr_predictions <- predict(model, as.matrix(X_train))
    te_predictions <- predict(model, as.matrix(X_test))
    
    tr_pred_labels <- class_labs(tr_predictions)
    te_pred_labels <- class_labs(te_predictions)
    
    res_data <- rbind(cbind(y_train ,tr_pred_labels, rep("train", length(tr_predictions))),
          cbind(y_test, te_pred_labels, rep("test", length(te_predictions))))
    all_predictions <- rbind(all_predictions, res_data)
    
  }
  
  train_predictions <- all_predictions[all_predictions[,3] == "train",]
  test_predictions <- all_predictions[all_predictions[,3] == "test",]
  
  te_accuracy <- round(sum(test_predictions[,1] == test_predictions[,2]) / nrow(test_predictions), 2)
  tr_accuracy <- round(sum(train_predictions[,1] == train_predictions[,2]) / nrow(train_predictions), 2)
  
  c(tr_accuracy, te_accuracy)
}

shape_model_lgb <-function(merged_data, prefix, m){
  input_ones <- merged_data[which(merged_data$TFlabel == 1), ] 
  input_zeros <- merged_data[which(merged_data$TFlabel == 0), ]  
  set.seed(100) 
  input_ones_tr_rows <- sample(1:nrow(input_ones), 0.7*nrow(input_ones)) 
  input_zeros_tr_rows <- sample(1:nrow(input_zeros), 0.7*nrow(input_zeros)) 
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
    learning_rate = 0.1,
    num_threads = 20
  )
  
  # Train a LightGBM model
  model <- lgb.train(params = params, data = train_data)
  
  # evaluation metrics
  predictions <- predict(model, as.matrix(X_test))
  
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
  
  results<-as.data.frame(c(prefix,c(unlist(strsplit(all[m],"_"))),accuracy))
  t(results)
}

pooled_spec_model <- function(merged_data, test_data, prefix){ #shape_pooled_spec_model
  input_ones <- merged_data[which(merged_data$TFlabel == 1), ] 
  input_zeros <- merged_data[which(merged_data$TFlabel == 0), ]  
  set.seed(100) 
  input_ones_tr_rows <- sample(1:nrow(input_ones), 0.7*nrow(input_ones)) 
  input_zeros_tr_rows <- sample(1:nrow(input_zeros), 0.7*nrow(input_ones)) 
  tr_ones <- input_ones[input_ones_tr_rows, ]  
  tr_zeros <- input_zeros[input_zeros_tr_rows, ]
  trainingData <- rbind(tr_ones, tr_zeros) 
  
  # Split data into predictors (X) and target (y)
  X_train <- trainingData[, -ncol(trainingData)]
  y_train <- trainingData[, ncol(trainingData)]
  X_train_matrix <- as.matrix(X_train)
  
  # prepare test data
  X_test <- as.matrix(test_data[,-ncol(test_data)])
  y_test <- as.matrix(test_data[,ncol(test_data)])
  
  # Create a LightGBM dataset from the matrix
  train_data <- lgb.Dataset(data = X_train_matrix, label = y_train)
  
  # Define LightGBM parameters
  params <- list(
    objective = "binary",
    metric = "binary_logloss",
    learning_rate = 0.1,
    num_threads = 20
  )
  
  # Train a LightGBM model
  model <- lgb.train(params = params, data = train_data)
  
  # evaluation metrics
  predictions <- predict(model, X_test)
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
  
  results <- as.data.frame(c(prefix, accuracy))
  t(results)
}

class_labs <- function(predictions){
  pred_labels <- lapply(1:length(predictions),function(x){
    class<-c()
    if(unlist(predictions[[x]][1]) >= 0.5){
      class<-1
    } else {
      class<-0
    }
    class
  })
  pred_labels<-unlist(pred_labels)
}

accuracy_parts <- function(input_ones, input_zeros){
  set.seed(100) 
  parts_ones <- split(input_ones, cut(sample(seq(nrow(input_ones))), breaks = num_parts, labels = FALSE))
  parts_zeroes <- split(input_zeros, cut(sample(seq(nrow(input_zeros))), breaks = num_parts, labels = FALSE))
  
  # create test and train tests
  all_predictions <- c()
  for(l in 1:num_parts){
    feature_names <- colnames(parts_ones[[l]])
    # test data specific to a fold
    test_data <- rbind(as.data.frame(parts_ones[l]), 
                       as.data.frame(parts_zeroes[l]))
    colnames(test_data) <- feature_names
    # create train data leaving one fold data
    train_data <- rbind(do.call(rbind, parts_ones[-l]), 
                        do.call(rbind, parts_zeroes[-l]))
    # prepare features and labels
    X_train <- train_data[, -ncol(train_data)]
    y_train <- train_data[, ncol(train_data)]
    X_train_matrix <- as.matrix(X_train)
    X_test <- test_data[,-ncol(test_data)]
    y_test <- test_data[,ncol(test_data)]
    
    # Create a LightGBM dataset from the matrix
    train_data <- lgb.Dataset(data = X_train_matrix, label = y_train)
    
    # Define LightGBM parameters
    params <- list(
      objective = "binary",
      metric = "binary_logloss",
      learning_rate = 0.1,
      num_threads = 20
    )
    
    # Train a LightGBM model
    model <- lgb.train(params = params, data = train_data)
    
    # generate predictions
    tr_predictions <- predict(model, as.matrix(X_train))
    te_predictions <- predict(model, as.matrix(X_test))
    
    tr_pred_labels <- class_labs(tr_predictions)
    te_pred_labels <- class_labs(te_predictions)
    
    data_res <- rbind(cbind(y_train ,tr_pred_labels, rep("train", length(tr_predictions))),
          cbind(y_test, te_pred_labels, rep("test", length(te_predictions))))
    all_predictions <- rbind(all_predictions, data_res)
    
  }
  
  train_predictions <- all_predictions[all_predictions[,3] == "train",]
  test_predictions <- all_predictions[all_predictions[,3] == "test",]
  
  te_accuracy <- round(sum(test_predictions[,1] == test_predictions[,2]) / nrow(test_predictions), 2)
  tr_accuracy <- round(sum(train_predictions[,1] == train_predictions[,2]) / nrow(train_predictions), 2)
  
  c(tr_accuracy, te_accuracy)
}

seq_model_parts <- function(merged_data, num_parts){
  input_ones <- merged_data[which(merged_data$TFlabel == 1), ] 
  input_zeros <- merged_data[which(merged_data$TFlabel == 0), ] 
  
  set.seed(100) 
  parts_ones <- split(input_ones, cut(sample(seq(nrow(input_ones))), breaks = num_parts, labels = FALSE))
  parts_zeroes <- split(input_zeros, cut(sample(seq(nrow(input_zeros))), breaks = num_parts, labels = FALSE))
  
  # create test and train tests
  all_predictions <- c()
  for(l in 1:num_parts){
    print(l)
    feature_names <- colnames(parts_ones[[l]])
    # test data specific to a fold
    test_data <- rbind(as.data.frame(parts_ones[l]), 
                       as.data.frame(parts_zeroes[l]))
    colnames(test_data) <- feature_names
    # create train data leaving one fold data
    train_data <- rbind(do.call(rbind, parts_ones[-l]), 
                        do.call(rbind, parts_zeroes[-l]))
    # prepare features and labels
    X_train <- train_data[, -ncol(train_data)]
    y_train <- train_data[, ncol(train_data)]
    X_train_matrix <- as.matrix(X_train)
    X_test <- test_data[,-ncol(test_data)]
    y_test <- test_data[,ncol(test_data)]
    
    # Create a LightGBM dataset from the matrix
    train_data <- lgb.Dataset(data = X_train_matrix, label = y_train)
    
    # Define LightGBM parameters
    params <- list(
      objective = "binary",
      metric = "binary_logloss",
      learning_rate = 0.1,
      num_threads = 20
    )
    
    # Train a LightGBM model
    model <- lgb.train(params = params, data = train_data)
    
    # generate predictions
    tr_predictions <- predict(model, as.matrix(X_train))
    te_predictions <- predict(model, as.matrix(X_test))
    
    tr_pred_labels <- class_labs(tr_predictions)
    te_pred_labels <- class_labs(te_predictions)
    print(c("processed ", l))
    
    res <- rbind(cbind(y_train ,tr_pred_labels, rep("train", length(tr_predictions))),
          cbind(y_test, te_pred_labels, rep("test", length(te_predictions))))
    all_predictions <- rbind(all_predictions, res)
    
  }
  
  train_predictions <- all_predictions[all_predictions[,3] == "train",]
  test_predictions <- all_predictions[all_predictions[,3] == "test",]
  
  te_accuracy <- round(sum(test_predictions[,1] == test_predictions[,2]) / nrow(test_predictions), 2)
  tr_accuracy <- round(sum(train_predictions[,1] == train_predictions[,2]) / nrow(train_predictions), 2)
  
  c(tr_accuracy, te_accuracy)
  
}


# list of TF names
tf_names <- read.table("encode_tf_list", header=F, stringsAsFactors = F)

# SHAPE TF specific TF specific
for(x in 1:nrow(tf_names)){
  tf_seq <- read.table(paste("random_shuffle_", tf_names[x,],sep=""),
                       header = F, stringsAsFactors = F)
  file_name <- paste("random_shuffle_", tf_names[x,], "_shape_alphabet", sep="")
  
  # converted_r <- assign_clusters(tf_seq, paste(tf_names[x,]), paste(file_name))
  converted_r <- read.table(paste(file_name), header = T, stringsAsFactors = F)
  converted_bs <- read.table(paste("tfbs_", tf_names[x,],"_shape_alphabet",sep=""),
                             header = T, stringsAsFactors = F)
  for(m in 4:21){
    mod_res <-c ()
    cluster_names<-clusters[1:num_clust[m]]
    # create TF specific data
    tf_data <- converted_bs[,m+2]
    r_data <- converted_r[,m+2]
    
    tf_seq_list <- c("tf_data", "r_data")
    TF_labels <- c(1,0)
    
    merged_data <- create_merged_data(tf_seq_list, TF_labels, cluster_names, 20)
    
    mod_res <- seq_model_parts(merged_data, 5)
    mod_res <- t(as.data.frame(c(paste(tf_names[x,]), all[m], mod_res)))
    print(mod_res)
    write.table(mod_res,"TFspecTr_TFspecTe_shape_cv5_2",row.names = F,
                col.names = F,quote=F, sep="\t", append = T)
  }
  
}


# make a model on pooled shape data
pooled_converted_r <- read.table("pooled_random_shuffle_shape_alphabet", 
                                 header = T, stringsAsFactors = F)
pooled_converted_bs <- read.table("pooled_tfbs_shape_alphabet", 
                                  header=F, stringsAsFactors = F)
colnames(pooled_converted_bs) <- colnames(pooled_converted_r)

for(m in 4:21){
  cluster_names<-clusters[1:num_clust[m]]
  
  # create pooled data
  pooled_tf_data <- pooled_converted_bs[,m+2]
  pooled_r_data <- pooled_converted_r[,m+2]
  
  seq_list <- c("pooled_tf_data", "pooled_r_data")
  TF_labels <- c(1,0)
  
  pooled_data <- create_merged_data(seq_list, TF_labels, cluster_names, 20)
  
  for(x in 1:nrow(tf_names)){
    pooled_spec_res <- c()
    converted_r <- read.table(paste("random_shuffle_", tf_names[x,], "_shape_alphabet", sep=""),
                              header = T, stringsAsFactors = F)
    converted_bs <- read.table(paste("tfbs_", tf_names[x,],"_shape_alphabet",sep=""),
                               header = T, stringsAsFactors = F)
    # create TF specific data
    tf_data <- converted_bs[,m+2]
    r_data <- converted_r[,m+2]
    
    tf_seq_list <- c("tf_data", "r_data")
    TF_labels <- c(1,0)
    
    tf_shape_data <- create_merged_data(tf_seq_list, TF_labels, cluster_names, 20)
    
    pooled_spec_res <- pooled_cv_model(pooled_data, tf_shape_data, 5)
    pooled_spec_res <- t(as.data.frame(c(tf_names[x,], all[m], pooled_spec_res)))
    
    write.table(pooled_spec_res,"pooledTr_specTe_shape_cv5_1",row.names = F,
                col.names = F,quote=F, sep="\t", append = T)
  }
  
}


# SHAPE POOLED POOLED
for(m in 4:21){
  cluster_names<-clusters[1:num_clust[m]]
  
  # create pooled data
  pooled_tf_data <- pooled_converted_bs[,m+2]
  pooled_r_data <- pooled_converted_r[,m+2]
  
  seq_list <- c("pooled_tf_data", "pooled_r_data")
  TF_labels <- c(1,0)
  
  pooled_data <- create_merged_data(seq_list, TF_labels, cluster_names, 20)
  
  pooled_spec_res <- seq_model_parts(pooled_data, 5)
  pooled_spec_res <- t(as.data.frame(c("pooled_shape", all[m], pooled_spec_res)))
  
  write.table(pooled_spec_res,"pooledTr_pooledTe_shape_cv5_1",row.names = F,
              col.names = F,quote=F, sep="\t", append = T)
  
}








# fasta file file pooled fasta sequences
tfbs_fasta <- read.table("tf_5k_merge.fasta", header=F)
tfbs_seq <- tfbs_fasta[seq(2, nrow(tfbs_fasta), 2),]

# fasta file for pooled background sequences
# backg_fasta <- read.table("background_5k_merge.fasta", header=F)
back_seq <- read.table("pooled_shuffled", header=F, stringsAsFactors = F)

nuc <- c("A", "C", "G", "T")
seq_list <- c("tfbs_seq", "back_seq")
tf_labels <- c(1, 0)

# encode both datasets and prepare pooled set of tfbs and 
# random background for training
pooled_merged_data <- create_merged_data_seq(seq_list, tf_labels)
mod_res <- seq_model_parts(pooled_merged_data, 5)
mod_res <- c(paste(tf_names[x,]), "sequence_pooled", mod_res)







# SEQUENCE pooled TF specific
for(x in 1:nrow(tf_names)){
  r_seq <- read.table(paste("random_shuffle_", tf_names[x,],sep=""),
                      header = F, stringsAsFactors = F)
  tf_fasta <- read.table(paste(tf_names[x,],"5k.fasta",sep=""),
                         header = F, stringsAsFactors = F)
  tf_seq <- as.data.frame(tf_fasta[seq(2, nrow(tf_fasta), 2), ])
  
  tf_seq_list <- c("tf_seq", "r_seq")
  TF_labels <- c(1,0)
  
  tf_merged_data <- create_merged_data_seq(tf_seq_list, TF_labels)
  
  # model accuracy
  acc_res <- pooled_cv_model(pooled_merged_data, tf_merged_data, 5)
  mod_res <- t(as.data.frame(c(paste(tf_names[x,]), "sequence", acc_res)))
  write.table(mod_res,"pooledTr_TFspecTe_sequence_cv5_1",row.names = F,
              col.names = F,quote=F, sep="\t", append = T)
  
}


# SEQUENCE TF specific TF specific
for(x in 1:nrow(tf_names)){
  r_seq <- read.table(paste("random_shuffle_", tf_names[x,],sep=""),
                       header = F, stringsAsFactors = F)
  tf_fasta <- read.table(paste(tf_names[x,],"5k.fasta",sep=""),
                       header = F, stringsAsFactors = F)
  tf_seq <- as.data.frame(tf_fasta[seq(2, nrow(tf_fasta), 2), ])
  
  tf_seq_list <- c("tf_seq", "r_seq")
  TF_labels <- c(1,0)
  
  merged_data <- create_merged_data_seq(tf_seq_list, TF_labels)
  
  mod_res <- seq_model_parts(merged_data, 5)
  mod_res <- t(as.data.frame(c(paste(tf_names[x,]),  mod_res)))
  write.table(mod_res,"specTr_TFspecTe_sequence_cv5_1",row.names = F,
              col.names = F,quote=F, sep="\t", append = T)
  
}















