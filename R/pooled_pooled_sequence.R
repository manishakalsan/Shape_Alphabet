
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
seq_tf_test_data <- function(tf_id, sequ_len){
  nuc <- c("A", "C", "G", "T")
  
  test_tfbs_data<- read.table(paste(tf_id, "5k.fasta", sep=""), header = F)
  test_backg_data <- read.table(paste("background_", tf_id, ".fasta", sep=""), header = F)
  
  # extract clustering method specific columns
  tfbs_df <- as.data.frame(test_tfbs_data[seq(2, nrow(test_tfbs_data), 2), ])
  back_df <- as.data.frame(test_backg_data[seq(2, nrow(test_backg_data), 2), ])
  
  # one hot encoding
  data_names <- c("tfbs_df", "back_df")
  data_labels <- c(1,0)
  
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
      one_seq_code<-oneHot(cl_data[loci,])
      c(one_seq_code)
    }, mc.cores=20)
    coded_df <- do.call(rbind, coded)
    coded_df <- cbind(coded_df, rep(data_labels[dn], nrow(coded_df)))
    coded_df
  })
  tf_backg_df <- do.call(rbind, coded_df)
  
  coln<-c()
  for(one_pos in 1:as.numeric(nchar(tfbs_df[1,1]))){
    for(one_letter in 1:length(nuc)){
      coln<-c(coln,paste("pos",one_pos,nuc[one_letter],sep="_"))
    }
  }
  colnames(tf_backg_df)<-c(coln,"TFlabel")
  tf_backg_df
}
seq_reg_pooled <- function(merged_data, test_data, prefix){
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
    learning_rate = 0.1
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
  roc <- roc(y_test, predictions)
  auc_roc <- round(auc(roc),2)
  
  results<-as.data.frame(c(prefix, rmse, mae, accuracy, auc_roc))
  t(results)
}

# fasta file file pooled fasta sequences
tfbs_fasta <- read.table("tf_5k_merge.fasta", header=F)
tfbs_seq <- tfbs_fasta[seq(2, nrow(tfbs_fasta), 2),]

# fasta file for pooled background sequences
backg_fasta <- read.table("background_5k_merge.fasta", header=F)
back_seq <- backg_fasta[seq(2, nrow(backg_fasta), 2),]

nuc <- c("A", "C", "G", "T")
seq_list <- c("tfbs_seq", "back_seq")
tf_labels <- c(1, 0)

# encode both datasets and prepare pooled set of tfbs and 
# random background for training
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
    nuc_coded <- oneHot(seq_data[x,])
    c(nuc_coded)
  }, mc.cores = 20)
  encoded <- do.call(rbind, encoded)
  encoded <- cbind(encoded, rep(tf_labels[x], nrow(encoded)))
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

results<-c("shape",rmse, mae, log_loss, accuracy, auc_roc)
t(results)


