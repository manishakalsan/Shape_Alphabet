
# take the filtered data table and use the files associated to each TF 
# and select one file for each TF based on its presence in the file 
# containing the list of datasets selected based on filetype criterion from encode

setwd("~/Desktop/shape_clusters/hepg2/")
dataset_list <- readxl::read_xlsx("Desktop/shape_clusters/hepg2/encode_filtered.xlsx")
selected_criteria_files <- read.table("Desktop/shape_clusters/hepg2/encode_datasets.txt", skip=1)

library(dplyr)

# create one row for each bed file of each TF
selected_files <- lapply(1:nrow(dataset_list), function(x){
     # find list of all files associated with an entry based on a dataset: 
     # this could be redundant in terms of a TF as one TF can have 
     # multiple datasets or studies
     file_names <- unlist(strsplit(gsub("/files/|/", "", dataset_list$Files[x]), ","))
     
     # all files with the selection criteria with the file names
     selected_datasets_names <- as.data.frame(selected_criteria_files[grep(paste(file_names, collapse="|"), 
                                                                           selected_criteria_files[,1]),1])
     
     new_df <- dataset_list %>%
          slice(rep(x, nrow(selected_datasets_names)))
     new_df <- cbind(new_df, selected_datasets_names)
     
     colnames(new_df) <- c(colnames(new_df)[-ncol(new_df)], "bed_file")
     
     new_df
     
})
selected_files <- do.call(rbind, selected_files)

# randomly select one bed file per TF
set.seed(173682)

list_tf <- unique(selected_files$`Target of assay`)

tf_bed <- lapply(1:length(list_tf), function(x){
     bed_df <- selected_files[which(selected_files$`Target of assay` == list_tf[x]),]
     
     # randomly select one entry per TF
     new_df <- bed_df[sample(1:nrow(bed_df), 1),]
     
     new_df
})
tf_bed <- do.call(rbind, tf_bed)
write.table(tf_bed, "selected_bed", row.names = F, sep="\t", quote=F)
































