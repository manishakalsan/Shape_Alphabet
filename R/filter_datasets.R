
# take the filtered data table and use the files associated to each TF 
# and select one file for each TF based on its presence in the file 
# containing the list of datasets selected based on filetype criterion from encode

setwd("~/Desktop/shape_clusters/")

file_list <- system("ls *experiment*", intern=T)
file_list <- file_list[-grep("tsv", file_list)]

# dataset_list <- readxl::read_xlsx("Desktop/shape_clusters/hepg2/encode_filtered.xlsx")
# selected_criteria_files <- read.table("Desktop/shape_clusters/hepg2/encode_datasets.txt", skip=1)

library(dplyr)

# for each cell line data, we need to use the tsv file to select the datasets 
# in selected criteria file to generate the list of files (from selected criteria) 
# for each TF and use it to download files from encode database

all_tf_num <- c()
for(i in 1:length(file_list)){
     cell_line <- strsplit(file_list[i], "_")[[1]][1]
     
     dataset_list <- read.table(file_list[i], sep="\t", header=T, skip=1)
     selected_criteria_files <- read.table(paste(cell_line, "download_file_list.txt", sep="_"), skip=1)
     
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
     
     list_tf <- unique(selected_files$`Target.of.assay`)
     
     tf_bed <- lapply(1:length(list_tf), function(x){
          bed_df <- selected_files[which(selected_files$`Target.of.assay` == list_tf[x]),]
          
          # randomly select one entry per TF
          new_df <- bed_df[sample(1:nrow(bed_df), 1),]
          
          new_df
     })
     tf_bed <- do.call(rbind, tf_bed)
     tf_num <- c(cell_line, length(list_tf))
     all_tf_num <- rbind(all_tf_num, tf_num)
     write.table(tf_bed, paste(cell_line, "selected_bed_list", sep="_"), row.names = F, sep="\t", quote=F)
     
}
colnames(all_tf_num) <- c("cell_line", "unique_TFs")
write.table(all_tf_num, "unique_cell_line_tfs", row.names = F, sep="\t", quote=F)

