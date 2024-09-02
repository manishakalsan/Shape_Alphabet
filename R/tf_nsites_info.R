setwd("/home/manisha/manisha-extra/shape_states")

# list of TF names
tf_names <- read.table("encode_tf_list", header=F, stringsAsFactors = F)

options(scipen=999)

tf_info <- lapply(1:nrow(tf_names), function(x){
     nseq <- system(paste("wc -l ", tf_names[x,], "_tfbs.fasta", sep=""), intern = T)
     nseq <- as.numeric(strsplit(nseq," ")[[1]][1]) / 2
     
     c(tf_names[x,], nseq)
     
})
tf_info <- do.call(rbind, tf_info)
write.table(tf_info, "tf_nsites_info", sep="\t", quote=F, row.names = F, col.names = F)

# calculate minimum number of sites for pooled datasets
min_sites <- tf_info[which(tf_info[,2] == min(tf_info[,2])),]

