setwd("/home/manisha/manisha-extra/shape_states")
# setwd("~/Desktop/shape_clusters/hepg2/raw_bed/")

# list of TF names
tf_names <- read.table("encode_tf_list_500thresh", header=F, stringsAsFactors = F)
tf_names <- as.data.frame(tf_names[,1])
# prepare alphabet files for TFBS and shuffled datasets using shape and random shape

options(scipen=999)
nsites <- as.numeric(500)

for(x in 1:nrow(tf_names)){
     
     # read the tf binding site sequences and used them to generate cluster names
     # on the basis of cluster labels
     set.seed(125637)
     tf_fasta <- read.table(paste(tf_names[x,],"_tfbs.fasta",sep=""),
                            header = F, stringsAsFactors = F)
     tf_seq <- as.data.frame(tf_fasta[seq(2,nrow(tf_fasta), 2),])
     
     if(nrow(tf_seq) > nsites){
          tf_seq <- as.data.frame(tf_seq[sample(1:nrow(tf_seq), nsites),])
     }
     
     write.table(tf_seq, "pooled_tfbs_500", row.names = F, col.names = F, quote=F, append = T)
     
     r_fasta <- read.table(paste("shuffled_", tf_names[x,], ".fasta", sep=""),
                           header = F, stringsAsFactors = F)
     r_seq <- as.data.frame(r_fasta[seq(2,nrow(r_fasta), 2),])
     if(nrow(r_seq) > nsites){
          r_seq <- as.data.frame(r_seq[sample(1:nrow(r_seq), nsites),])
     }
     
     write.table(tf_seq, "pooled_shuffled_500", row.names = F, col.names = F, quote=F, append = T)
     
     
}

