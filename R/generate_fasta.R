# use the tf and bed id file to extract summit positions 
# from each bed file and use it to create fasta files for each TF

setwd("/home/manisha/Desktop/shape_clusters/hepg2/raw_bed/")
tf_id <- read.table("raw_bed/tf_bed_ids", header=T)

options(scipen = 999)

flank_length <- 25

for(i in 1:nrow(tf_id)){
     # read narrowpeak format file
     bed_file <- read.table(paste(tf_id[i,2], "bed", sep="."), header=F)
     
     # prepare coordinates
     chr_num <- bed_file[,1]
     summit_pos <- bed_file[,2] + bed_file[,10]
     start_pos <- as.numeric(summit_pos - flank_length)
     end_pos <- as.numeric(summit_pos + flank_length + 1)
     seq_id <- c(paste(as.character(tf_id[i, 1]), 1:length(start_pos), sep="_"))
     
     new_bed <- data.frame(
          chr = as.character(chr_num),
          start = as.numeric(start_pos),
          end = as.numeric(end_pos),
          seq_id = seq_id,
          stringsAsFactors = FALSE  
     )
     
     write.table(new_bed, file=paste(tf_id[i,1],"_tfbs.bed",sep=""), row.names = F, col.names = F,quote=F,sep="\t")
     system(paste("bedtools getfasta -fi hg38.fa "," -bed ",tf_id[i,1],"_tfbs.bed"," -fo ",tf_id[i,1],"_tfbs.fasta", sep=""))
     
}
