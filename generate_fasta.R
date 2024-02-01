# tf_id <- c("304-", "305-", "3-", "461-")
tf_id <- read.table("motif_prefix_list", header=F)

options(scipen = 999)

for(i in 1:nrow(tf_id)){
  list_files <- system(paste("ls ", tf_id[i,],"*calls.all.bed", sep = ""), intern=T)
  
  bed_file<-list_files[-grep(".RC", list_files)]
  # prefix <- strsplit(bed_file, "-calls")[[1]][1]
  prefix <- tf_id[i,]
  
  bed_data<-c()
  if(file.exists(bed_file)){
    bed_data <- read.delim(bed_file, header=F) # to check the length of binding sites
  }
  
  #check length of binding sites
  bs_len <- bed_data[-1,3] - bed_data[-1,2] 
  print(bs_len[1])
  bed_data<-c()
  bs_data<-read.csv(paste(prefix,"-calls.csv",sep=""))
  head(bs_data)
  bs_data$motif <- rep(prefix, nrow(bs_data))
  
  # create data according to motif centered 20bp regions
  end <- as.data.frame(format(as.numeric(bs_data[,3] + bs_len[1] + 2), scientific=F),stringsAsFactors = F)
  start <- as.data.frame(as.numeric(bs_data[,3] - 2))
  
  names(bs_data)<-NULL
  new <- cbind(as.character(bs_data[,2]),
               as.numeric(start[,1]),
               as.numeric(end[,1]),
               as.character(bs_data[,8]))
  
  mid <- as.data.frame(round((as.numeric(new[,3]) + as.numeric(new[,2]))/2 , 0))
  new_start <- as.data.frame(mid - 10)
  new_end <- as.data.frame(mid + 10)
  new_bed <- cbind(as.character(bs_data[,2]),
                   as.numeric(new_start[,1]),
                   as.numeric(new_end[,1]),
                   as.character(bs_data[,8]))
  
  write.table(new_bed, file=paste(prefix,"_tfbs.bed",sep=""), row.names = F,col.names = F,quote=F,sep="\t")
  system(paste("bedtools getfasta -fi mm10.fa "," -bed ",prefix,"_tfbs.bed"," -fo ",prefix,"_tfbs.fasta", sep=""))
  
  
}

