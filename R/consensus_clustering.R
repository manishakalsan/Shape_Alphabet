library(ConsensusClusterPlus)

setwd("/media/manisha/baa3e1ca-d715-461f-a686-43a6d4a8a35e/home/manisha/dynaseq_data")
dyn_data <- readRDS("dictionary_ensemble_5bin_5mer.rds")

med_data <- apply(dyn_data,1,mad)
d = dyn_data[rev(order(med_data)),]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results = ConsensusClusterPlus(t(d), maxK=20, reps=100, pItem=0.8, pFeature=1,
                               title="pam_euclidean_shape", clusterAlg="pam", 
                               distance="euclidean", writeTable = TRUE,
                               seed=1262118388.71279, plot="pdf")
resICL = calcICL(results,title="pam_euclidean_shape", plot="pdf")


