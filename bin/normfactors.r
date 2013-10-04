# normfactors.r
# Tom Cooke
# 2013-09-03

library(plyr)

samples <- read.table('samples', colClasses=c('character'))[,1]
#samples <- samples[1:6]

counts <- data.frame(seq(0,980,20), seq(20,1000,20))

for (i in 1:length(samples)){
  filename <- paste(samples[i], '.rs_mapping_stats', sep='')
  hist <- read.table(filename, header=T)
  hist <- hist[which(hist$insert_size <= 1000&hist$insert_size > 0),]
  hist$bin <- cut(hist$insert_size, breaks=seq(0,1000,20))
  hist2 <- ddply(hist, .drop=F, .(bin), summarize, counts=sum(reads))
  counts <- cbind(counts, hist2$counts)
}
names(counts) <- c('lower_bin_limit', 'upper_bin_limit', samples)
mean <- rowMeans(counts[,samples])
normfactors <- counts[,samples] / mean
normfactors <- cbind(counts[,1:2], normfactors)
write.table(normfactors, 'normfactors.txt', row.names=F, col.names=T, sep='\t', quote=F)
#write.table(normfactors, 'normfactors.1-6.txt', row.names=F, col.names=T, sep='\t', quote=F)
