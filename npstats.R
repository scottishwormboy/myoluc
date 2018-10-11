setwd("~/Dropbox/Datafiles/Bats/myoluc/npstats")
#this path may change
#read in 1000 x 1000SNPs files
samp1.1000 <- read.table("rand_stats/sample1.1m.1000.stats",sep="\t",header=TRUE)
samp2.1000 <- read.table("rand_stats/sample2.1m.1000.stats",sep="\t",header=TRUE)
samp3.1000 <- read.table("rand_stats/sample3.1m.1000.stats",sep="\t",header=TRUE)
samp4.1000 <- read.table("rand_stats/sample4.1m.1000.stats",sep="\t",header=TRUE)
samp5.1000 <- read.table("rand_stats/sample5.1m.1000.stats",sep="\t",header=TRUE)
samp6.1000 <- read.table("rand_stats/sample6.1m.1000.stats",sep="\t",header=TRUE)

samp.quants <- rbind(quantile(samp1.1000$Watterson,c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),
quantile(samp2.1000$Watterson,c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),
quantile(samp3.1000$Watterson,c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),
quantile(samp4.1000$Watterson,c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),
quantile(samp5.1000$Watterson,c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)),
quantile(samp6.1000$Watterson,c(0.01,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99))
)



