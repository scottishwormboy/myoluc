#myoluc_genepop.R
#below is based on 6 populations, edited from previous

setwd("~stevep11/Datafiles/Bats/myoluc/pops/counts")
raw.counts <- list()
for(i in 1:6) raw.counts[[i]] <- read.table(paste("sample",i,".counts.gz",sep=""),stringsAsFactors=FALSE,col.names=c("chr","pos","seq","depth","pileup","ref","alt"),sep="\t",comment.char="",quote="")[,c(1,2,4,6,7)]

sapply(raw.counts,nrow) #will be slightly different for each
for(i in 1:6) raw.counts[[i]][,"snp.id"] <- paste(raw.counts[[i]][,"chr"],raw.counts[[i]][,"pos"],sep="_")
for(i in 1:6) raw.counts[[i]][,"vary.cnt"] <- raw.counts[[i]][,"depth"] - rowSums(raw.counts[[i]][,c("ref","alt")])

#filter.counts <- matrix(TRUE,nrow=nrow(raw.counts[[1]]),ncol=length(raw.counts))
tmp <- character(0)
for(i in 1:6) tmp<-c(tmp,raw.counts[[i]]$snp.id)
big.counts <- data.frame(snp.id=unique(tmp),stringsAsFactors=FALSE)
for(i in 1:6){
  big.counts <- merge(big.counts,raw.counts[[i]][,c("snp.id","ref","alt","vary.cnt")],all.x=TRUE)
  names(big.counts)[names(big.counts) %in% c("ref","alt","vary.cnt")] <- paste(c("ref","alt","vary.cnt"),i,sep="_")
}

tst <- big.counts[1:10000,]
tst$all.zero <- apply(tst[,grep("vary",names(tst))] ,1, function(X){all(X==0)})


big.counts$all.zero <- apply(big.counts[,grep("vary",names(big.counts))] ,1, function(X){all(X==0)})
#table(big.counts$all.zero) #c. 1.1m FALSE 21m TRUE
tmp.alt <- rowSums(big.counts[,grep("alt",names(big.counts))])
tmp.ref <- rowSums(big.counts[,grep("ref",names(big.counts))])
big.counts$pop.freq <- tmp.alt/(tmp.alt+tmp.ref)
big.counts$freq.filt <- big.counts$pop.freq>0.05 & big.counts$pop.freq<0.95
big.counts$total.depth <- tmp.alt + tmp.ref

tmp.quan <- quantile(big.counts$total.depth,c(0.1,0.9),na.rm=TRUE)

big.counts$depth.filt <- big.counts$total.depth>=tmp.quan[1] & big.counts$total.depth<tmp.quan[2] & !is.na(big.counts$total.depth)

filt.counts <- big.counts[big.counts$all.zero & big.counts$freq.filt & big.counts$depth.filt,grep("vary",names(big.counts),invert=TRUE)] #13.5m SNPs
write.table(filt.counts,"filt.counts.w6.tsv",quote=FALSE,row.names=FALSE,sep="\t")
save(big.counts,filt.counts,file="counts.w6.Rdata")
rm(raw.counts)

#Thomas: Check I have the samples correct#
sample.map <- data.frame(no=1:6,site=c("UPMI","UPMI","PA","NY","NY","PA"),time=c("dur","pre","post","post","pre","pre"),selection=c(0.5,0,1,1,0,0))
sample.map$no.fac <- factor(paste(sample.map$site,sample.map$no,sep="_"))

#read in glm function edited from Fasciola analysis
source("~stevep11/Datafiles/Bats/labbook/bat_glm.R")

tst.snp <- head(filt.counts,40)
test.list <- apply(tst.snp,1,snpGLM,m.sample=sample.map,m.formula=~selection)
names(test.list) <- tst.snp[,"snp.id"]
#test.df <- cbind(tst.snp[,1],t(sapply(test.list,function(X){as.numeric(X$coef)})),sapply(test.list,function(X){as.numeric(X$LRT)}),sapply(test.list,function(X){as.numeric(X$pval)}))
test.mat <- cbind(t(sapply(test.list,function(X){as.numeric(X$coef)})),sapply(test.list,function(X){as.numeric(X$LRT)}),sapply(test.list,function(X){as.numeric(X$pval)}))
colnames(test.mat) <- c("Intercept","selection","LRT","pval")
rownames(test.mat) <- tst.snp[,"snp.id"]

#too many high LRT/low pvalues. Need to produce a randomisation
lrt.null <- numeric(0)
sel.null <- numeric(0)

for(i in 1:100){
  tst.snp <- filt.counts[seq(1+i,nrow(filt.counts),length=1000),]
  rand.map <- sample.map[sample(1:6),]
  test.list <- apply(tst.snp,1,snpGLM,m.sample=rand.map,m.formula=~selection)
  lrt.null <- c(lrt.null,sapply(test.list,function(X){as.numeric(X$LRT)}))
  sel.null <- c(sel.null,sapply(test.list,function(X){as.numeric(X$coef)[2]}))
}

wgs.sel.list <- apply(filt.counts,1,snpGLM,m.sample=sample.map,m.formula=~selection)
names(wgs.sel.list) <- filt.counts[,"snp.id"]
wgs.sel.mat <- cbind(t(sapply(wgs.sel.list,function(X){as.numeric(X$coef)})),sapply(wgs.sel.list,function(X){as.numeric(X$LRT)}),sapply(wgs.sel.list,function(X){as.numeric(X$pval)}))
colnames(wgs.sel.mat) <- c("Intercept","selection","LRT","pval")
save(wgs.sel.mat,file="wgs.sel.mat.w6.Rdata")

wgs.sel.df <- as.data.frame(wgs.sel.mat)
wgs.sel.df$chr <- sub("_.*","",rownames(wgs.sel.df))
wgs.sel.df$pos <- as.numeric(sub(".*_","",rownames(wgs.sel.df)))

tmp1 <- wgs.sel.df[grep("^G",wgs.sel.df$chr),]
tmp2 <- wgs.sel.df[grep("^G",wgs.sel.df$chr,invert=TRUE),]
wgs.srt.sel.df <- rbind(tmp1[order(tmp1$chr,tmp1$pos),],tmp2[order(tmp2$chr,tmp2$pos),])
save(wgs.srt.sel.df,file="wgs.srt.sel.w6.df.Rdata")
rm(tmp1,tmp2)
gc()

#write functions for heterozygosity and Fst
#from p141 Weir
calcHety <- function(ref, alt){
  if(any(is.na(c(ref,alt)))) return(list(H=NA,VarH=NA))
  n.haps <- sum(c(ref,alt))
  if(n.haps<=0) return(list(H=NA,VarH=NA))
  n <-  n.haps/2
  p <- ref/n.haps
  H <- 2*p*(1-p)
  list(H=H,VarH=(H*(1-H))/n)
}

#from p166 Weir
calcFst <- function(ref,alt){
  if(any(is.na(c(ref,alt)))) return(Fst=NA)
  if(length(ref)!=length(alt)) return(Fst=NA)
  if(length(ref)<2 | length(alt)<2) return(Fst=NA)
  r <- length(ref)
  p <- alt/(ref+alt)
  #p.hat <- sum(alt)/(sum(ref)+sum(alt))
  p.hat <- mean(p)
  n <- (alt+ref)
  #Fst <- sum((p-p.hat)^2)/((r-1)*p.hat*(1-p.hat))
  #Fst <- sum((p-p.hat)^2)/(r*p.hat*(1-p.hat)) #note r not r-1
  pt.a <- sum(n*(p-p.hat)^2)/((r-1)*mean(n))
  pt.b <- p.hat*(1-p.hat)
  Fst <- pt.a/pt.b
  Fst
}

#try again, this works, but doesn't weight by sampling effort
calcFst2 <- function(ref,alt){
  if(any(is.na(c(ref,alt)))) return(Fst=NA)
  if(length(ref)!=length(alt)) return(Fst=-100)
  if(length(ref)<2 | length(alt)<2) return(Fst=100)
  r <- length(ref)
  p <- alt/(ref+alt)
  H <- 2*p*(1-p)
  Hs <- mean(H)
  Ht <- 2*mean(p)*(1-mean(p))
  Fst <- (Ht-Hs)/Ht
  Fst
}

Het.list <- apply(filt.counts[,-1],1,function(X){calcHety(X[grep('ref',names(X))],X[grep('alt',names(X))])})
Fst.list <- apply(filt.counts[,-1],1,function(X){calcFst2(X[grep('ref',names(X))],X[grep('alt',names(X))])})

popgen.df <- data.frame(filt.counts$snp.id,chr=sub("_.*","",filt.counts$snp.id),pos=as.numeric(sub(".*_","",filt.counts$snp.id)),Fst=Fst.list,stringsAsFactors=FALSE)
popgen.df[,c("H1","H2","H3","H4","H5","H6")] <- t(sapply(Het.list,function(X){X$H}))
popgen.df[,c("Hse1","Hse2","Hse3","Hse4","Hse5","Hse6")] <- t(sapply(Het.list,function(X){sqrt(X$VarH)})) #check this is se not sd, I think its the variance of the estimate

tmp1 <- popgen.df[grep("^G",popgen.df$chr),]
tmp2 <- popgen.df[grep("^G",popgen.df$chr,invert=TRUE),]
popgen.srt.df <- rbind(tmp1[order(tmp1$chr,tmp1$pos),],tmp2[order(tmp2$chr,tmp2$pos),])
rm(tmp1,tmp2)
save(popgen.srt.df,file="popgen.srt.w6.df.Rdata")

popgen.quant.fst <- quantile(popgen.srt.df$Fst,c(0.5,0.75,0.9,0.95,0.99,0.999))
#GL429776 looks to have a lot of high Fst values
#check, BTG1 possible candidate, GL429776:14598505:14600859, or is this too far away?
#KITLG, GL429776:11131621:11179955
#from fst region is 9.6 - 13.2Mb
#http://string-db.org/cgi/network.pl?taskId=goYFPEFOwXbs

#also backed up by glm
quantile(wgs.sel.df$LRT[wgs.sel.df$chr == "GL429776"],c(0.5,0.75,0.9,0.95,0.99,0.999))
quantile(wgs.sel.df$LRT,c(0.5,0.75,0.9,0.95,0.99,0.999))

#provide a moving window function. eg median over 100 SNPs, move 10 SNPs at a time
#take a chromosome and return chr, mid position, median value
movingWindowCalc <- function(in.df,chr="chr",pos="pos",jump=10,sample=100,value){
  out.df <- data.frame(chr=character(0),pos=numeric(0),value=numeric(0))
  if(in.df[1,chr] != in.df[nrow(in.df),chr]) warning("different chromosomes into moving window")
  if(nrow(in.df) <= sample){
    #window smaller than sample
    tmp.out <- data.frame(chr=in.df[1,chr],pos=mean(range(in.df[,pos])),value=median(in.df[,value]),stringsAsFactors=FALSE)
    out.df <- tmp.out 
    }else{
    for(i in seq(1,nrow(in.df)-sample,jump)){
      tmp <- in.df[seq(i,i+sample-1),c(chr,pos,value)]
      tmp.out <- data.frame(chr=tmp[1,chr],pos=mean(range(tmp[,pos])),value=median(tmp[,value],na.rm=TRUE),stringsAsFactors=FALSE)
      out.df <- rbind(out.df,tmp.out)
      }
   }
  names(out.df) <- c(chr,pos,value)
  return(out.df)
}

tst <- movingWindowCalc(popgen.srt.df[1:1000,],value="Fst")
tst <- movingWindowCalc(popgen.srt.df[popgen.srt.df$chr=="GL429776",],value="Fst",jump=100,sample=1000)
#might be a 5Mb region, 50 genes?

fst.1000.window <- data.frame(chr=character(0),pos=numeric(0),Fst=numeric(0))
for(scaf in unique(popgen.srt.df$chr)){
  fst.1000.window <- rbind(fst.1000.window,movingWindowCalc(popgen.srt.df[popgen.srt.df$chr %in% scaf,],value="Fst",jump=100,sample=1000))
}

ct <- 0
lrt.1000.window <- data.frame(chr=character(0),pos=numeric(0),LRT=numeric(0))
for(scaf in unique(wgs.srt.sel.df$chr)){
  lrt.1000.window <- rbind(lrt.1000.window,movingWindowCalc(wgs.srt.sel.df[wgs.srt.sel.df$chr %in% scaf,],value="LRT",jump=100,sample=1000))
  ct <- ct+1
  if(ct %% 100 == 0) cat(scaf,"\n")
}
save(fst.1000.window,file="fst.w6.1000.windows.Rdata")
save(lrt.1000.window,file="lrt.w6.1000.windows.Rdata")
#could also generate 5%, 25%, 50%, 75%, 95%, for a pretty picture?

#list of number of SNPs
snps.in.chr.by <- by(popgen.srt.df,popgen.srt.df$chr,nrow)
save(snps.in.chr.by,file="snps.in.chr.w6.by.Rdata")

#pair-wise comparison of populations NY pre/post UPMI pre/dur NY pre/UPMI pre

#  no site time selection no.fac
#  1 UPMI  dur       0.5 UPMI_1
#  2 UPMI  pre       0.0 UPMI_2
#  3   PA post       1.0   PA_3
#  4   NY post       1.0   NY_4
#  5   NY  pre       0.0   NY_5
#  6   PA  pre       0.0   PA_6

Fst.list.NY <- apply(filt.counts[,-1],1,function(X){calcFst2(X[grep('ref',names(X))[c(4,5)]],X[grep('alt',names(X))[c(4,5)]])})
Fst.list.UPMI <- apply(filt.counts[,-1],1,function(X){calcFst2(X[grep('ref',names(X))[c(1,2)]],X[grep('alt',names(X))[c(1,2)]])})
Fst.list.NY_UPMI <- apply(filt.counts[,-1],1,function(X){calcFst2(X[grep('ref',names(X))[c(2,5)]],X[grep('alt',names(X))[c(2,5)]])})
Fst.list.NY_PA <- apply(filt.counts[,-1],1,function(X){calcFst2(X[grep('ref',names(X))[c(3,4)]],X[grep('alt',names(X))[c(3,4)]])})
Fst.list.PA <- apply(filt.counts[,-1],1,function(X){calcFst2(X[grep('ref',names(X))[c(3,6)]],X[grep('alt',names(X))[c(3,6)]])})


Fst.pairs <- data.frame(filt.counts$snp.id,chr=sub("_.*","",filt.counts$snp.id),pos=as.numeric(sub(".*_","",filt.counts$snp.id)),NY=Fst.list.NY,UPMI=Fst.list.UPMI,NY_UPMI=Fst.list.NY_UPMI,NY_PA=Fst.list.NY_PA,PA=Fst.list.PA,stringsAsFactors=FALSE)
tmp1 <- Fst.pairs[grep("^G",Fst.pairs$chr),]
tmp2 <- Fst.pairs[grep("^G",Fst.pairs$chr,invert=TRUE),]
Fst.srt.pairs <- rbind(tmp1[order(tmp1$chr,tmp1$pos),],tmp2[order(tmp2$chr,tmp2$pos),])

ct <- 0
fst.pr.1000.window <- data.frame(chr=character(0),pos=numeric(0),NY=numeric(0),UPMI=numeric(0),NY_UPMI=numeric(0))
for(scaf in unique(Fst.srt.pairs$chr)){
  ct <- ct+1
  if(ct %% 100 == 0) cat(scaf,"\n")
  tmp.NY <- movingWindowCalc(Fst.srt.pairs[Fst.srt.pairs$chr %in% scaf,],value="NY",jump=100,sample=1000)
  tmp.UPMI <- movingWindowCalc(Fst.srt.pairs[Fst.srt.pairs$chr %in% scaf,],value="UPMI",jump=100,sample=1000)
  tmp.NY_UPMI <- movingWindowCalc(Fst.srt.pairs[Fst.srt.pairs$chr %in% scaf,],value="NY_UPMI",jump=100,sample=1000)
  tmp.NY_PA <- movingWindowCalc(Fst.srt.pairs[Fst.srt.pairs$chr %in% scaf,],value="NY_PA",jump=100,sample=1000)
  tmp.PA <- movingWindowCalc(Fst.srt.pairs[Fst.srt.pairs$chr %in% scaf,],value="PA",jump=100,sample=1000)
  tmp <- cbind(tmp.NY,tmp.UPMI[,"UPMI"],tmp.NY_UPMI[,"NY_UPMI"],tmp.NY_PA[,"NY_PA"],tmp.PA[,"PA"])
  fst.pr.1000.window <- rbind(fst.pr.1000.window,tmp)  
}
#should have got scaf from the sorted version, remember to sort out
names(fst.pr.1000.window) <- c("chr","pos","NY","UPMI","NY_UPMI","NY_PA","PA")
save(fst.pr.1000.window,file="fst.pr.1000.w6.window.Rdata")


movingWindowCalc2 <- function(in.df,chr="chr",pos="pos",jump=10,sample=100,value){
  #for hety, colMeans
  out.df <- in.df[0,c(chr,pos,value)]
  #out.df[,value] <- numeric(0)
  if(in.df[1,chr] != in.df[nrow(in.df),chr]) warning("different chromosomes into moving window")
  if(nrow(in.df) <= sample){
    #window smaller than sample
    tmp.out <- data.frame(chr=in.df[1,chr],pos=mean(range(in.df[,pos])),stringsAsFactors=FALSE)
    tmp.out[,value] <- colMeans(in.df[,value],na.rm=TRUE)
    out.df <- tmp.out 
    }else{
    for(i in seq(1,nrow(in.df)-sample,jump)){
      tmp <- in.df[seq(i,i+sample-1),c(chr,pos,value)] #in for loop
      tmp.out <- data.frame(chr=tmp[1,chr],pos=mean(range(tmp[,pos])),stringsAsFactors=FALSE)
      tmp.out[,value] <- colMeans(tmp[,value],na.rm=TRUE)
      out.df <- rbind(out.df,tmp.out)
      }
   }
  names(out.df) <- c(chr,pos,value)
  return(out.df)
}


het.1000.window <- data.frame(chr=character(0),pos=numeric(0),H1=numeric(0),H2=numeric(0),H3=numeric(0),H4=numeric(0),H5=numeric(0),H6=numeric(0))
ct <- 0
for(scaf in unique(popgen.srt.df$chr)){
  het.1000.window <- rbind(het.1000.window,movingWindowCalc2(popgen.srt.df[popgen.srt.df$chr %in% scaf,],value=c("H1","H2","H3","H4","H5","H6"),jump=100,sample=1000))
  ct <- ct+1
  if(ct %% 100 == 0) cat(scaf,"\n")
}
save(het.1000.window,file="het.1000.windows.w6.Rdata")

#extract data for GL429776
GL429776.cnt <- filt.counts[sub("_.*","",filt.counts$snp.id) == "GL429776",]
GL429776.cnt$chr <- sub("_.*","",GL429776.cnt$snp.id)
GL429776.cnt$pos <- as.numeric(sub(".*_","",GL429776.cnt$snp.id))
GL429776.cnt <- GL429776.cnt[order(GL429776.cnt$pos),]
save(GL429776.cnt,file="GL429776.cnt.Rdata")

#Added 8/4/18 get files ready for SelEstim

#filt.counts.w6 <- read.table(filt.counts,"filt.counts.w6.tsv",sep="\t",header=TRUE)
selest <- filt.counts[,grep("ref_|alt_",names(filt.counts))]
source("/pub34/stevep11/bin/SelEstim_1.1.7/R/SelEstim.R")

#randomly swap reference allele
tmp.s <- sample(c(0,1),nrow(selest),replace=TRUE)
selest[tmp.s==1,] <- selest[tmp.s==1,c(t(matrix(c(seq(2,ncol(selest),2),seq(1,ncol(selest),2)),ncol(selest)/2,2)))] #probably an easier way. Check if reusing code.
