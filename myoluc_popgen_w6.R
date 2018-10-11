#analysis of myoluc bats, some plots genome scan etc
setwd("~/Dropbox/Datafiles/Bats/myoluc/popgen")
myoluc.fai <- read.table("myoluc2.0.fa.fai",stringsAsFactors=FALSE,col.names=c("chr","length","V3","V4","V5"))
myoluc.fai[,3:5] <- NULL
myoluc.fai$cum.start <- c(0,cumsum(myoluc.fai$length))[-nrow(myoluc.fai)]
myoluc.fai$scaf.no <- seq(1,nrow(myoluc.fai))
myoluc.fai$bin <- 1
myoluc.fai$bin[seq(1,nrow(myoluc.fai),2)] <- 2 #this is to plot colours later on

load("fst.w6.1000.windows.Rdata")
load("lrt.w6.1000.windows.Rdata")
load("snps.in.chr.by.Rdata")
fst.1000.window$g.pos <- 0
for(scaf in unique(fst.1000.window$chr)){
  fst.1000.window$g.pos[fst.1000.window$chr %in% scaf] <- fst.1000.window$pos[fst.1000.window$chr %in% scaf] + myoluc.fai$cum.start[myoluc.fai$chr %in% scaf]
}
fst.1000.window$pass.filt <- TRUE
fst.1000.window$pass.filt[fst.1000.window$chr %in% names(snps.in.chr.by)[snps.in.chr.by<1000]] <- FALSE
fst.1000.window$bin <- 1
fst.1000.window$bin[fst.1000.window$chr %in% myoluc.fai$chr[myoluc.fai$bin==2]] <- 2
plot(Fst~g.pos,data=fst.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.12))
tmp.filt <- fst.1000.window$pass.filt & fst.1000.window$Fst>2*mean(fst.1000.window$Fst)
points(Fst~g.pos,data=fst.1000.window[tmp.filt,],col=c("darkblue","darkgreen")[fst.1000.window$bin[tmp.filt]],pch=".",cex=1.5)

#focus on GL429776
plot(Fst~pos,data=fst.1000.window[fst.1000.window$chr == "GL429776",],pch=".",ylim=c(0,0.12),cex=2)
text(0,0.1,"GL429776",pos=4)


lrt.1000.window$g.pos <- 0
for(scaf in unique(lrt.1000.window$chr)){
  lrt.1000.window$g.pos[lrt.1000.window$chr %in% scaf] <- lrt.1000.window$pos[lrt.1000.window$chr %in% scaf] + myoluc.fai$cum.start[myoluc.fai$chr %in% scaf]
}
lrt.1000.window$pass.filt <- TRUE
lrt.1000.window$pass.filt[lrt.1000.window$chr %in% names(snps.in.chr.by)[snps.in.chr.by<1000]] <- FALSE
plot(LRT~g.pos,data=lrt.1000.window[lrt.1000.window$pass.filt,],pch=".")


#pair Fst data
load("fst.pr.1000.w6.window.Rdata")
tmp1 <- fst.pr.1000.window[grep("^G",fst.pr.1000.window$chr),]
tmp2 <- fst.pr.1000.window[grep("^G",fst.pr.1000.window$chr,invert=TRUE),]
fst.pr.1000.window <- rbind(tmp1,tmp2) #sort and write over
names(fst.pr.1000.window) <- c("chr","pos","NY","UPMI","NY_UPMI","NY_PA","PA")
fst.pr.1000.window$g.pos <- 0
for(scaf in unique(fst.pr.1000.window$chr)){
  fst.pr.1000.window$g.pos[fst.pr.1000.window$chr %in% scaf] <- fst.pr.1000.window$pos[fst.pr.1000.window$chr %in% scaf] + myoluc.fai$cum.start[myoluc.fai$chr %in% scaf]
}

pdf(file="pair.fst.plot.pdf",width=12,height=4)
plot(NY~g.pos,data=fst.pr.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.12))
plot(UPMI~g.pos,data=fst.pr.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.12))
plot(PA~g.pos,data=fst.pr.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.12))
plot(NY_UPMI~g.pos,data=fst.pr.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.12))
plot(NY_PA~g.pos,data=fst.pr.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.12))
dev.off()

plot(PA~pos,data=fst.pr.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="black",ylim=c(0,0.12))
text(x=0,y=0.1,labels="GL429776",pos=4)
load("het.1000.windows.w6.Rdata")
het.1000.window$g.pos <- 0
for(scaf in unique(het.1000.window$chr)){
  het.1000.window$g.pos[het.1000.window$chr %in% scaf] <- het.1000.window$pos[het.1000.window$chr %in% scaf] + myoluc.fai$cum.start[myoluc.fai$chr %in% scaf]
}

colMeans(het.1000.window[,grep('^H',names(het.1000.window))])
#       H1        H2        H3        H4        H5        H6 
#0.2049284 0.1814257 0.2117057 0.1767890 0.1864314 0.1388724 


pdf(file="het.plot.pdf",width=12,height=4)
plot(H2~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
text(x=0,y=0.25,labels="UPMI_pre",pos=4)
plot(H1~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
text(x=0,y=0.25,labels="UPMI_dur",pos=4)
plot(H5~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
text(x=0,y=0.25,labels="NY_pre",pos=4)
plot(H4~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
text(x=0,y=0.25,labels="NY_post",pos=4)
plot(H6~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
text(x=0,y=0.25,labels="PA_pre",pos=4)
plot(H3~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
text(x=0,y=0.25,labels="PA_post",pos=4)
dev.off()


plot(H3~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))
plot(H6~g.pos,data=het.1000.window[fst.1000.window$pass.filt,],pch=".",col=c("lightblue","lightgreen")[fst.1000.window$bin[fst.1000.window$pass.filt]],ylim=c(0,0.3))

plot(H3~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="red",ylim=c(0,0.3))
points(H6~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="black")
points(H1~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="green")
points(H2~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="blue")

#PA
plot(H3-mean(het.1000.window[,"H3"])~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="red",ylim=c(-0.15,0.15))
points(H6-mean(het.1000.window[,"H6"])~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="black")

#UPMI
plot(H1-mean(het.1000.window[,"H1"])~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="green",ylim=c(-0.15,0.15))
points(H2-mean(het.1000.window[,"H2"])~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="black")

#NY
plot(H4-mean(het.1000.window[,"H4"])~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="blue",ylim=c(-0.15,0.15))
points(H5-mean(het.1000.window[,"H5"])~g.pos,data=het.1000.window[fst.1000.window$pass.filt&fst.1000.window$chr=="GL429776",],pch=".",col="black")


#PA has loss of heterozygosity in same place as Fst spike

#  no site time selection no.fac
#  1 UPMI  dur       0.5 UPMI_1
#  2 UPMI  pre       0.0 UPMI_2
#  3   PA post       1.0   PA_3
#  4   NY post       1.0   NY_4
#  5   NY  pre       0.0   NY_5
#  6   PA  pre       0.0   PA_6

#examine GL429776
load("GL429776.cnt.Rdata")
GL429776.cnt$NY_pre <- GL429776.cnt$alt_5/(GL429776.cnt$alt_5+GL429776.cnt$ref_5)
GL429776.cnt$NY_post <- GL429776.cnt$alt_4/(GL429776.cnt$alt_4+GL429776.cnt$ref_4)
GL429776.cnt$UPMI_pre <- GL429776.cnt$alt_2/(GL429776.cnt$alt_2+GL429776.cnt$ref_2)
GL429776.cnt$UPMI_dur <- GL429776.cnt$alt_1/(GL429776.cnt$alt_1+GL429776.cnt$ref_1)
GL429776.cnt$PA_pre <- GL429776.cnt$alt_6/(GL429776.cnt$alt_6+GL429776.cnt$ref_6)
GL429776.cnt$PA_post <- GL429776.cnt$alt_3/(GL429776.cnt$alt_3+GL429776.cnt$ref_3)

GL429776.diff <- GL429776.cnt[GL429776.cnt$pos>8957547 & GL429776.cnt$pos<13473506 & GL429776.cnt$pop.freq>=0.1 & GL429776.cnt$pop.freq <= 0.9,]

library(vegan)
GL429776.dist <- dist(t(GL429776.diff[,c("NY_pre","NY_post","UPMI_pre","UPMI_dur","PA_pre","PA_post")]),method="euc")
GL429776.mds.euc <- metaMDS(GL429776.dist) #warning about low stress
GL429776.mds.bray <- metaMDS(t(GL429776.diff[,c("NY_pre","NY_post","UPMI_pre","UPMI_dur","PA_pre","PA_post")]))
plot(GL429776.mds.euc,type="t")

GL429776.pc <- princomp(t(GL429776.diff[,c("NY_pre","NY_post","UPMI_pre","UPMI_dur","PA_pre","PA_post")]))

GL429776.back <- GL429776.cnt[(GL429776.cnt$pos<8957547 | GL429776.cnt$pos>13473506) & GL429776.cnt$pop.freq>=0.1 & GL429776.cnt$pop.freq <= 0.9,]
GL429776.bdist <- dist(t(GL429776.back[,c("NY_pre","NY_post","UPMI_pre","UPMI_dur","PA_pre","PA_post")]),method="euc")
GL429776.bmds.euc <- metaMDS(GL429776.bdist) 
plot(GL429776.bmds.euc,type="t")

library(ape)
plot(as.phylo(GL429776.bdist))

tst <- princomp(GL429776.diff[,c("NY_pre","NY_post","UPMI_pre","UPMI_dur","PA_pre","PA_post")])