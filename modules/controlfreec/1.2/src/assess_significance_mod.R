#!/usr/bin/env Rscript

library(rtracklayer)

args <- commandArgs()

dataTable <-read.table(args[5], header=TRUE);
ratio<-data.frame(dataTable)

dataTable <- read.table(args[4], header=FALSE)
cnvs<- data.frame(dataTable)

ratio$Ratio[which(ratio$Ratio==-1)]=NA

cnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))
ratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio,median_ratio=ratio$MedianRatio)

overlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)
normals <- setdiff(ratio.bed,cnvs.bed)
normals <- subsetByOverlaps(ratio.bed,normals)

#mu <- mean(score(normals),na.rm=TRUE)
#sigma<- sd(score(normals),na.rm=TRUE)

#hist(score(normals),n=500,xlim=c(0,2))
#hist(log(score(normals)),n=500,xlim=c(-1,1))

#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])
#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))
#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)

#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])
#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))
#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)

numberOfCol=length(cnvs)

for (i in c(1:length(cnvs[,1]))) {
  values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))
  #wilcox.test(values,mu=mu)
  W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)
	if(class(resultw)=="try-error") return(list("statistic"=NA,"parameter"=NA,"p.value"=NA,"null.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultw}
  KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)
	if(class(resultks)=="try-error") return(list("statistic"=NA,"p.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultks}
  #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)
  #	if(class(resultks)=="try-error") NA) else resultks
  cnvs[i,numberOfCol+1]=W(values,normals)$p.value
  cnvs[i,numberOfCol+2]=KS(values,normals)$p.value

  # added to record log2 median ratios and num markers, as downstream tools like gistic2 expect it
  med_ratio=unique(mcols(subsetByOverlaps(ratio.bed,cnvs.bed[i]))$median_ratio)
  num_markers=length(mcols(subsetByOverlaps(ratio.bed,cnvs.bed[i]))$median_ratio)
  if (length(med_ratio) > 1){
    print(paste0("more than one median ratio for segment: ", i, " only saving the first one"))
    print(cnvs.bed[i])
    med_ratio=med_ratio[1]
  } else if(length(med_ratio) == 0){
    print(paste0("no median ratio for segment: ", i))
    print(cnvs.bed[i])
    med_ratio=NA
    num_markers=0
  }
  cnvs[i,numberOfCol+3]=num_markers
  cnvs[i,numberOfCol+4]=log2(med_ratio)
  }

if (numberOfCol==5) {
  names(cnvs)=c("chr","start","end","copy number","status","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue", "num_markers", "log2_median_ratio")
}
if (numberOfCol==7) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue", "num_markers", "log2_median_ratio")
}
if (numberOfCol==9) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","somatic/germline","precentageOfGermline","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue", "num_markers", "log2_median_ratio")
}
write.table(cnvs, file=paste(args[4],".p.value.txt",sep=""),sep="\t",quote=F,row.names=F)


