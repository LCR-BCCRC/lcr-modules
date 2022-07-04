

library(optparse)
library("MutationTimeR")
library(stringr)

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, help="VCF file with mutation data", metavar="character"),
  make_option(c("-c", "--copynumber"), type="character", default=NULL, help="Battenberg subclones output file with copy number data", metavar="character"),
  make_option(c("-l", "--cluster"), type="character", default=NULL, help="Battenberg subclones output file with copy number data", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplename = opt$samplename
vcf_file = opt$vcf
subclones_file = opt$copynumber
cluster_file = opt$cluster
output_dir = opt$output

.checkfile = function(infile) {
  if (!file.exists(infile)) {
    stop(paste("File", infile, "does not exist", sep=""))
  }
}

.checkfile(vcf_file)
.checkfile(cluster_file)
.checkfile(subclones_file)

##################################################
# Process input files
##################################################
# Parse the input file and obtain the required data for this run
vcf_mtr<-readVcf(file =vcf_file)
cnv_mtr<-read.delim(file = subclones_file,header = TRUE)
clust_mtr<-read.delim(file = cluster_file,header = TRUE)

##################################################
# define output files
##################################################
time_segment= file.path(output_dir, paste(samplename,"_copy_number_segments_molecular_time.table", sep = ""))
count_variant= file.path(output_dir, paste(samplename,"_variant_type_counts.table", sep = ""))
type_variant= file.path(output_dir, paste(samplename,"_variants_annotatted.table", sep = ""))
time_segment_plot= file.path(output_dir, paste(samplename,"_copy_number_segments_molecular_time.png", sep = ""))

bb<-makeGRangesFromDataFrame(cnv_mtr,keep.extra.columns=TRUE)

### run mtr main command ##
mt <- mutationTime(vcf_mtr, bb, clusters=clust_mtr, n.boot=10)


### get table of variant_type counts
muttaion_type<-data.frame(table(mt$V$CLS))
rr<-data.frame(apply(muttaion_type,2,function(x)gsub('\\s+', '',x)))
rr$Var1<-gsub("\\[|]","",rr$Var1)
colnames(rr)<-c("variant_type","frequency")
write.table(rr,file = count_variant, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



## time of copy number segments
time_mtr<-data.frame((mt$T))
#cbind(mcols(bb_mm),time_mtr)
time_tab<-cbind((data.frame(bb)),time_mtr)
write.table(time_tab,file = time_segment, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


## variants and their annotations ##
vcf_annotatted <- addMutTime(vcf_mtr, mt$V)
vcf_annotatted_tab<-data.frame(info(vcf_annotatted))
vcf_annotatted_tab_chr<-data.frame("chrom"=(str_split_fixed(rownames(vcf_annotatted_tab), ":", 2)[,1]))
vcf_annotatted_tab_pos<-data.frame("pos"=(str_split_fixed(str_split_fixed(rownames(vcf_annotatted_tab), ":", 2)[,2],"_",2)[,1]))
vcf_annotatted_tabmerged<-cbind(vcf_annotatted_tab_chr,vcf_annotatted_tab_pos,vcf_annotatted_tab)
vcf_annotatted_tabmerged<-data.frame(apply(vcf_annotatted_tabmerged,2,function(x)gsub('\\s+', '',x)))
vcf_annotatted_tabmerged$CLS<-gsub("\\[|]","",vcf_annotatted_tabmerged$CLS)
#rownames(vcf_annotatted_tabmerged) <-NULL

write.table(vcf_annotatted_tabmerged, file = type_variant,col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

vcf <- addMutTime(vcf_mtr, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

#data.frame(mt$T)
#mcols(bb_mm) <- data.frame(cbind(mcols(bb_mm),mt$T))

png(file = time_segment_plot)
plotSample(vcf,bb)
dev.off()
