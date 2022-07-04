
library(optparse)

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
  make_option(c("-m", "--maf"), type="character", default=NULL, help="VCF file with mutation data", metavar="character"),
  make_option(c("-p", "--cellularity"), type="character", default=NULL, help="Battenberg rho and psi output file", metavar="character"),
  make_option(c("-c", "--copynumber"), type="character", default=NULL, help="Battenberg subclones output file with copy number data", metavar="character"),
  make_option(c("-l", "--cluster"), type="character", default=NULL, help="Battenberg subclones output file with copy number data", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplename = opt$samplename
maf_file = opt$maf
cellularity_file = opt$cellularity
subclones_file = opt$copynumber
cluster_file = opt$cluster
output_dir = opt$output

.checkfile = function(infile) {
  if (!file.exists(infile)) {
    stop(paste("File", infile, "does not exist", sep=""))
  }
}

.checkfile(maf_file)
.checkfile(cellularity_file)
.checkfile(subclones_file)

##################################################
# Process input files
##################################################
# Parse the input file and obtain the required data for this run
maf<-read.delim(file =maf_file, header =TRUE,sep = "\t")
cellularity <- read.delim(file = cellularity_file,header = TRUE)
cnv<-read.delim(file = subclones_file,header = TRUE)
clust<-read.delim(file = cluster_file,header = TRUE)

##################################################
# define output files
##################################################
out_vcf= file.path(output_dir, paste(samplename,"_maf_convertedTo.vcf", sep = ""))
out_clust= file.path(output_dir, paste(samplename,"_cluster_mtr.txt", sep = ""))
out_cnv= file.path(output_dir, paste(samplename,"_cnv_mtr.txt", sep = ""))

##############################################
# convert maf to required vcf format
##############################################

mafTovcf = function(maf_file= maf_file,outfile = out_vcf){
  maf1<-maf[maf$VARIANT_CLASS == "SNV",]
  maf1$Chromosome<-gsub("chr","",maf1$Chromosome)
  if (nrow(maf1)==0) {
    print ("no variant in the maf file")
   } else {

  maf2<-maf1[,c("Chromosome","vcf_pos", "Reference_Allele","Tumor_Seq_Allele2","t_ref_count", "t_alt_count","t_depth","Variant_Type")]
  maf2$t_alt_count_m<-paste("t_alt_count=",maf2$t_alt_count, sep = "")
  maf2$t_ref_count_m<-paste("t_ref_count=",maf2$t_ref_count, sep = "")
  maf2$FORMAT<-rep("AD:DP", nrow(maf2))
  maf2$SAMPLE<-paste(paste(".",maf2$t_alt_count,sep = ","),maf2$t_depth,sep = ":")
  maf2$INFO<-paste(maf2$t_alt_count_m,maf2$t_ref_count_m, sep =";")

  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
  #1	47921661	.	N	A	.	.	t_alt_count=33;t_ref_count=59	AD:DP	.,33:92

  maf2$ID<-rep(".",nrow(maf2))
  maf2$QUAL<-rep(".",nrow(maf2))
  maf2$FILTER<-rep(".",nrow(maf2))

  maf4<-maf2[,c("Chromosome","vcf_pos","ID","Reference_Allele","Tumor_Seq_Allele2","QUAL","FILTER", "INFO","FORMAT","SAMPLE")]
  colnames(maf4)<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

  out_vcf=out_vcf
   cat("##fileformat=VCFv4.1",
      "\n##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allelic depths (number of reads in each observed allele)\">",
      "\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">",
      "\n##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Variant filters\">",
      "\n##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description=\"Tumour ref count\">",
      "\n##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description=\"Tumour alt count\">",
      "\n",
      #file = shared_samples[j])
      file = out_vcf)

      write.table(maf4, file = outfile,append=TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

   } ## if statement 

}

mafTovcf(maf_file=maf_file, outfile=out_vcf)

##############################################
# convert maf to required vcf format
##############################################
#purity<-read.table(file = "PD26403d--PD26403b_cellularity_ploidy.txt", header = TRUE)
#cnv<-read.delim(file = "PD26403d--PD26403b_subclones.txt", header = TRUE, sep = "\t")

makeCNV =function(subclones_file=subclones_file , cellularity_file = cellularity_file, outfile = out_cnv){
 
 cnv1<-cnv[,c("chr","startpos","endpos","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A")]
 cnv1$chr<-gsub("chr","",cnv1$chr)
 #cnv1$width<-cnv1$endpos-cnv1$startpos

 out_res<-NULL
 for (i in 1:nrow(cnv)){

   cnv_foc<-cnv1[i,] 
   if (is.na(cnv_foc$nMaj2_A)==TRUE){
   cnv_foc$width<-cnv_foc$endpos-cnv_foc$startpos
   cnv_foc$strand<-rep("*",nrow(cnv_foc))
   cnv_foc$clonal_frequency = rep(cellularity$cellularity)
   cnv3<-cnv_foc[,c("chr","startpos","endpos","width","strand", "nMaj1_A","nMin1_A","clonal_frequency")]
   colnames(cnv3)<-c("seqnames","start","end","width","strand", "major_cn","minor_cn", "clonal_frequency")
   } ## forst if 

   if (is.na(cnv_foc$nMaj2_A)==FALSE){

   cnv_foc$width<-cnv_foc$endpos-cnv_foc$startpos
   cnv_foc$strand<-rep("*",nrow(cnv_foc))
   cnv_foc$width<-cnv_foc$endpos-cnv_foc$startpos
   
   cnv_foc_sub<-cnv_foc[rep(1, 1), ]
   
   cnv_foc$clonal_frequency = (cellularity$cellularity*cnv_foc$frac1_A)
   cnv_foc_clon<-cnv_foc[,c("chr","startpos","endpos","width","strand", "nMaj1_A","nMin1_A","clonal_frequency")]
   colnames(cnv_foc_clon)<-c("seqnames","start","end","width","strand", "major_cn","minor_cn", "clonal_frequency")

   cnv_foc_sub$clonal_frequency = (cellularity$cellularity - cnv_foc_clon$clonal_frequency)
   cnv3_sub<-cnv_foc_sub[,c("chr","startpos","endpos","width","strand", "nMaj2_A","nMin2_A","clonal_frequency")]
   colnames(cnv3_sub)<-c("seqnames","start","end","width","strand", "major_cn","minor_cn", "clonal_frequency")

   cnv3<-rbind(cnv_foc_clon,cnv3_sub)
   colnames(cnv3)<-c("seqnames","start","end","width","strand", "major_cn","minor_cn", "clonal_frequency")
  
 }

  out_res<-rbind(out_res,cnv3)
  write.table(out_res, file =outfile , row.names = FALSE, sep = "\t", quote = FALSE)
 }

}  ## fucnction
makeCNV(subclones_file=subclones_file, cellularity_file = cellularity_file, outfile=out_cnv)

##############################################
# correct cluster number
##############################################

#clust<-read.table(file = "PD26403d_2000iters_1000burnin_bestClusterInfo.txt", header = TRUE)
# pattern
##   cluster n_ssms proportion
## 1       0   5518       0.54
## 2       1    845       0.23

makeCluster<-function(cluster_file,cellularity_file,outfile){
 clust$cluster<-seq(from = 0, to = nrow(clust)-1)
 clust$proportion = ifelse(clust$cluster==0,cellularity$cellularity,cellularity$cellularity*clust$location)
 #clust$proportion<-purity$cellularity
 clust2<-clust[,c("cluster","no.of.mutations","proportion")]
 colnames(clust2)<-c("cluster","n_ssms","proportion")
 write.table(clust2, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
}

makeCluster(cluster_file=cluster_file, cellularity_file = cellularity_file, outfile=out_clust)





