# INSTALLATION GUIDE:

# (1) use conda:
# -download ichorcna env
# -download cmake env
	
# (2) git clone hmmcopy_utils into envs/ichorcna/hmmcopy_utils/
	# cd hmmcopy_utils
	# cmake .
	# make

# (3) git clone git@github.com:broadinstitute/ichorCNA.git into envs/ichorcna/ichorCNA/
	# install R dependencies (R-3.6.0 +):
		## install from CRAN
		# install.packages("plyr") 
		## install packages from
 		# install.packages("BiocManager")
 		# BiocManager::install("HMMcopy")  
 		# BiocManager::install("GenomeInfoDb")  
 		# BiocManager::install("GenomicRanges")  
	# install R package
		# cd envs/ichorcna/
		# bin/R CMD INSTALL ichorCNA  

