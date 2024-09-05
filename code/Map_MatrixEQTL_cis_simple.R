# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
library(tidyverse)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
snpPrefix <- args[1]
# snpLocPrefix <- args[2]
geneLocs <- args[2]
expressionBase <- args[3]
covBase <- args[4]
errorCovariance_file <- args[5]
outputBase <- args[6]
cisDistance <- args[7]
Npermutations <- as.numeric(args[8])
#NB: need to change code for using coarse vs fine annotation


for (chrom in 1:22){
  SNP_file_name <- paste0(snpPrefix, "chr", chrom, "_for_matrixeqtl.snps")
  snps_location_file_name <- paste0(snpPrefix, "chr", chrom, "_for_matrixeqtl.snploc")
  gene_location_file_name <- paste0(geneLocs, "_chr", chrom,".locs")
  expression_file_name <- paste0(expressionBase, "_chr", chrom,"_cpm-inorm.bed")
  # for chromosome-level covariates, use the following
  covariates_file_name <- paste0(covBase, "_chr", chrom,"_cpm-inorm.bed.covs")
  # # for sample-level covariates, use thie following:
  # covariates_file_name <- paste0(covBase, "_cpm-inorm.bed.ruvs.ruvs") 
  output_file_name_cis <- paste0(outputBase, "_chr", chrom,"_cpm-inorm.cis")
  
  useModel = modelLINEAR;
  output_file_name_tra = tempfile();
  pvOutputThreshold_cis = 1;
  errorCovariance <- as.matrix(read.table(errorCovariance_file,sep='\t'))
  
  cisDist = as.numeric(cisDistance);
  print(cisDist)
  print(paste0("chromosome ", chrom))
  ## Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);

  ## Load gene expression data
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
 
  # subset the error covariance file (weights matrix) to include just the
  # individuals present in the expression file (and in the same order)
   fam <- read_table("/project2/gilad/umans/oxygen_eqtl/data/relatedness/YRI_plink.fam", col_names = FALSE)$X1
   # query <- read.table(file = covariates_file_name)[1,-1] %>% unlist() %>% unname()
   query <- colnames(gene)
   errorCovariance <- errorCovariance[match(query, fam), match(query, fam)]
  #errorCovariance <- numeric()
  
  query <- gene$columnNames
  # subset and reorder the snp table to match the gene expression data
  snpnames <- snps$columnNames
  snps$ColumnSubsample(match(query, snpnames));
  
  # reorder the covariate to match the gene expression data
  cvrtnames <- cvrt$columnNames
  cvrt$ColumnSubsample(match(query, cvrtnames));
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  #Include all pvalues in output so that qvalues can be calculated. Filter after.
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = NULL,
    pvOutputThreshold     = 0,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  print('done with real pass')
  
  
  me$cis$eqtls %>% write.table(file = output_file_name_cis, sep='\t', quote = F, row.names = F)
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

  ### Run permutations
  print("Running permutations...")
  for (i in 1:Npermutations){
    # Permute column labels for both (using same seed for randomization)
    # technically I am permuting the column data, and preserving the column labels. That way, I do not have to do the same for the (huge) genotype file, which would take more computational time
    set.seed(i)

    permCols <- sample(query, length(query))
    snpsPerm = snps$Clone()
    snpsPerm$ColumnSubsample(match(permCols, query))
    snpsPerm$columnNames <- query
    
    
    
    
    #Calculate Pvalues from permutated data
    permuted = Matrix_eQTL_main(
      snps = snpsPerm,
      gene = gene,
      cvrt = cvrt,
      output_file_name     = NULL,
      pvOutputThreshold     = 0,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = F,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = 1,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);
    
    print(paste('done with permutation pass', i))
    write.table(permuted$cis$eqtls, paste0(output_file_name_cis, "_permutation_", i, ".out"), row.names = FALSE, quote = FALSE)
  }
}









