#!/bin/bash

parentDirect=$1
phenoBasename=$2
snpPrefix=$3
genelocs=`basename $4 .locs`
GRM=$5
cisDist=$6
permutations=$7
celltype=$8
#for i in ${phenoDirect}/*.bed
#for i in ${phenoDirect}/expressiontable_matrixeqtl_geschwind_coarse*
#name=expressiontable_matrixeqtl_geschwind_coarse_control10_oRG_cpm-inorm.bed

snpFile=${parentDirect}/snps/${snpPrefix}/
snplocFile=${parentDirect}/snps/${snpPrefix}/
geneLocs=${parentDirect}/snps/${genelocs}
phenoFile=${parentDirect}/expression/${phenoBasename}_${celltype}
covFile=${parentDirect}/covariates/${phenoBasename}_${celltype}
cisOut=${parentDirect}/output/${phenoBasename}_${celltype}
logFile=${parentDirect}/log/${phenoBasename}_${celltype}


#echo $base $phenoFile $logFile
#echo $snpFile $snplocFile $geneLocs $phenoFile $covFile $cisOut $celltype
sbatch --partition=gilad --time=01:00:00 --mem=5G --wrap "Rscript Map_MatrixEQTL_cis_simple.R $snpFile $geneLocs $phenoFile $covFile $GRM $cisOut $cisDist $permutations &> $logFile" 

#echo $snpFile $snplocFile ${geneLocs} $phenoFile $covFile $GRM ${parentDirect}/output/$cisOut ${parentDirect}/output/$cisQQ ${parentDirect}/output/$permutedOut ${parentDirect}/output/$permutedQQ $cisDist $permutations
