####################################
###### FUNCTIONS ###################
####################################


`%not_in%` <- purrr::negate(`%in%`)

add_vireo <- function(obj, dir) {
  donor_ids <- read.table(paste0(dir, 'vireo/donor_ids.tsv'), header = TRUE, stringsAsFactors = FALSE)
  
  obj <- AddMetaData(obj, donor_ids$donor_id,   col.name ='vireo.individual')
  obj <- AddMetaData(obj, donor_ids$prob_max,   col.name ='vireo.prob.singlet')
  
  ambient.file <- paste0(dir, 'vireo/prop_ambient.tsv')
  if (file.exists(ambient.file)) {
    prop <- read.table(paste0(dir, 'vireo/prop_ambient.tsv'), header = TRUE, stringsAsFactors = FALSE)
    
    inds <- unique(donor_ids$best_singlet)
    donor_ids$prop_donor <- 0
    for (i in 1:length(inds)) {
      w <- which(donor_ids$best_singlet==inds[i])
      donor_ids$prop_donor[w] <- prop[w,inds[i]]
    }
    obj <- AddMetaData(obj, donor_ids$prop_donor, col.name ='vireo.prop.donor')
  }
  obj
}


generate.pseudobulk <- function(object, labels, assay="RNA", slot="counts") {
  factorlist <- list()
  for (i in labels) factorlist[[i]] <- unique(object@meta.data[,i])
  meta <- expand.grid(factorlist, stringsAsFactors = FALSE)
  rownames(meta) <- apply(meta, 1, function(x) paste0(x, collapse = '.'))
  
  extrameta <- data.frame(cycles = rep(NA, nrow(meta)),
                          cdna = rep(NA, nrow(meta)),
                          fpassage = rep(NA, nrow(meta)),
                          ffpassage = rep(NA, nrow(meta)),
                          ncells = rep(NA, nrow(meta)),
                          mt = rep(NA, nrow(meta)),
                          umi = rep(NA, nrow(meta)),
                          gene = rep(NA, nrow(meta))
  )
  rownames(extrameta) <- rownames(meta)
  
  # build the output matrix
  n <- nrow(meta)
  out <- matrix(nrow=dim(object[[assay]])[1], ncol=n, data=0)
  rownames(out) <- rownames(object[[assay]])
  colnames(out) <- rownames(meta)
  # ncells <- c()
  # ncounts <- c()
  total.cells <- dim(object[[assay]])[2]
  for (i in 1:n)
  {
    #prog(i,n)
    cells <- 1:total.cells
    for (j in names(meta)) {
      keep  <- which(object@meta.data[[j]] == meta[i,j])
      cells <- cells[cells %in% keep]
    }
    # ncells[i] <- length(cells)
    extrameta[i, "ncells"] <- length(cells)
    
    #some other thing to measure
    if (length(cells)==1) {
      out[,i] <- slot(object[[assay]], slot)[,cells]
      extrameta[i, "cycles"] <- object@meta.data$cycles[cells]
      extrameta[i, "cdna"] <- object@meta.data$cDNA[cells]
      extrameta[i, "mt"] <- object@meta.data$percent.mt[cells]
      extrameta[i, "fpassage"] <- object@meta.data$Fpassage[cells]
      extrameta[i, "ffpassage"] <- object@meta.data$FFpassage[cells]
      extrameta[i, "umi"] <- object@meta.data$nCount_RNA[cells]
      extrameta[i, "gene"] <- object@meta.data$nFeature_RNA[cells]
    } else {
      out[,i] <- Matrix::rowSums(slot(object[[assay]], slot)[,cells])
      # meta[i, "ncounts"] <- sum(slot(object[[assay]], slot)[,cells])/max(sum(slot(object[[assay]], slot)[,cells]))
      extrameta[i, "cycles"] <- mean(object@meta.data$cycles[cells])
      extrameta[i, "cdna"] <- mean(object@meta.data$cDNA[cells])
      extrameta[i, "mt"] <- mean(object@meta.data$percent.mt[cells])
      extrameta[i, "fpassage"] <- mean(object@meta.data$Fpassage[cells])
      extrameta[i, "ffpassage"] <- mean(object@meta.data$FFpassage[cells])
      extrameta[i, "umi"] <- median(object@meta.data$nCount_RNA[cells])
      extrameta[i, "gene"] <- median(object@meta.data$nFeature_RNA[cells])
    }
  }
  # meta$ncells <- ncells
  # meta$ncounts <- ncounts/max(ncounts)
  #add that something else as metadata
  return(list(counts=out, meta=cbind(meta, extrameta)))
}


filter.pseudobulk <- function(pseudobulk, threshold = 0) {
  w <- which(pseudobulk$meta$ncells > threshold)
  pseudobulk$counts <- pseudobulk$counts[,w]
  pseudobulk$meta <- pseudobulk$meta[w,]
  pseudobulk
}




scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

combined_across_chroms <- function(results_directory, results_basename, condition, celltype, output_basename="results_combined_"){
  if(file.exists(paste0(results_directory, results_basename, condition, "_", celltype, "_chr1_cpm-inorm.cis"))){
    m = read.table(paste0(results_directory, results_basename, condition, "_", celltype, "_chr1_cpm-inorm.cis"), head=TRUE, stringsAsFactors=FALSE)
  }else{
    print(paste0(results_directory, results_basename, condition, "_", celltype, "_chr1_cpm-inorm.cis does not exist!"))
  }
  for (i in 2:22){
    if(file.exists(paste0(results_directory, results_basename, condition, "_", celltype, "_chr", i, "_cpm-inorm.cis"))){
      m = rbind(m, read.table(paste0(results_directory, results_basename, condition, "_", celltype, "_chr", i, "_cpm-inorm.cis"), head=TRUE, stringsAsFactors=FALSE))
    }else{
      print(paste0(results_directory, results_basename, condition, "_", celltype, "_chr", i, "_cpm-inorm.cis does not exist!"))
    }
  }
  write.table(m, paste0(results_directory, output_basename, condition, "_", celltype, "_nominal.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, append = FALSE)  
  print(paste0("done with ", celltype, " ", condition))
  rm(m)
}

plot_cisqtl_cvrt <- function(expression.data, genotype.data, snp, gene, covariate_file){
  test.genotypes <- genotype.data[which(genotype.data$ID==snp),]
  test.expression <- expression.data[which(expression.data$ID==gene),]
  cvrt <- read_table(file = covariate_file, col_names = TRUE, show_col_types = FALSE) %>%
    pivot_longer(-id) %>% 
    pivot_wider(names_from=id, values_from=value)
  
  for_regression <- left_join(test.expression, cvrt, by=c("individual"="name"))
  for_regression$resids <- lm(formula = as.formula(paste("expression ~ ", paste(colnames(cvrt)[-1], collapse= "+"))), data = for_regression)$residuals
  combined <- left_join(for_regression, test.genotypes, by="individual")

  ggplot(combined, mapping = aes(x=factor(snp), y=resids)) + geom_point(alpha=0.5, position = position_jitter(width = 0.2, height = 0)) + ggtitle(gene) + xlab(snp) + theme_light()
}

make_boxplot <- function(celltype, condition, testsnp, testgene){
  expression <- read_table(paste0("data/MatrixEQTL/expression/combined_fine_quality_filter20_032024/expressiontable_matrixeqtl_combined_fine_", condition, "_", celltype, "_chr1_cpm-inorm.bed")) %>% pivot_longer(cols=starts_with("NA"), names_to = "individual", values_to = "expression") 
  for (i in 2:22){
    expression <- rbind(expression, read_table(paste0("data/MatrixEQTL/expression/combined_fine_quality_filter20_032024/expressiontable_matrixeqtl_combined_fine_", condition, "_", celltype, "_chr",i, "_cpm-inorm.bed"), show_col_types = FALSE) %>% pivot_longer(cols=starts_with("NA"), names_to = "individual", values_to = "expression"))
  }
  
  covariate_file <- paste0("/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/covariates/combined_fine_quality_filter20_032024/expressiontable_matrixeqtl_combined_fine_", condition, "_", celltype, "_chr1_cpm-inorm.bed.covs")
  # uniquely weighted
  plot_cisqtl_cvrt(expression.data = expression, genotype.data = genotypes, snp = testsnp, gene=testgene, covariate_file = covariate_file) + ggtitle(paste0(testgene, " ", celltype, " ", condition))
}



####################################
###### STYLE CONVENTIONS ###########
####################################


manual_palette_fine <- c("#56beac", "#8DD3C7", "#c4e8e2",
                         "#8f88bf", "#BEBADA",
                         "#FFFFB3", 
                         "#FB8072", 
                         "#fc9016", "#fed8ae", # "#FDB462"
                         "#478EBF", "#80B1D3", "#a6c8e0",
                         "#AFAFAF", "#323232", "#D9D9D9",  "#858585" ,"#5B5B5B",  
                         "#B3DE69", 
                         "#FCCDE5", "#BC80BD")
names(manual_palette_fine) <- c("RGcycling" , "RG", "CorticalHem", 
                                "IPcycling", "IP", "GliaProg", 
                                "Immature", "Glut", "GlutNTS", 
                                "NeuronOther", "Cajal", "MidbrainDA", 
                                "Inh",  "InhThalamic", "InhGNRH", 
                                "InhSST", "InhMidbrain", "VLMC",
                                "Oligo", "Choroid")

fine.order <- c("RGcycling" , "RG", "CorticalHem",
                "IPcycling", "IP",
                "GliaProg", "Oligo",
                "Immature", "Glut", "GlutNTS", 
                "NeuronOther", "MidbrainDA",  "Cajal",
                "Inh",  "InhThalamic",
                "InhGNRH",  "InhSST", "InhMidbrain", 
                "Choroid",
                "VLMC")


manual_palette_coarse <- c( "#8DD3C7", 
                            "#BEBADA",
                            "#FFFFB3", 
                            "#FB8072", 
                            "#fc9016",
                            "#478EBF", "#80B1D3",
                            "#AFAFAF",
                            "#B3DE69", 
                            "#BC80BD")
names(manual_palette_coarse) <- c( "RG", "IP", "Glia", 
                                   "Immature", "Glut",  "NeuronOther", 
                                   "Cajal",  "Inh",   "VLMC", "Choroid")

coarse.order <- c("RG", "IP", "Glia", "Immature", "Glut", "NeuronOther", "Cajal","Inh",  "VLMC", "Choroid")


class_colors <- c("class1"= "lightblue",
                  "class2"="pink",
                  "class3"="red", 
                  "class4"="darkblue")


