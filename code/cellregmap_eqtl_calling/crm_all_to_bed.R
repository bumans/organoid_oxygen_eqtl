library(tidyverse)
library(vroom)
library(mashr)

qtl_loc <- snakemake@input[['crm_hits']]
bim_file <- snakemake@input[['bim_file']]
gtf_loc <- snakemake@input[['gtf_loc']]

bed_file <- snakemake@output[['bedfile']]

## Get a mapping from rsID to chromosome positions
### Note that some variants will be missing from this file bc they only pass MAF threshold 
### after subsetting donors in one of the celltypes that didn't contain everyone
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed

# Get a mapping from HGNC to ENSG
gencode <- vroom(gtf_loc) %>% 
  select(c(hgnc, ensg))

# Merge tests from all three trajectories
qtls <- vroom(qtl_loc)
qtls_bed <- qtls %>%
  filter(VARIANT_ID != ".") %>%
  select(c(VARIANT_ID, GENE_HGNC)) %>%
  distinct() %>%
  inner_join(bim, by=c("VARIANT_ID"="variant_id")) %>%
  inner_join(gencode, by=c("GENE_HGNC"="hgnc")) %>%
  dplyr::rename(GENE_ENSG=ensg) %>%
  relocate(`#CHR`, START, END, GENE_ENSG, GENE_HGNC, VARIANT_ID) %>%
  mutate(EB_CONTEXT="CRM") %>%
  write_tsv(bed_file)
