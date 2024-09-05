import pandas as pd

def list_cellregmap_output_files(wildcards):
    test_genes = set(pd.read_csv(f"/project2/gilad/umans/oxygen_eqtl/topicqtl/mash_and_equivalent_fine_reharmonized.bed", sep="\t")['GENE_HGNC'])
    return [f"/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/fasttopics_fine_15_topics.{g}.cellregmap.tsv" for g in test_genes]

def list_cellregmap_beta_files(wildcards):
    iqtl_genes = set(pd.read_csv(f"/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.sighits.tsv", sep="\t")['GENE_HGNC'])
    return [f"/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/fasttopics_fine_15_topics.{g}.cellregmap.betas.tsv" for g in iqtl_genes]

rule make_genotype_bed:
    resources:
        mem_mb=100000
    input:
        genotypes="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full_combined.vcf.gz",
        inds="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/sample_list"
    output:
        expand("yri_maf0.1_all.hg38.hwe.{out}", out=['bed', 'bim', 'fam'])
    params:
        prefix="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full_combined.vcf.gz"
    shell:
        "code/cellregmap_eqtl_calling/make_genotype_bed.sh {input.genotypes} {input.inds} {params.prefix}"

rule make_kinship:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        genotypes="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full_combined.vcf.gz",
        inds="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/sample_list"
    output:
        expand("/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/yri_kinship.{ext}", ext=['rel', 'rel.id', 'log', 'nosex'])#,
        pruned_bed=expand("/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/pruned_bed.{ext}", ext=['bed', 'log'])
    params:
        kinship_prefix="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/yri_kinship",
        bed_prefix="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/pruned_bed"
    shell:
        "code/cellregmap_eqtl_calling/make_kinship.sh {input.genotypes} {input.inds} {params.kinship_prefix} {params.bed_prefix}"

rule wrangle_kinship:
    resources:
        mem_mb=10000,
        time="1:00:00"
    input:
        kinship="data/genotypes/yri_kinship.rel",
        kinship_id="data/genotypes/yri_kinship.rel.id"
    output:
        kinship_tsv="data/genotypes/yri_kinship.tsv"
    conda: 
        "../environments/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/wrangle_kinship.py"

rule compress_expression:
    resources:
        mem_mb=30000,
        time="1:00:00"
    input:
        pseudocell_data="/project2/gilad/umans/oxygen_eqtl/output/topics_pseudocell_counts_nocontrol21_normalized.tsv"
    output:
        exp="/project2/gilad/umans/oxygen_eqtl/topicqtl/topics_pseudocell_counts_nocontrol21_normalized.nc",
    conda: 
        "../environments/cellregmap.yml"
    script:
        "../code/cellregmap_eqtl_calling/compress_expression.py"


rule run_interaction_test_fasttopics:
    resources:
        partition = "gilad",
        mem_mb = 20000,
        time = "23:00:00"
    input:
        test_eqtl_file="/project2/gilad/umans/oxygen_eqtl/topicqtl/mash_and_equivalent_fine_reharmonized.bed",
        sample_mapping_file = "/project2/gilad/umans/oxygen_eqtl/topicqtl/pseudocell_metadata_r20_harmonized.tsv",
        genotype_file="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full/yri_maf0.1_all.hg38.bed",
        kinship_file = "/project2/gilad/umans/oxygen_eqtl/topicqtl/kinship_tsv",
        exp = "/project2/gilad/umans/oxygen_eqtl/topicqtl/topics_pseudocell_counts_nocontrol21_normalized.nc", 
        cell_contexts = "/project2/gilad/umans/oxygen_eqtl/topicqtl/pseudocell_loadings_k15.tsv"
    output:
        out="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/fasttopics_fine_15_topics.{g}.cellregmap.tsv"
    conda:
        "../environments/cellregmap.yml"
    script:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/code/cellregmap_eqtl_calling/single_gene_interaction_test.py"

rule merge_interaction_test_fasttopics:
    resources:
        partition = "gilad",
        mem_mb = 10000,
        time = "10:00:00"
    input:
        unpack(list_cellregmap_output_files)
    output:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.tsv"
    params:
        tempfile="cellregmap_tempfile.txt"
    shell:
        """
        # Combine input files into a single file
        cat {input} > {params.tempfile}
        
        # Filter lines to unique rows and save to output
        echo -e "GENE_HGNC\tVARIANT_ID\tP_CELLREGMAP" > {output}
        sort -u {params.tempfile} >> {output}
        
        # Remove temp file
        rm {params.tempfile}
        """

rule cellregmap_mtc:
    resources:
        mem_mb = 10000,
        time = "10:00"
    input:
        eqtls="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.tsv"
    output:
        all_tests="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.mtc.tsv",
        top_tests="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.tophits.tsv",
        sig_hits="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.sighits.tsv"
    conda:
        "../environments/r-mashr-updated.yml"
    script:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/code/cellregmap_eqtl_calling/cellregmap_mtc.R"

rule estimate_effect_sizes:
    resources:
        mem_mb = 20000,
        time = "24:00:00"
    input:
        interaction_eqtl_file="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.sighits.tsv",
        sample_mapping_file = "/project2/gilad/umans/oxygen_eqtl/topicqtl/pseudocell_metadata_r20_harmonized.tsv",
        genotype_file="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full/yri_maf0.1_all.hg38.bed",
        kinship_file = "/project2/gilad/umans/oxygen_eqtl/topicqtl/kinship_tsv",
        exp = "/project2/gilad/umans/oxygen_eqtl/topicqtl/topics_pseudocell_counts_nocontrol21_normalized.nc",
        cell_contexts = "/project2/gilad/umans/oxygen_eqtl/topicqtl/pseudocell_loadings_k15.tsv" 
    output:
        out="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/fasttopics_fine_15_topics.{g}.cellregmap.betas.tsv"
    conda:
        "../environments/cellregmap.yml"
    script:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/code/cellregmap_eqtl_calling/single_gene_effect_estimates.py"

rule merge_effect_size_estimates:
    resources:
        mem_mb = 10000,
        time = "10:00:00"
    input:
        unpack(list_cellregmap_beta_files)
    output:
        "temp/cellregmap_eqtl_calling.fasttopics_15_topics.cellregmap.DUMMYDONEFILE.tsv"
    shell:
        """
        echo done > {output}
        """

rule crm_to_bed:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        crm_hits="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.sighits.tsv",
        bim_file="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full/yri_maf0.1_all.hg38.hwe.bim",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        bedfile="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/crm-signif_variant_gene_pairs.bed"
    conda: "../environments/r-mashr-updated.yml"
    script:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/code/static_eqtl_followup/crm_to_bed.R"
        

rule cellregmap_mtc_all_signif:
    resources:
        mem_mb = 10000,
        time = "10:00"
    input:
        eqtls="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.tsv"
    output:
        all_qtls="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.all_signif.fdr0.1.tsv"
    conda:
        "../environments/r-mashr-updated.yml"
    script:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/code/cellregmap_eqtl_calling/cellregmap_mtc_all.R"

rule crm_full_to_bed:
    resources:
        mem_mb=50000,
        time="15:00"
    input:
        crm_hits="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/all_genes_merged_fine_fasttopics_15_topics.cellregmap.all_signif.fdr0.1.tsv",
        bim_file="/project2/gilad/umans/oxygen_eqtl/data/MatrixEQTL/snps/YRI_genotypes_maf10hwee-6_full/yri_maf0.1_all.hg38.hwe.bim",
        gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        bedfile="/project2/gilad/umans/oxygen_eqtl/topicqtl/outputs/crm-all-signif_variant_gene_pairs.bed"
    conda: "../environments/r-mashr-updated.yml"
    script:
        "/project2/gilad/umans/oxygen_eqtl/topicqtl/code/cellregmap_eqtl_calling/crm_all_to_bed.R"

