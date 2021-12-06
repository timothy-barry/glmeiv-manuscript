#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=2

source set_gasp_args.sh
export NXF_OPTS="-Xms500M -Xmx2G" # limit NF to 2 GB of memory

#######################
# Run Nextflow pipeline
#######################
rm -f trace.txt

# nextflow run $glmeiv_nf_pipeline \
nextflow run https://github.com/timothy-barry/glmeiv-pipeline -r main \
--pairs $gRNA_gene_pairs \
--covariate_matrix $covariate_matrix \
--pair_pod_size $pair_pod_size \
--result_dir $result_dir \
--gene_pod_size $gene_pod_size \
--gene_odm $gene_odm \
--gene_metadata $gene_metadata \
--m_offsets $m_offsets \
--m_fam_str $m_fam \
--gRNA_pod_size $gRNA_pod_size \
--gRNA_odm $gRNA_odm \
--gRNA_metadata $gRNA_metadata \
--g_offsets $g_offsets \
--g_fam_str $g_fam \
-with-trace -w $work_dir
