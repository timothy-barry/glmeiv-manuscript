#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 48:00:00
#SBATCH --ntasks-per-node=2

source set_xie_args.sh
export NXF_OPTS="-Xms500M -Xmx2G" # limit NF to 2 GB of memory

#######################
# Run Nextflow pipeline
#######################
rm -f trace.txt

# echo args
clear
echo covariate_matrix: $covariate_matrix
echo result_dir: $result_dir
echo gene_odm: $gene_odm
echo gene_metadata: $gene_metadata
echo m lib size: $m_offsets
echo m family: $m_fam
echo gRNA odm: $gRNA_odm
echo gRNA metadata: $gRNA_metadata
echo g lib size: $g_offsets
echo g fam: $g_fam
echo work dir: $work_dir
echo profile: $prof
echo pairs: $gRNA_gene_pairs
echo gene pod size: $gene_pod_size
echo gRNA pod size: $gRNA_pod_size
echo pair pod size: $pair_pod_size

# run the pipeline
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
 --result_file_name $result_file_name \
 -with-trace -w $work_dir $prof
