#!/bin/bash

fasta_stats="bash /home/fad19/scratch/scripts/random_scripts/reads_simple_stats.sh"
salmonella_dir="/home/fad19/scratch/protein_aligner/paper_results/salmonella_sequences"
filtered_dir="/home/fad19/scratch/protein_aligner/paper_results/salmonella_sequences/filtered_alignments_50_percent_length"
filtered_50="/home/fad19/scratch/protein_aligner/paper_results/salmonella_sequences/filtered_alignments_50_percent_length/panpa_50_comparisons"
filtered_70="/home/fad19/scratch/protein_aligner/paper_results/salmonella_sequences/filtered_alignments_50_percent_length/panpa_70_comparisons"


####################### raw
echo "Salmonella amino acid sequences stats"
$fasta_stats $salmonella_dir/salmonella_refseq_annotation_aa.fasta

echo ""
echo "Salmonella DNA sequences stats"
$fasta_stats $salmonella_dir/salmonella_refseq_annotation_dna.fasta

echo ""
echo "BWA raw alignments number"
sort $salmonella_dir/salmonella_refseq_annotation_dna_ecoli_ref_genome_bwa.sam | uniq | awk '{if ($2 != 4) {counter += 1}} END {print counter}'

echo ""
echo "GraphAligner raw alignments number"
sort $salmonella_dir/salmonella_refseq_annotation_dna_ecoli_pangenome.gaf | uniq | wc -l 


echo ""
echo "PanPA raw alignments number with index k 5 w 5 index seed limit 0 seed limit 3"
sort $salmonella_dir/salmonella_refseq_annotation_aa_stats_alignment_k_5_w_5_seed_lim_3.gaf | uniq | wc -l 


echo ""
echo "PanPA raw alignments number with index k 5 w 5 index seed limit 0 seed limit 10"
sort $salmonella_dir/salmonella_refseq_annotation_aa_stats_alignment_k_5_w_5_seed_lim_10.gaf | uniq | wc -l 


####################### filtered at 50% alignment length
echo ""
echo "BWA alignments number filtered 50 percent alignment length"
sort $filtered_dir/salmonella_refseq_annotation_dna_ecoli_ref_genome_bwa_filtered_length_50_percent.sam | uniq | awk '{if ($2 != 4) {counter += 1}} END {print counter}'

echo ""
echo "GraphAligner alignments number filtered 50 percent alignment length"
sort $filtered_dir/salmonella_refseq_annotation_dna_ecoli_pangenome_filtered_length_50_percent.gaf | uniq | wc -l 


echo ""
echo "PanPA alignments number with index k 5 w 5 index seed limit 0 seed limit 3 filtered 50 percent alignment length"
sort $filtered_dir/salmonella_refseq_annotation_aa_stats_alignment_k_5_w_5_seed_lim_3_filtered_length_50_percent.gaf | uniq | wc -l 


echo ""
echo "PanPA alignments number with index k 5 w 5 index seed limit 0 seed limit 10 filtered 50 percent alignment length"
sort $filtered_dir/salmonella_refseq_annotation_aa_stats_alignment_k_5_w_5_seed_lim_10_filtered_length_50_percent.gaf | uniq | wc -l 

####################### filtered at 50% alignment length and 50% alignment id

echo ""
echo "BWA alignments number filtered 50 percent alignment length, 50 percent alignment id"
sort $filtered_50/salmonella_refseq_annotation_dna_ecoli_ref_genome_bwa_filtered_id_50_length_50_percent.sam| uniq | awk '{if ($2 != 4) {counter += 1}} END {print counter}'

echo ""
echo "GraphAligner alignments number filtered 50 percent alignment length, 50 percent alignment id"
sort $filtered_50/salmonella_refseq_annotation_dna_ecoli_pangenome_id_50_filtered_length_50_percent.gaf | uniq | wc -l 


echo ""
echo "PanPA alignments number with index k 5 w 5 index seed limit 0 seed limit 10 filtered 50 percent alignment length, 50 percent alignment id"
sort $filtered_50/salmonella_refseq_annotation_aa_stats_alignment_k_5_w_5_seed_lim_10_id_50_filtered_length_50_percent.gaf | uniq | wc -l 


####################### filtered at 70% alignment length and 50% alignment id

echo ""
echo "BWA alignments number filtered 50 percent alignment length, 70 percent alignment id"
sort $filtered_70/salmonella_refseq_annotation_dna_ecoli_ref_genome_bwa_filtered_id_70_length_50_percent.sam | uniq | awk '{if ($2 != 4) {counter += 1}} END {print counter}'

echo ""
echo "GraphAligner alignments number filtered 50 percent alignment length, 70 percent alignment id"
sort $filtered_70/salmonella_refseq_annotation_dna_ecoli_pangenome_id_70_filtered_length_50_percent.gaf | uniq | wc -l 


echo ""
echo "PanPA alignments number with index k 5 w 5 index seed limit 0 seed limit 10 filtered 70 percent alignment length, 50 percent alignment id"
sort $filtered_70/salmonella_refseq_annotation_aa_stats_alignment_k_5_w_5_seed_lim_10_id_70_filtered_length_50_percent.gaf | uniq | wc -l 
