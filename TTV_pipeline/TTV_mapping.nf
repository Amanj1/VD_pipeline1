#!/usr/bin/env nextflow

/*
How to run:
nextflow -C discovery.nf.config run discovery.nf --fastq_files preprocessing -profile hamlet
*/

params.fastq_dir='preprocessing'
params.contigs_dir='preprocessing'

reads_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,unpaired}.fq.gz",size:3)
contigs_files = Channel.fromFilePairs("${params.contigs_dir}/**/*_filt_contigs.fa",size:1)

reads_files.into{
 r_counter_in;
 r_BWA_in;
 r_blast_in
}
 
contigs_files.into{
 asm_counter_in;
 asm_blastn_in
}

/**************************************************************************************************
 CONTIGS
**************************************************************************************************/

process contig_counter{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/contigs", mode:'link'

  input:
  set sample_id, contig_file from asm_counter_in

  output:
  set sample_id, "${sample_id}_nr_of_seq_contigs.txt" into contig_nr_seq

  script:
  """
  cat ${contig_file[0]} | grep ">" | wc -l > ${sample_id}_nr_of_seq_contigs.txt
  """
}

/**************************************************************************************************
 CONTIGS MAGICBLAST
**************************************************************************************************/

process contig_magicblast{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/contigs/magicblast", mode:'link'

  input:
  set sample_id, contig_file from asm_blastn_in

  output:
  set sample_id, "${sample_id}_mapped_results_contigs.magicblast" into contig_blast_out

  script:
  """
  magicblast -db picornavirus_lineage -query ${contig_file[0]} -outfmt tabular -no_unaligned -num_threads ${task.cpus} -splice false -score 125 > ${sample_id}_mapped_results_contigs.magicblast
  """
}

contig_blast_out.into{
 contig_blast_counter;
 contig_blast_selection
}

process contig_magicblast_counter{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/contigs/magicblast", mode:'link'

  input:
  set sample_id, mapped_seq from contig_blast_counter

  output:
  set sample_id, "${sample_id}_nr_of_seq_contigs_magicblast.txt" into contig_blast_nr_seq

  script:
  """
  cat ${mapped_seq} | awk '{print \$1}' | awk '!seen[\$x]++' | wc -l > ${sample_id}_nr_of_seq_contigs_magicblast.txt
  """
}


process contig_magicblast_selection{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/contigs/magicblast", mode:'link'

  input:
  set sample_id, blast from contig_blast_selection

  output:
  set sample_id, "${sample_id}_occurrences_of_matched_ref_magicblast.txt", "${sample_id}_selected_contig_magicblast.txt" into contig_magicblast_selection_out

  script:
  """
  awk '{print \$2 "\t" \$1 "\t" \$13}' ${blast} | tail -n +2 | tail -n +2 | tail -n +2| datamash -sg 2 max 3 -f | awk '{print \$1 "\t" \$2 "\t" \$3}' > ${sample_id}_selected_contig_magicblast.txt
  awk '{print \$2 "\t" \$1 "\t" \$13}' ${blast} | tail -n +2 | tail -n +2 | tail -n +2| datamash -sg 2 max 3 -f | awk '{print \$1 "\t" \$2 "\t" \$3}' | datamash -g 1 count 2 > ${sample_id}_occurrences_of_matched_ref_magicblast.txt
  """
}

contig_combine_2_magicblast_channels = contig_magicblast_selection_out.combine(contig_blast_nr_seq, by: 0)
contig_combine_3_magicblast_channels = contig_combine_2_magicblast_channels.combine(contig_nr_seq, by: 0)


process contig_magicblast_final_sampleResults{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/contigs/magicblast", mode:'link'

  input:
  set sample_id, occur_of_matched_ref, selected_contig, nr_of_tot_contig_mapped, nr_of_tot_contig_fasta from contig_combine_3_magicblast_channels

  output:
  set sample_id, "${sample_id}_sample_results_magicblast_contig.txt" into contig_magicblast_sample_results

  script:
  """
  #!/bin/bash
  
  touch ${sample_id}_sample_results_magicblast_contig.txt
  
  nr_of_tot_fq=\$(head -n 1 ${nr_of_tot_contig_fasta})
  nr_of_tot_map=\$(head -n 1 ${nr_of_tot_contig_mapped})

  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_occur=\$(echo \$line |awk '{print \$2}')
		ref_name=\$(blastdbcmd -entry \$ref_id -db picornavirus_lineage -range 1-1 | tr -d '\n' | sed 's/.\$//')
		ref_name=\${ref_name// /_}
		ref_name=\${ref_name//>/}
		ref_name=\${ref_name//:1-1_/}
		ref_name=\${ref_name//\$ref_id/}
		ref_occur_perTTV=\$(echo -e \$ref_occur/\$nr_of_tot_map "\n" scale=2 | bc -l)
		ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_of_tot_fq "\n" scale=2 | bc -l)
		TTV_perALL=\$(echo -e \$nr_of_tot_map/\$nr_of_tot_fq "\n" scale=2 | bc -l)
		echo -e "${sample_id}\t\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$nr_of_tot_map\t\$nr_of_tot_fq" >> ${sample_id}_sample_results_magicblast_contig.txt
  done < ${occur_of_matched_ref}
  """
}

contig_magicblast_sample_results.into{
 contig_magicblast_sample_results1;
 contig_magicblast_sample_results2
}

process contig_magicblast_final_Results{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/contigs/magicblast", mode:'link'

  input:
  set sample_id, "all_res" from contig_magicblast_sample_results1

  output:
  file "${sample_id}_final_results_magicblast_contig.tsv"  into contig_magicblast_final_sample

  script:
  """
  #!/bin/bash
  cat ${all_res} > tmp.txt
  awk '{print \$1 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1 sum 3 -f | datamash sum 4,5 > max_values_TTV_and_allSeq.txt
  TTV_occur=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$1}')
  nr_reads_fastq=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$2}')
  rm max_values_TTV_and_allSeq.txt
  echo -e "Ref_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_contigs" > ${sample_id}_final_results_magicblast_contig.tsv
  cat tmp.txt | awk '{print \$2 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' | datamash -sg 1,2 sum 3,4,5 > ${sample_id}_results_summed_colmn.2.3.7.8.9.magicblast.txt
  //rm tmp.txt
  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_name=\$(echo \$line |awk '{print \$2}')
        ref_occur=\$(echo \$line |awk '{print \$3}')
        ref_occur_perTTV=\$(echo -e \$ref_occur/\$TTV_occur "\n" scale=2 | bc -l)
        ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        TTV_perALL=\$(echo -e \$TTV_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        echo -e "\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$TTV_occur\t\$nr_reads_fastq" >> ${sample_id}_final_results_magicblast_contig.tsv
  done < ${sample_id}_results_summed_colmn.2.3.7.8.9.magicblast.txt
  """
}

process contig_magicblast_concat_allResults{
  tag {"All_samples"}

  publishDir "${params.publish_base_dir}/all/contigs/magicblast", mode:'link'

  input:
  file "sample_res" from contig_magicblast_sample_results2.map{it[1]}.collect()

  output:
  file "All_results_combined_magicblast_contig.tsv" into contig_magicblast_sample_results_concat

  script:
  """
  #!/bin/bash
  touch All_results_combined_magicblast.tsv
  echo -e "Sample_id\tRef_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_contigs" >> All_results_combined_magicblast_contig.tsv
  for sample_file in ${sample_res}
  do
        cat \$sample_file >> All_results_combined_magicblast_contig.tsv
  done
  """
}

process contig_magicblast_final_allResults{
  tag {"All_samples"}

  publishDir "${params.publish_base_dir}/all/contigs/magicblast", mode:'link'

  input:
  file "all_res" from contig_magicblast_sample_results_concat

  output:
  file "All_final_results_magicblast_contig.tsv"  into contig_magicblast_final_all

  script:
  """
  #!/bin/bash
  tail -n +2 ${all_res} > tmp.txt
  awk '{print \$1 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1 sum 3 -f | datamash sum 4,5 > max_values_TTV_and_allSeq.txt
  TTV_occur=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$1}')
  nr_reads_fastq=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$2}')
  rm max_values_TTV_and_allSeq.txt
  echo -e "Ref_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_contigs" > All_final_results_magicblast_contig.tsv
  awk '{print \$2 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1,2 sum 3,4,5 > All_results_summed_colmn.2.3.7.8.9.magicblast.txt
  //rm tmp.txt
  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_name=\$(echo \$line |awk '{print \$2}')
        ref_occur=\$(echo \$line |awk '{print \$3}')
        ref_occur_perTTV=\$(echo -e \$ref_occur/\$TTV_occur "\n" scale=2 | bc -l)
        ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        TTV_perALL=\$(echo -e \$TTV_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        echo -e "\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$TTV_occur\t\$nr_reads_fastq" >> All_final_results_magicblast_contig.tsv
  done < All_results_summed_colmn.2.3.7.8.9.magicblast.txt
  """
}

/**************************************************************************************************
 READS
**************************************************************************************************/

process reads_counter{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from r_counter_in

  output:
  set sample_id, "${sample_id}_nr_of_seq_reads.txt" into r_nr_seq

  script:
  """
  zcat reads_1.fq.gz > seq.fq
  cat seq.fq | grep "@" | wc -l > ${sample_id}_nr_of_seq_reads.txt
  rm seq.fq
  """
}

r_nr_seq.into{
 r_nr_seq_bwa_in;
 r_nr_seq_blast_in
}

/**************************************************************************************************
 READS MAGICBLAST
**************************************************************************************************/

process reads_magicblast{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/magicblast", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from r_blast_in

  output:
  set sample_id, "${sample_id}_mapped_result_reads.magicblast" into r_blast_out

  script:
  """
  zcat reads_1.fq.gz > seq1.fq
  zcat reads_2.fq.gz > seq2.fq
  magicblast -db picornavirus_lineage -query seq1.fq -query_mate seq2.fq -infmt fastq -outfmt tabular -no_unaligned -num_threads ${task.cpus} -splice false -score 35 > ${sample_id}_mapped_result_reads.magicblast
  rm seq1.fq seq2.fq

  """
}

r_blast_out.into{
 r_blast_counter_in;
 r_blast_selection_in
}

process reads_magicblast_counter{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/magicblast", mode:'link'

  input:
  set sample_id, mblast from r_blast_counter_in

  output:
  set sample_id, "${sample_id}_nr_of_seq_reads_magicblast.txt" into r_blast_counter_out

  script:
  """
    cat ${mblast} | awk '{print \$1}' | awk '!seen[\$0]++' | wc -l > ${sample_id}_nr_of_seq_reads_magicblast.txt
  """
}

process reads_magicblast_selection{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/magicblast", mode:'link'

  input:
  set sample_id, blast from r_blast_selection_in

  output:
  set sample_id, "${sample_id}_occurrences_of_matched_ref_magicblast.txt", "${sample_id}_selected_reads_magicblast.txt" into r_magicblast_selection_out

  script:
  """
  awk '{print \$2 "\t" \$1 "\t" \$13}' ${blast} | tail -n +2 | tail -n +2 | tail -n +2| datamash -sg 2 max 3 -f | awk '{print \$1 "\t" \$2 "\t" \$3}' > ${sample_id}_selected_reads_magicblast.txt
  awk '{print \$2 "\t" \$1 "\t" \$13}' ${blast} | tail -n +2 | tail -n +2 | tail -n +2| datamash -sg 2 max 3 -f | awk '{print \$1 "\t" \$2 "\t" \$3}' | datamash -g 1 count 2 > ${sample_id}_occurrences_of_matched_ref_magicblast.txt
  """
}

combine_2_magicblast_channels = r_magicblast_selection_out.combine(r_blast_counter_out, by: 0)
combine_3_magicblast_channels = combine_2_magicblast_channels.combine(r_nr_seq_blast_in, by: 0)
//combine_3_magicblast_channels.println()


process reads_magicblast_final_sampleResults{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/magicblast", mode:'link'

  input:
  set sample_id, occur_of_matched_ref, selected_reads, nr_of_tot_reads_mapped, nr_of_tot_reads_fastq from combine_3_magicblast_channels

  output:
  set sample_id, "${sample_id}_sample_results_magicblast_reads.txt" into r_magicblast_sample_results

  script:
  """
  #!/bin/bash
  
  touch ${sample_id}_sample_results_magicblast_reads.txt
  
  nr_of_tot_fq=\$(head -n 1 ${nr_of_tot_reads_fastq})
  nr_of_tot_map=\$(head -n 1 ${nr_of_tot_reads_mapped})

  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_occur=\$(echo \$line |awk '{print \$2}')
		ref_name=\$(blastdbcmd -entry \$ref_id -db picornavirus_lineage -range 1-1 | tr -d '\n' | sed 's/.\$//')
		ref_name=\${ref_name// /_}
		ref_name=\${ref_name//>/}
		ref_name=\${ref_name//:1-1_/}
		ref_name=\${ref_name//\$ref_id/}
		ref_occur_perTTV=\$(echo -e \$ref_occur/\$nr_of_tot_map "\n" scale=2 | bc -l)
		ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_of_tot_fq "\n" scale=2 | bc -l)
		TTV_perALL=\$(echo -e \$nr_of_tot_map/\$nr_of_tot_fq "\n" scale=2 | bc -l)
		echo -e "${sample_id}\t\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$nr_of_tot_map\t\$nr_of_tot_fq" >> ${sample_id}_sample_results_magicblast_reads.txt
  done < ${occur_of_matched_ref}
  """
}


r_magicblast_sample_results.into{
 r_magicblast_sample_results1;
 r_magicblast_sample_results2
}

process reads_magicblast_final_Results{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/magicblast", mode:'link'

  input:
  set sample_id, "all_res" from r_magicblast_sample_results1

  output:
  file "${sample_id}_final_results_magicblast_reads.tsv"  into r_magicblast_final_sample

  script:
  """
  #!/bin/bash
  cat ${all_res} > tmp.txt
  awk '{print \$1 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1 sum 3 -f | datamash sum 4,5 > max_values_TTV_and_allSeq.txt
  TTV_occur=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$1}')
  nr_reads_fastq=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$2}')
  rm max_values_TTV_and_allSeq.txt
  touch All_final_results_magicblast_reads.tsv
  echo -e "Ref_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_reads" >> ${sample_id}_final_results_magicblast_reads.tsv
  cat tmp.txt | awk '{print \$2 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' | datamash -sg 1,2 sum 3,4,5 > ${sample_id}_results_summed_colmn.2.3.7.8.9.magicblast.txt
  rm tmp.txt
  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_name=\$(echo \$line |awk '{print \$2}')
        ref_occur=\$(echo \$line |awk '{print \$3}')
        ref_occur_perTTV=\$(echo -e \$ref_occur/\$TTV_occur "\n" scale=2 | bc -l)
        ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        TTV_perALL=\$(echo -e \$TTV_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        echo -e "\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$TTV_occur\t\$nr_reads_fastq" >> ${sample_id}_final_results_magicblast_reads.tsv
  done < ${sample_id}_results_summed_colmn.2.3.7.8.9.magicblast.txt
  """
}

process reads_magicblast_concat_allResults{
  tag {"All_samples"}

  publishDir "${params.publish_base_dir}/all/reads/magicblast", mode:'link'

  input:
  file "sample_res" from r_magicblast_sample_results2.map{it[1]}.collect()

  output:
  file "All_results_combined_magicblast_reads.tsv" into r_magicblast_sample_results_concat

  script:
  """
  #!/bin/bash
  touch All_results_combined_magicblast_reads.tsv
  echo -e "Sample_id\tRef_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_reads" >> All_results_combined_magicblast_reads.tsv
  for sample_file in ${sample_res}
  do
        cat \$sample_file >> All_results_combined_magicblast_reads.tsv
  done
  """
}

process reads_magicblast_final_allResults{
  tag {"All_samples"}

  publishDir "${params.publish_base_dir}/all/reads/magicblast", mode:'link'

  input:
  file "all_res" from r_magicblast_sample_results_concat

  output:
  file "All_final_results_magicblast_reads.tsv"  into r_magicblast_final_all

  script:
  """
  #!/bin/bash
  tail -n +2 ${all_res} > tmp.txt
  awk '{print \$1 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1 sum 3 -f | datamash sum 4,5 > max_values_TTV_and_allSeq.txt
  TTV_occur=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$1}')
  nr_reads_fastq=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$2}')
  rm max_values_TTV_and_allSeq.txt
  touch All_final_results_magicblast_reads.tsv
  echo -e "Ref_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_reads" >> All_final_results_magicblast_reads.tsv
  awk '{print \$2 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1,2 sum 3,4,5 > All_results_summed_colmn.2.3.7.8.9.magicblast.txt
  rm tmp.txt
  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_name=\$(echo \$line |awk '{print \$2}')
        ref_occur=\$(echo \$line |awk '{print \$3}')
        ref_occur_perTTV=\$(echo -e \$ref_occur/\$TTV_occur "\n" scale=2 | bc -l)
        ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        TTV_perALL=\$(echo -e \$TTV_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        echo -e "\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$TTV_occur\t\$nr_reads_fastq" >> All_final_results_magicblast_reads.tsv
  done < All_results_summed_colmn.2.3.7.8.9.magicblast.txt
  """
}

/**************************************************************************************************
 READS BWA
**************************************************************************************************/

process reads_bwa{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/BWA", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from r_BWA_in

  output:
  set sample_id, "${sample_id}_reads_bwaSamToBed.txt" into r_bwaBed_out

  script:
  """
  zcat reads_1.fq.gz > seq1.fq
  zcat reads_2.fq.gz > seq2.fq
  bwa mem -t ${task.cpus} ${params.bwa_idx}/picornavirus_subtree_lineage_sequences.fasta seq1.fq seq2.fq > ${sample_id}_reads_bwa_aln-pe.sam
  bedtools bamtobed -i ${sample_id}_reads_bwa_aln-pe.sam > ${sample_id}_reads_bwaSamToBed.txt
  gzip ${sample_id}_reads_bwa_aln-pe.sam
  """
}

r_bwaBed_out.into{
 r_bwaBed_counter_in;
 r_bwaBed_selection_in
}

process reads_bwa_counter{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/BWA", mode:'link'

  input:
  set sample_id, bwaSamToBed from r_bwaBed_counter_in

  output:
  set sample_id, "${sample_id}_nr_of_seq_reads_bwaBed.txt" into r_bwaBed_counter_out

  script:
  """
    cat ${bwaSamToBed} | awk '{print \$4 "\t" \$5}' | awk '(\$2 > 24)' | awk '{print \$1}' | sed 's/.\$//' | awk '!seen[\$0]++' | wc -l > ${sample_id}_nr_of_seq_reads_bwaBed.txt
  """
}

process reads_bwa_selection{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/BWA", mode:'link'

  input:
  set sample_id, bwaSamToBed from r_bwaBed_selection_in

  output:
  set sample_id, "${sample_id}_occurrences_of_matched_ref_bwaBed.txt", "${sample_id}_selected_reads_bwaBed.txt" into r_bwaBed_selection_out

  script:
  """
    awk '{print \$1 "\t" \$5 "\t" \$4}' ${bwaSamToBed} | sed 's/.\$//' | awk '!seen[\$0]++' | datamash -W -g 3 max 2 -f | awk '{print \$1 "\t" \$3 "\t" \$4}' | awk '(\$3 > 24)' > ${sample_id}_selected_reads_bwaBed.txt
    awk '{print \$1 "\t" \$5 "\t" \$4}' ${bwaSamToBed} | sed 's/.\$//' | awk '!seen[\$0]++' | datamash -W -g 3 max 2 -f | awk '{print \$1 "\t" \$3 "\t" \$4}' | awk '(\$3 > 24)' | datamash -sg 1 count 1 > ${sample_id}_occurrences_of_matched_ref_bwaBed.txt
  """
}

combine_2_bwa_channels = r_bwaBed_selection_out.combine(r_bwaBed_counter_out, by: 0)
combine_3_bwa_channels = combine_2_bwa_channels.combine(r_nr_seq_bwa_in, by: 0)
//combine_3_bwa_channels.println()


process reads_bwa_final_sampleResults{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/BWA", mode:'link'

  input:
  set sample_id, occur_of_matched_ref, selected_reads, nr_of_tot_reads_mapped, nr_of_tot_reads_fastq from combine_3_bwa_channels

  output:
  set sample_id, "${sample_id}_sample_results_bwaBed_reads.txt" into r_bwaBed_sample_results

  script:
  """
  #!/bin/bash
  
  touch ${sample_id}_sample_results_bwaBed_reads.txt
  
  nr_of_tot_fq=\$(head -n 1 ${nr_of_tot_reads_fastq})
  nr_of_tot_map=\$(head -n 1 ${nr_of_tot_reads_mapped})

  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_occur=\$(echo \$line |awk '{print \$2}')
		ref_name=\$(blastdbcmd -entry \$ref_id -db picornavirus_lineage -range 1-1 | tr -d '\n' | sed 's/.\$//')
		ref_name=\${ref_name// /_}
		ref_name=\${ref_name//>/}
		ref_name=\${ref_name//:1-1_/}
		ref_name=\${ref_name//\$ref_id/}
		ref_occur_perTTV=\$(echo -e \$ref_occur/\$nr_of_tot_map "\n" scale=2 | bc -l)
		ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_of_tot_fq "\n" scale=2 | bc -l)
		TTV_perALL=\$(echo -e \$nr_of_tot_map/\$nr_of_tot_fq "\n" scale=2 | bc -l)
		echo -e "${sample_id}\t\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$nr_of_tot_map\t\$nr_of_tot_fq" >> ${sample_id}_sample_results_bwaBed_reads.txt
  done < ${occur_of_matched_ref}
  """
}

r_bwaBed_sample_results.into{
 r_bwaBed_sample_results1;
 r_bwaBed_sample_results2
}

process reads_bwa_final_Results{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads/BWA", mode:'link'

  input:
  set sample_id, "all_res" from r_bwaBed_sample_results1

  output:
  file "${sample_id}_final_results_bwaBed_reads.tsv"  into r_bwaBed_final_sample

  script:
  """
  #!/bin/bash
  cat ${all_res} > tmp.txt
  awk '{print \$1 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1 sum 3 -f | datamash sum 4,5 > max_values_TTV_and_allSeq.txt
  TTV_occur=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$1}')
  nr_reads_fastq=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$2}')
  rm max_values_TTV_and_allSeq.txt
  touch All_final_results_bwaBed_reads.tsv
  echo -e "Ref_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_reads" >> ${sample_id}_final_results_bwaBed_reads.tsv
  cat tmp.txt | awk '{print \$2 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' | datamash -sg 1,2 sum 3,4,5 > ${sample_id}_results_summed_colmn.2.3.7.8.9.bwaBed.txt
  rm tmp.txt
  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_name=\$(echo \$line |awk '{print \$2}')
        ref_occur=\$(echo \$line |awk '{print \$3}')
        ref_occur_perTTV=\$(echo -e \$ref_occur/\$TTV_occur "\n" scale=2 | bc -l)
        ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        TTV_perALL=\$(echo -e \$TTV_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        echo -e "\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$TTV_occur\t\$nr_reads_fastq" >> ${sample_id}_final_results_bwaBed_reads.tsv
  done < ${sample_id}_results_summed_colmn.2.3.7.8.9.bwaBed.txt
  """
}


process reads_bwa_concat_allResults{
  tag {"All_samples"}

  publishDir "${params.publish_base_dir}/all/reads/BWA", mode:'link'

  input:
  file "sample_res" from r_bwaBed_sample_results2.map{it[1]}.collect()

  output:
  file "All_results_combined_bwaBed_reads.tsv" into r_bwaBed_sample_results_concat

  script:
  """
  #!/bin/bash
  touch All_results_combined_bwaBed_reads.tsv
  echo -e "Sample_id\tRef_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_reads" >> All_results_combined_bwaBed_reads.tsv
  for sample_file in ${sample_res}
  do
        cat \$sample_file >> All_results_combined_bwaBed_reads.tsv
  done
  """
}

process reads_bwa_final_allResults{
  tag {"All_samples"}

  publishDir "${params.publish_base_dir}/all/reads/BWA", mode:'link'

  input:
  file "all_res" from r_bwaBed_sample_results_concat

  output:
  file "All_final_results_bwaBed_reads.tsv"  into r_bwaBed_final_all

  script:
  """
  #!/bin/bash
  tail -n +2 ${all_res} > tmp.txt
  awk '{print \$1 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1 sum 3 -f | datamash sum 4,5 > max_values_TTV_and_allSeq.txt
  TTV_occur=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$1}')
  nr_reads_fastq=\$(head -n 1 max_values_TTV_and_allSeq.txt | awk '{print \$2}')
  rm max_values_TTV_and_allSeq.txt
  touch All_final_results_bwaBed_reads.tsv
  echo -e "Ref_id\tRef_name\tRef_occur/TTV_occur\tRef_occur/ALL\tTTV/ALL\tRef_occur\tTTV_occur\tAll_reads" >> All_final_results_bwaBed_reads.tsv
  awk '{print \$2 "\t" \$3 "\t" \$7 "\t" \$8 "\t" \$9}' tmp.txt | datamash -sg 1,2 sum 3,4,5 > All_results_summed_colmn.2.3.7.8.9.bwaBed.txt
  rm tmp.txt
  while read -r line; do
        ref_id=\$(echo \$line |awk '{print \$1}')
        ref_name=\$(echo \$line |awk '{print \$2}')
        ref_occur=\$(echo \$line |awk '{print \$3}')
        ref_occur_perTTV=\$(echo -e \$ref_occur/\$TTV_occur "\n" scale=2 | bc -l)
        ref_occur_perALL=\$(echo -e \$ref_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        TTV_perALL=\$(echo -e \$TTV_occur/\$nr_reads_fastq "\n" scale=2 | bc -l)
        echo -e "\$ref_id\t\$ref_name\t\$ref_occur_perTTV\t\$ref_occur_perALL\t\$TTV_perALL\t\$ref_occur\t\$TTV_occur\t\$nr_reads_fastq" >> All_final_results_bwaBed_reads.tsv
  done < All_results_summed_colmn.2.3.7.8.9.bwaBed.txt
  """
}




