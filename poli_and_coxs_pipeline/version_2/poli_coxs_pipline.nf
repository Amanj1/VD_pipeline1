#!/usr/bin/env nextflow

params.reads_dir='preprocessing'

fq_reads = Channel.fromFilePairs("${params.reads_dir}/*/*{1,2}.fq.gz",size: 2){ (it.name =~ /P[0-9]{3,5}_[0-9]{3,5}/)[0]}


//fq_reads.println()
//fasta_contigs.println()



fq_reads.into{
    BWA_in;
	blastn1_in;
	blastn2_in;
	extract_reads
}


process BWA_reads{
  tag{"${sample_id}"}
  publishDir "${params.publish_base_dir}/${sample_id}/BWA_reads/", mode:'link'
   
  input:
  set sample_id, 'reads_*.fq.gz' from BWA_in

  output:
  set sample_id, "${sample_id}_reads_bwaSamToBed.txt", "${sample_id}_reads_bwa_aln-pe.sam.gz" into bwa_out

  script:

  """
 zcat reads_1.fq.gz > read1.fq
 zcat reads_2.fq.gz > read2.fq
 bwa mem -t "${task.cpus}" \$bwaIndex/ref_contigs.fa read1.fq read2.fq > ${sample_id}_reads_bwa_aln-pe.sam
 bedtools bamtobed -i ${sample_id}_reads_bwa_aln-pe.sam > ${sample_id}_reads_bwaSamToBed.txt
 gzip ${sample_id}_reads_bwa_aln-pe.sam
 rm read1.fq read2.fq
  """
}

bwa_out.into{
    bwa_extraction_data;
	bwa_collection
}

process blastn_read1{
  tag{"${sample_id}"}
  publishDir "${params.publish_base_dir}/${sample_id}/blastn_reads/", mode:'link'
   
  input:
  set sample_id, 'reads_*.fq.gz' from blastn1_in

  output:
  set sample_id, "${sample_id}_read1.blastnTabular" into blastn1_out

  script:

  """
  zcat reads_1.fq.gz > read1.fq
  sed -n '1~4s/^@/>/p;2~4p' read1.fq > tmp.fa
  blastn -db ref_contigs -query tmp.fa -num_threads "${task.cpus}" -outfmt 6 >  ${sample_id}_read1.blastnTabular
  rm tmp.fa read1.fq
  """
}

blastn1_out.into{
	blastn1_extraction_data;
	blastn1_collection
}

process blastn_read2{
  tag{"${sample_id}"}
  publishDir "${params.publish_base_dir}/${sample_id}/blastn_reads/", mode:'link'
   
  input:
  set sample_id, 'reads_*.fq.gz' from blastn2_in

  output:
  set sample_id, "${sample_id}_read2.blastnTabular" into blastn2_out

  script:

  """
  zcat reads_2.fq.gz > read2.fq
  sed -n '1~4s/^@/>/p;2~4p' read2.fq > tmp.fa
  blastn -db ref_contigs -query tmp.fa -num_threads "${task.cpus}" -outfmt 6 >  ${sample_id}_read2.blastnTabular
  rm tmp.fa read2.fq
  """
}

blastn2_out.into{
	blastn2_extraction_data;
	blastn2_collection
}

combine_blastn1_with_blastn2 = blastn1_extraction_data.combine(blastn2_extraction_data, by: 0)
all_methods_results = bwa_extraction_data.combine(combine_blastn1_with_blastn2, by: 0)
all_extract_reads = all_methods_results.combine(extract_reads, by: 0)


//combine_BWA_extraction.println()
//combine_blastn_extraction.println()
//all_extract_reads.println()


process Reads_extraction{
  tag{"${sample_id}"}
  publishDir "${params.publish_base_dir}/${sample_id}/extracted_reads/fastq/", mode:'link'
   
  input:
  set sample_id, bwaBed, bwaSam, blastn_r1, blastn_r2, 'reads_*.fq.gz' from all_extract_reads

  output:
  set sample_id, "${sample_id}_extracted_reads_1.fq.gz", "${sample_id}_extracted_reads_2.fq.gz" into reads_extraction_out

  script:

  """
  zcat reads_1.fq.gz > read1.fq
  zcat reads_2.fq.gz > read2.fq
  cat ${blastn_r1} ${blastn_r2} | awk '{print \$1}' | awk '!seen[\$0]++' > blastn_extracted_reads_list.txt
  awk '{print \$4}' ${bwaBed} | sed 's/.\$//' | sed 's/.\$//' | awk '!seen[\$0]++' > bwa_extracted_reads_list.txt
  cat blastn_extracted_reads_list bwa_extracted_reads_list | awk '!seen[\$0]++' > extracted_reads_list.txt  
  seqtk subseq read1.fq extracted_reads_list.txt > ${sample_id}_extracted_reads_1.fq
  seqtk subseq read2.fq extracted_reads_list.txt > ${sample_id}_extracted_reads_2.fq
  gzip ${sample_id}_extracted_reads_1.fq
  gzip ${sample_id}_extracted_reads_2.fq
  rm extracted_reads_list.txt read1.fq read2.fq blastn_extracted_reads_list.txt bwa_extracted_reads_list.txt
  """
}

process Extracted_reads_fq_to_fa{
  tag{"${sample_id}"}
  publishDir "${params.publish_base_dir}/${sample_id}/extracted_reads/fasta/", mode:'link'
   
  input:
  set sample_id, read1, read2 from reads_extraction_out

  output:
  set sample_id, "${sample_id}_extracted_reads_1.fa.gz", "${sample_id}_extracted_reads_2.fa.gz" into reads_fqTofa

  script:

  """
  zcat ${read1} > read1.fq
  zcat ${read2} > read2.fq
  #Convert ILLUMINA 1.3+ FASTQ to FASTA and mask bases with quality lower than 20 to lowercases
  seqtk seq -aQ64 -q20 read1.fq > ${sample_id}_extracted_reads_1.fa
  seqtk seq -aQ64 -q20 read2.fq > ${sample_id}_extracted_reads_2.fa
  gzip ${sample_id}_extracted_reads_1.fa
  gzip ${sample_id}_extracted_reads_2.fa
  rm read1.fq read2.fq
  """
}



