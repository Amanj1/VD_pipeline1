# New VD pipeline - files

Each sample will have two folders

* reads - read-based analysis
* megahit -  assembly-based analysis (megahit is the assembler's name)

## Read-based analysis

3 tools are used for the read-based analysis

### Metaphlan2

Tool from Segata lab, nucleotide-based, that characterizes microbiome based on marker genes. Not very sensitive, but it is the least prone to False Positives

Outputs:

* **Microbial profile**: `<sample_name>_metaphlan2.tsv`.

  Relative abundance normalized to 100% for each microbe


### Kraken2:

Tool from John Hopkins, nucleotide based. More sensitive to false positives (reports reads mapped to each genome). Reads were mapped to the default database with Bacteria, Virus and Fungi.


Outputs:

* **Microbial profile**: `<sample_name>_kraken2_report.txt`

  Columns:
  * Percent of total sequences covered by this taxon(this level + all children)
  * Number of sequences covered by this taxon(this level + all children)
  * Number of sequences *assigned directly* to this taxon
  * Rank code:(U)nclassified (D)omain, (K)ingdom, ... etc
  * NCBI taxonomic id
  * Scientiifc name

* **Read classification**: `<sample_name>_kraken2.txt`

  This file is really big, specifies how each read was classified. Not so important

### FastViromeExplorer

A virus-specific mapping tool.

Outputs:

* **Virome profile**: `<sample_name>_fastviromeexplorer_abundance.tsv`

* **Read classification**: `<sample_name>_fastviromeexplorer.sam`

## Assembly-based analysis (Megahit)

Contents in the `megahit/` folder

### Assembly results

We create longer contigs using an assembler called Megahit.
* The raw results of the assembly are in `1_assembly/final.contigs.fa` (all contigs)
* Contigs greater than 500bp are in `2_filt_contigs/contigs_filt.fa`

Some QC data is available in `2_filt_contigs/<sample_name>_flagstat.txt`, which tells what proportion of the reads actually map to the contigs

The analysis was performed with

### Gene prediction - FragGeneScan

The results are in the `orfs/` folder. Files are

* `<sample_id>_megahit_fgs_orfs.faa`: Aminoacid sequence of predicted gene seqs
* `<sample_id>_megahit_fgs_orfs.fna`: Nucleotide sequence of predicted gene seqs
* `<sample_id>_megahit_fgs_orfs.gff` GFF file with the contig and coordinates of where the gene was predicted
* `<sample_id>_megahit_fgs_orfs.gff`: Fasta like file with the contig and coordinates of where the genes were predicted

### Classification - Kraken

The Microbial profile and read classification outputs are the same as in the reads. However, the Contig classification file is more relevant because you can see the direct classification of each contig.

Also there is an extra file with the unmapped contigs to kraken2.

* **Microbial profile**: `<sample_name>_megahit_kraken2_report.txt`
* **Contig classification**: `<sample_name>_megahit_kraken2.txt`
* **Unmapped contigs**: `<sample_name>_megahit_kraken2_unmapped.fa`


### Classification - Diamond

Diamond is a blastx replacement. Performs translated search against a protein database (in this case, NR). Diamond does not issue a microbial profile, just maps the contigs to the database, and returns the reads.

Outputs

* **Contig alignments**: `<sample_name>_megahit_diamond.tsv`

Column headers are in `diamond_blast_cols.txt`


* **Unmapped contigs**: `<sample_name>_megahit_diamond_unmapped.fa`

* **Irrelevant**: `<sample_name>_megahit_diamond.daa`

### Classification - Virfinder

It is a classification tool that gives a score of how likely a contig has a viral origin

Output

* **Contig viral scores**: `<sample_name>_megahit_virfinder.csv`
  Columns:
  * Contig name
  * Contig length
  * Score: Between 0-1. The higher(closer to 1), the more likely it is a virus
  * p-value: Of how diffferent is the sequence from prokaryotic hosts
