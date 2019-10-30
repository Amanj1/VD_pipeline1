#!/bin/bash
sample_id=$2
assembler=$3

cut -f 2 $1| awk '!seen[$0]++' > tmp.txt
count=1
echo "Accession_nr	Title	Organism" > tmp.tsv
echo "Accession_numbers_not_found_in_nr_db" > acc_nr_not_found.txt
echo "Existing acession number for esearch" > acc_nr_exist_for_esearch.txt
for line in $(cat tmp.txt)
do
	title_blast=$(blastdbcmd -entry $line -db nr -range 1-1)
	title=$(echo $title_blast|awk -F $line '{print $2}'| awk -F '>' '{print $1}')
	if [ -z "$title" ]
	then
		if [ -z "$title_blast" ]
		then
			echo "$line" >> acc_nr_not_found.txt
		else
			echo "$line" >> acc_nr_exist_for_esearch.txt
		fi
	else
		organism=$(echo "$title"|awk -F '[' '{print $2}'|awk -F ']' '{print $1}')
		organism=${organism// /_}
		title=${title// /_}
		echo $line	$title	$organism >>tmp.tsv
	fi
done
tr ' ' '\t' <tmp.tsv >"${sample_id}_${assembler}_contigs_diamond_organism_names.tsv"
rm tmp.tsv
mv acc_nr_not_found.txt ${sample_id}_${assembler}_diamond_acc_nr_not_found.txt
mv acc_nr_exist_for_esearch.txt ${sample_id}_${assembler}_diamond_acc_nr_exist_for_esearch.txt
