Viral results consist of:
	Filtered results:
		Viruses
		Prokaryotic phages
		Other sequences
		Unclassified sequences:
			Some unclassified sequences could be prokaryotic cells
	
All results consist of:
	All possible results:
		Viruses
		Prokaryotic phages
		Other
		Unclassified
		Host contamination (sequences that were not removed in preprocessing pipeline)
		Eukaryotic species/cells
		Prokaryotic cells

Unknown sequences consist of:
	Everything that does not have any taxonomic information:
		Name: Unknown
		Division: Unknown
		Rank: Unknown
	TSV file:
		contains additional column with nucleotide sequences of contigs
	Html file:
		Does not contain column with nucleotide sequences

Methods used:
	Reads (only nucleotide based approach):
		Kraken2
		Metaphlan2
		FastViromeExplorer
	Contigs (Both nucleotide and protein based approach):
		Kraken2 (nucleotide based approach)
		Diamond blastx (nucleotide translated into protein sequences based approach)
			Everything in the tables with an accession number is retrieved from Diamond blastx
	HTML table (Internet connection required to fetch methods used in table format and display results):
		Search box with filtering based on query input from user
		Clickable headers for sorting based on values from selected column
		Show: select number of results to display on each page
	TSV table:
		Can be opend in spreadsheet software (excel) 
	JSON file:
		Web-based reading of file

Minimum reads lengths:
	reads 2x150bp (PE)
	contigs 500bp (Filtered in discovery pipeline)
