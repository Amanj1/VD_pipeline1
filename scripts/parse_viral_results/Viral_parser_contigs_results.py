#!/usr/bin/python
import sys
import os
import glob
import re
import string
import csv
import json
import collections
import shutil

# aliases
OrderedDict = collections.OrderedDict
# Pathway to current dir
cur_dir=os.getcwd() 
### Boolean for metadata 1 = true or 0 = false
Bool_metadata = 1

###FILTER ALL results or ONLY VIRUS results, 1 = true or 0 = false. True for all results
Bool_filter = 0

#### PATHWAYS ####
###File pathway for creating table
files_kraken2 = glob.glob('/home/amanj/Documents/Results_from_pipeline/Table_results_from_pipeline/Table_contigs/creating_script_for_contigs/discovery/P*/megahit/P*_kraken2_report.txt')
files_diamond_org = glob.glob('/home/amanj/Documents/Results_from_pipeline/Table_results_from_pipeline/Table_contigs/creating_script_for_contigs/discovery/P*/megahit/P*_full_list.tsv')
files_diamond_tax_data = glob.glob('/home/amanj/Documents/Results_from_pipeline/Table_results_from_pipeline/Table_contigs/creating_script_for_contigs/discovery/P*/megahit/P*_taxonomic_info.tsv')
files_diamond_evalue = glob.glob('/home/amanj/Documents/Results_from_pipeline/Table_results_from_pipeline/Table_contigs/creating_script_for_contigs/discovery/P*/megahit/P*_megahit_diamond.tsv')

if Bool_metadata == 1:
        ### File pathway for adding metadata
        #File pathway for 2nd tsv table with metadata
        SndTsv = glob.glob('/home/amanj/Documents/190710/file_for_script/TTV_miseq_metadata.txt')
        #Files for human material
        files_human = glob.glob('/home/amanj/Documents/190710/file_for_script/P*/humanrm_flagstat/*.flagstat')
        #Files for total sequences
        files_totalSeq = glob.glob('/home/amanj/Documents/190710/file_for_script/P*/raw_fastqc/P*/fastqc_data.txt')

def replaceOldDB(filename, container): ###Function for Methaphlan2 container
        skip = 0 #Skipping header
        d = {} #Dictonary container
        with open(filename) as f:
                for line in f:
                        (key, val) = line.split()
                        d[str(key)] = str(val)
        f.close()
        for c in container:
                if skip != 0 and c[1] in d:
                        c[1] = d[str(c[1])]
                skip = 1 # value after skipping header
        return container

def taxonLevelString(arr, select):
        strn = ""
        tmp = []
        #Select 1 for kraken and select 2 for methaphlan
        if select == 1:
                for a in arr:
                        if strn == "":
                                strn = strn + a
                        else:
                                strn = strn + "->" + a
        if select == 2:
                for a in arr:
                        if "f__" in a:
                                tmp.append(a.replace("f__", ""))
                        elif "g__" in a:
                                tmp.append(a.replace("g__", ""))
                        elif "s__" in a:
                                tmp.append(a.replace("s__", ""))
                        elif "t__" in a:
                                tmp.append(a.replace("t__", ""))
                for a in tmp:
                        if strn == "":
                                strn = strn + a
                        else:
                                strn = strn + "->" + a
        return strn
                
def CreateFile(fileName):
        name = "Parse_viral_results_output/Combined_Results_"+fileName+".txt"
        f= open(name,"w+")
        f.close() 
        return None

def WriteToFile(resultArray, fileName):
        name = "Combined_Results_"+fileName+".txt"
        tmp = ""
        f=open(name, "a+")
        for r in resultArray:
                tmp = str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t"+str(r[3])+"\t"+str(r[4])+"\t"+str(r[5])+"\t"+str(r[6])+"\t"+str(r[7])+"\t"+str(r[8])+"\n"
                f.write(tmp)
        f.close() 
        return None

def ireplace(old, new, text):
        index_l = text.lower().index(old.lower())
        return text[:index_l] + new + text[index_l + len(old):]

def parseKraken2(fileInput):
        """
        Columns in each row:
        Column 1: percentage of reads in the clade/taxon in Column 6
        Column 2: number of reads in the clade.
                Column 3: number of reads in the clade but not further classified.
        Column 4: code indicating the rank of the classification: 
                                                                        (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, (S)pecies).
        Column 5: NCBI taxonomy ID.
        
        Example run:
                ['35.82', '1461', '1461', 'U', '0', 'unclassified']
                ['64.18', '2618', '35', 'R', '1', 'root']
                ['59.08', '2410', '11', 'R1', '131567', 'cellular', 'organisms']
                ['58.69', '2394', '176', 'D', '2', 'Bacteria']
                ['51.61', '2105', '85', 'P', '1224', 'Proteobacteria']
        
        Source: http://sepsis-omics.github.io/tutorials/modules/kraken/
        """
        #Filter bacteria from results
        check_level = 1 # If level =1 stop appending to container
        
        addClasses = ['F','G','S']
        ClassLabel = ['Family', 'Genus', 'Species', 'Sub_Species', 'Sub_Genus', 'Sub_Family']
        #taxonLevel [Classification Char]
        taxonLevel = []
        #taxonName [Classification Level Name]
        taxonName = []
        
        container = []
        tmp = []
        tmpstr = ""
        run = 1
        f = open(fileInput)
        while run != 0:
                line = f.readline().split()
                if len(line) == 0:
                        break
                """
                FIX pattern some results have G1, S1,S2,S3 instead of G and S only
                regular expression match 
                pattern = re.compile("^([A-Z][0-9]+)+$")
                pattern.match(string)
                """
                pattern = re.compile("D\Z")
                if pattern.match(str(line[3])):
                        if "virus" not in str(line[5].lower()) or "phage" not in str(line[5].lower()):
                                if Bool_filter == 0:
                                    check_level = 1
                                else:
                                    check_level = 0
                        if "virus" in str(line[5].lower()) or "phage" in str(line[5].lower()):
                                check_level = 0
                if check_level != 1:    
                        if any(str(it) in str(line[3]) for it in addClasses):
                                taxonLevel.append(line[3])
                                tmp = line[5:len(line)]
                                line = (line[0:5])
                                for n in tmp:
                                        if tmpstr == "":
                                                tmpstr = tmpstr + n
                                        else:
                                                tmpstr = tmpstr + '_' + n 
                                line.insert(0, tmpstr)
                                taxonName.append(tmpstr)
                                tmpstr = ""
                                for n in range(len(addClasses)):
                                        pattern = re.compile("S[0-9]")
                                        if pattern.match(str(line[4])):
                                                line.insert(0, ClassLabel[3])
                                                container.append(line)
                                                break
                                        if str(line[4]) == str(addClasses[n]):
                                                line.insert(0, ClassLabel[n])
                                                container.append(line)
                                                break              
                                #print(line)
        ###################################################
        tmpArr = []
        ListofArr = []
        tmpArrChar = []
        ListofArrChar = []
        for i in range(len(taxonLevel)):
                if i == len(taxonLevel)-1:
                        ListofArr.append(tmpArr)
                        ListofArrChar.append(tmpArrChar)
                if tmpArr == []:
                        tmpArr.append(taxonName[i])
                        tmpArrChar.append(taxonLevel[i])
                elif taxonLevel[i] == "F":
                        ListofArr.append(tmpArr)
                        tmpArr = []
                        ListofArrChar.append(tmpArrChar)
                        tmpArrChar = []
                        tmpArr.append(taxonName[i])
                        tmpArrChar.append(taxonLevel[i])
                else:
                        tmpArr.append(taxonName[i])
                        tmpArrChar.append(taxonLevel[i])                
        #########################################################
        tmpArr = []
        tmpArrChar = []
        check_level = 0
        for c in container:
                for i in range(len(ListofArrChar)):
                        if c[1] in ListofArr[i] and check_level == 0:
                                c.append(ListofArrChar[i])
                                c.append(ListofArr[i])
                                check_level = 1
                check_level = 0
        Genus = ""
        Species = ""
        ######################################################
        for c in container:
                tmpArr = []
                check_level = 0
                for i in range(len(c[7])):
                        if c[5] == "F":
                                tmpArr.append(c[1])
                                c[8] = taxonLevelString(tmpArr, 1)
                                c.remove(c[7])
                                break
                        if c[5] == "G":
                                if c[7][i] == "F":
                                        tmpArr.append(c[8][i])
                                        tmpArr.append(c[1])
                                        c[8] = taxonLevelString(tmpArr, 1)
                                        c.remove(c[7])
                                        break
                        pattern = re.compile("S\Z")
                        if pattern.match(c[5]):
                                check_level = 1
                                if c[7][i] == "F":
                                        tmpArr.append(c[8][i])
                                if c[7][i] == "G":
                                        Genus = c[8][i]
                                if c[8][i] == c[1]:
                                        tmpArr.append(Genus)
                                        tmpArr.append(c[1])
                                        c[8] = taxonLevelString(tmpArr, 1)
                                        c.remove(c[7])
                                        break
                        pattern = re.compile("S[0-9]")
                        if pattern.match(str(c[5])) and check_level == 0: 
                                if c[7][i] == "F":
                                        tmpArr.append(c[8][i])
                                if c[7][i] == "G":
                                        Genus = c[8][i]
                                pattern = re.compile("S\Z")
                                if pattern.match(c[7][i]):
                                        Species = c[8][i]
                                if c[8][i] == c[1]:
                                        tmpArr.append(Genus)
                                        tmpArr.append(Species)
                                        tmpArr.append(c[1])
                                        c[8] = taxonLevelString(tmpArr, 1)
                                        c.remove(c[7])
                                        break
        return container

def parseDiamond(fileOrg,fileTax,fileEvalue):
    container = []
    tmp = []
    unknownSeq = []
    tmpstr = ""
    run=1
    run2=1
    run3=1
    f = open(fileOrg)
    line = f.readline()
    while run != 0: #Find Organism name
        line = f.readline().split()
        tmp = []
        if len(line) == 2:
            line.append("unclassified")
        if not line or line is None:
            #print('Finished')
            run=0        
            run2=0        
            run3=0
        else:
            line[1] = line[1].replace(':1-1_','')            
            #print(line)
        f2=open(fileTax)
        line2 = f2.readline()
        while run2 != 0: # Find divison and rank
            line2 = f2.readline().split()
            if not line2 or line is None:
                #print("Unknown in run2")
                #print(line)
                #print(line2)
                #print('\n')
                unknownSeq.append(line)
                run2=0
                run3=0
            #elif line[2].lower() in line2[0].lower() or line2[0].lower() in line[2].lower():
            elif line2[0].lower() == line[2].lower():
                #print(line2)
                if len(line2) == 2:
                    line2.append("unclassified")
                tmp.append(line2[2])  # rank
                tmp.append(line2[0])  # Organism name
                tmp.append(line2[1])  # Divison
                tmp.append(line[0])   # Accession number
                tmp.append(line[1])   # Title
                f2.close()
                run2=0
        f3=open(fileEvalue)
        while run3 != 0: # Find e-value for accession number
            line3 = f3.readline().split()
            if not line3 or line is None:
                print("Unknown in run3")
                print(line)
                print(line3)
                run3=0
            elif line[0] in line3[1]:
                tmp.append(line3[9]) # pident
                tmp.append(line3[6]) # evalue
                #print(line3)
                f3.close()
                run3=0
        run2 = 1
        run3 = 1
        #print(tmp)
        #print('\n')
        #print('\n')
        # tmp = [rank, organism name, divison, accession number, title, pident, evalue]
        if len(tmp) != 0:
            container.append(tmp)
    container = container[:-1]
    #print(len(container))
    all_results = [container,unknownSeq]
    return all_results


def combineResults(Kraken, Diamond, ID):
        """
        #Create header 
        # Patient_ID    Name    Kraken  FastViromeExplorer      MetaPhlan       Classification  Full_Taxon_level
        Header = "Patient_ID\tName\tKraken_#n/home/amanj/Documents/190710
        /home/amanj/Documents/190714rReads\tFastViromeExplorer\tMetaPhlan\tClassification\tFull_Taxon_level(Family->Genus->Species->Sub_species)\n"
        CreateFile(ID)
        name = "Combined_Results_"+ID+".txt"
        f=open(name, "a+")
        f.write(Header)
        f.close() 
        """
        #tmp = [ID,"Name","Accession number","Kraken", "Diamond_p", "Diamond_e", "Class","Divison", "Level"]
        tmp = [ID,"-","-", "-", "-", "-", "-", "-", "-"]
        container = []
        check = 1
        
        for d in Diamond:
            tmp[1] = d[1] #Name
            tmp[6] = d[0] #Class
            tmp[2] = d[3] #Accession number
            tmp[7] = d[2].replace('_','') #Divison
            tmp[4] = d[5] #pident
            tmp[5] = d[6] #e-value
            for k in Kraken:
                if str(d[1].lower()) == str(k[1].lower()):
                    tmp[6] = k[0] #Class
                    tmp[3] = k[3] #Kraken nr reads
                    tmp[8] = k[7] #Taxon level from Kraken
            if Bool_filter == 1 or "virus" in str(tmp[7].lower()) or "unclassified" in str(tmp[7].lower()) or "phage" in str(tmp[7].lower()):
                container.append(tmp)
            tmp = [ID,"-","-", "-", "-", "-", "-", "-", "-"]            
    
        for k in Kraken:
            for c in container:
                if str(c[1].lower()) == str(k[1].lower()):
                    check=0
            if check == 1:
                #Name
                tmp[1] = k[1]
                #Class
                tmp[6] = k[0]
                #Kraken
                tmp[3] = k[3]
                #Taxon level from Kraken
                tmp[8] = k[7]
                container.append(tmp)
                tmp = [ID,"-","-", "-", "-", "-", "-", "-", "-"]
            check = 1
        #for c in container:
        #       print(c)
        WriteToFile(container, ID)
        return None

def CreateInitialTSVresults(files_diamond_org, files_diamond_tax, files_diamond_e, files_kraken2):
        Unknown = []
        for i in range(len(files_diamond_org)):
            fileName1 = files_diamond_org[i]
            fileName2 = files_diamond_tax[i]
            fileName3 = files_diamond_e[i]
            D_container = parseDiamond(fileName1, fileName2, fileName3)
            #Fetching Sample ID
            Sample_ID = fileName1.split('/').pop()
            count = 0
            tmp = ""
            for id in Sample_ID:
                if id == '_':
                    count=count+1
                if count < 2:
                    tmp = tmp+id
                    #print(tmp)
            Sample_ID = tmp    
            #print(Sample_ID)
            fileName = files_kraken2[i]
            K_container = parseKraken2(fileName)
            combineResults(K_container, D_container[0], Sample_ID)
            for unk_d in D_container[1]:
                unk_d.append(Sample_ID)
                Unknown.append(unk_d)
        
        all_results = glob.glob('Combined_Results_*.txt')
        all_results.sort()
        
        #Create header 
        #[ID,"Name","Accession number","Kraken", "Diamond_p", "Diamond_e", "Classification","Divison", "Full_Taxon_level"]
        Header = "sample_id\tname\taccession_number\tkraken_nr_reads\tpercentage_of_identical_matches\texpected_value\tclassification\tdivison\tfull_taxon_level:family_genus_species_subspecies\n"
        f = open("TMPresults.txt", "w")
        f.write(Header)
        for result in all_results:
                f1 = open(result, "r")
                f.write(f1.read())
                f1.close()
        f.close()
        Header = "Accession_nr\tTitle\tOrganism\tsample_id\n"
        f = open("Unknown_sequences_from_diamond.txt", "w")
        f.write(Header)
        for result in Unknown:
                f.write(str(result[0])+'\t'+str(result[1])+'\t'+str(result[2])+'\t'+str(result[3])+'\n')
        f.close()
        
        return None

def addMetaData(FstTsv, SndTsv, files_human, files_totalSeq):
        #Contianer for all new columns: List of lists
        Container_newColns = []
        #Temporary array with columns for one sample_id
        temp = []
        #Varibales for loop
        run = 1
        f = open(SndTsv[0])
        #Skipping header
        line = f.readline().split()
        
        while run != 0:
                line = f.readline().split()
                if len(line) != 0:
                        sequencing_id = line[1]
                        amplification = line[(len(line)-1)]
                        line.pop()
                        sample_id = "_".join(line[3:(len(line))])
                        sequencing_type = line[2]
                        temp = [sequencing_id, sample_id, amplification, sequencing_type]
                        Container_newColns.append(temp)
                else:
                        run = 0
        f.close()
        for i in range(len(Container_newColns)):
                files_human[i]
                f = open(files_human[i])
                lines=f.readlines()
                #adding number of reads mapped
                temp=lines[4].split()
                #Example ['5851303', '+', '0', 'mapped', '(90.43%', ':', 'N/A)'] -> '5851303'
                Container_newColns[i].append(temp[0])
        f.close()
        
        C_newCol = 0 #Counter for iteration of Container_newColns
        C_totSeq = 0 #Counter for iteration of files_totalSeq
        run = 1
        SumOfTotalSeq = 0
        check = 0
        #files_totalSeq contains 180 files and container_newcolns contains 48 lists
        # Multiple files for each sample_id
        # while loop adds upp all total sequences for each sample_id
        while run != 0:
                f = open(files_totalSeq[C_totSeq])
                lines=f.readlines()
                temp = lines[6].split()
                check = files_totalSeq[C_totSeq].find(Container_newColns[C_newCol][0])
                if check != -1:
                        SumOfTotalSeq = SumOfTotalSeq + int(temp[2])
                        C_totSeq = C_totSeq + 1
                        #print(SumOfTotalSeq)
                        if C_totSeq == len(files_totalSeq):
                                Container_newColns[C_newCol].append(str(SumOfTotalSeq))
                                run = 0
                else:
                        Container_newColns[C_newCol].append(str(SumOfTotalSeq))
                        SumOfTotalSeq = 0
                        C_newCol = C_newCol + 1
                        if C_newCol == len(Container_newColns):
                                run = 0
        f.close()
        #print(FstTsv[0])
        f = open(FstTsv[0])
        lines = f.readlines()
        #print(lines[2])
        f.close()
        #f.write('\n'.join(j + '\t'.join(i[1:])))
        with open('TMPresults_fixed_added_metadata.txt', 'w') as f: 
                tmp = lines[0].split() + ['name_of_pool', 'amplification', 'sequence_type', 'human_material', 'total_sequences']
                tmp = '\t'.join(tmp).rstrip('\n')
                #print(tmp)
                f.write(tmp.lower() + '\n')
                #print(Container_newColns)
                for i in Container_newColns:
                        #print(i[0])
                        for j in lines:
                                #print(j.split()[1])
                                if i[0] == j.split()[0]:
                                        tmp = j.rstrip(',').split() + i[1:]
                                        #print(len(i[1:]))
                                        #print(i[1:])
                                        tmp = '\t'.join(tmp).rstrip('\n')
                                        #print(tmp)
                                        #print(len(tmp.split()))
                                        f.write(tmp + '\n')
        return None

def TSV_file_into_JSON(src, dst, header):
        """
        https://stackoverflow.com/questions/48752209/python-script-tsv-conversion-to-json
        """
        data = []
        with open(src, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
                for row in reader:
                        if row[0] == 'sample_id':
                                print("\n")
                        else:
                                if row[0].strip()[0] == '#':  #
                                        continue
                                row = filter(None, row)
                                data.append(OrderedDict(zip(header, row)))

        with open(dst, 'w') as jsonfile:
                json.dump(data, jsonfile, indent=2)
        return None

def HTML_table(src, header, htmlEndStr):
        
        f = open(src)
        lines = f.readlines()
        f.close()
        
        with open('viral_results.html', 'w') as f: 
                f.write(header)
                for line in lines:
                        f.write(line)
                f.write(htmlEndStr)
        f.close()
        return None

def main():     
        #PATHWAY for creating initial TSV file
        #files_metaphlan2 = glob.glob('/home/amanj/Documents/Data/ttv_bonk/20190401/discovery/P*/reads/*metaphlan2.tsv')
        #files_fastviromeexplorer = glob.glob('/home/amanj/Documents/Data/ttv_bonk/20190401/discovery/P*/reads/*fastviromeexplorer_abundance.tsv')
        #files_kraken2 = glob.glob('/home/amanj/Documents/Data/ttv_bonk/20190401/discovery/P*/reads/*kraken2_report.txt')
        files_diamond_org.sort()
        files_diamond_tax_data.sort()
        files_kraken2.sort()
        CreateInitialTSVresults(files_diamond_org, files_diamond_tax_data, files_diamond_evalue, files_kraken2)
        #Output file TMPresults.txt
 
        ###################START_ADD_METADATA#############################
        if Bool_metadata == 1:
                #Previous viral results table from virom explorer fix
                #FstTsv = glob.glob('/home/amanj/Documents/190902/New_Script_all_in_one_parse_viral_data/TMPresultsFixedResults.txt')
                FstTsv = cur_dir+'/TMPresultsFixedResults.txt'
                #File pathway for 2nd tsv table with metadata
                #SndTsv = glob.glob('/home/amanj/Documents/190710/file_for_script/TTV_miseq_metadata.txt')
                #Files for human material aslo metadata
                #files_human = glob.glob('/home/amanj/Documents/190710/file_for_script/P*/humanrm_flagstat/*.flagstat')
                files_human.sort()
                #Files for total sequences also metadata
                #files_totalSeq = glob.glob('/home/amanj/Documents/190710/file_for_script/P*/raw_fastqc/P*/fastqc_data.txt')
                files_totalSeq.sort()
                #Function for creating new TSV containing metadata -> TMPresults_fixed_added_metadata.txt
                addMetaData(FstTsv, SndTsv, files_human, files_totalSeq)
        
        ###################CONVERT_TSV_file_INTO_JSON#############################
        if Bool_metadata == 1:
                src = 'TMPresults_fixed_added_metadata.txt'
        else:
                src = 'TMPresults.txt'
        dst = 'viral_results.json'
        header = ['sample_id', 'name','accession_number', 'kraken_nr_reads', 'percentage_of_identical_matches', 'expected_value', 'classification','divison', 'full_taxon_level:family_genus_species_subspecies']
        headerWithMetaData = [
                'sample_id', 'viral_name', 'kraken_nr_reads', 
                'fast_virome_explorer', 'metaphlan', 'classification', 'full_taxon_level:family_genus_species_subspecies',
                'name_of_pool', 'amplification', 'sequence_type', 'human_material', 'total_sequences'
        ]       
        if Bool_metadata == 1:
                TSV_file_into_JSON(src, dst, headerWithMetaData)
        else:
                TSV_file_into_JSON(src, dst, header)
                
        
        ###################CREATE_HTML_TABLE_FROM_JSON############################
        strHTMLheader = """
        
        <!doctype html>

        <html>
        <head>
        <meta charset="utf-8"/>
        <!-- 
        * Working with dynamic table
        * WebSite: http://www.dynatable.com/
        -->
        <!--    Bootstrap v2.3.2 -->
        <link rel="stylesheet" media="all" href="https://s3.amazonaws.com/dynatable-docs-assets/css/bootstrap-2.3.2.min.css" />
        <!-- Plugin styles -->
        <link rel="stylesheet" media="all" href="https://s3.amazonaws.com/dynatable-docs-assets/css/jquery.dynatable.css" />
        <!--  jQuery v3.0.0-beta1 -->
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.0.0-beta1/jquery.js"></script>
        <!-- JS Pluging -->
        <script type='text/javascript' src='https://s3.amazonaws.com/dynatable-docs-assets/js/jquery.dynatable.js'></script>
        <script type="text/javascript"> 
                $(document).ready( function() {
                $('#example').dynatable({
                dataset: {
                        records: JSON.parse($('#patients').text())
                }
                });
                });

        </script>               
        </head> 
                <body>
                        <div class = "container"  style="float:left;">
                                <table id="example" class="table table-striped table-bordered" cellspacing="0" width="100%">
                                        <thead>
                                                <tr>
                                                        <th>sample_id</th>
                                                        <th>name</th>
                                                        <th>accession_number</th>
                                                        <th>kraken_nr_reads</th>
                                                        <th>percentage_of_identical_matches</th>
                                                        <th>expected_value</th>
                                                        <th>classification</th>
                                                        <th>divison</th>
                                                        <th>full_taxon_level:family_genus_species_subspecies</th>
                                                </tr>
                                        </thead>

                                        <tfoot>
                                                <tr>
                                                        <th>sample_id</th>
                                                        <th>name</th>
                                                        <th>accession_number</th>
                                                        <th>kraken_nr_reads</th>
                                                        <th>percentage_of_identical_matches</th>
                                                        <th>expected_value</th>
                                                        <th>classification</th>
                                                        <th>divison</th>
                                                        <th>full_taxon_level:family_genus_species_subspecies</th>
                                                </tr>
                                        </tfoot>
                                         <tbody id="tbody">
                                        </tbody>
                                        <tbody>

                                        </tbody>
                                </table>
                         </div>
                <script id="patients">
                
        """
        
        strHTMLheaderWithMetaData = """
        
        <!doctype html>

        <html>
        <head>
        <meta charset="utf-8"/>
        <!-- 
        * Working with dynamic table
        * WebSite: http://www.dynatable.com/
        -->
        <!--    Bootstrap v2.3.2 -->
        <link rel="stylesheet" media="all" href="https://s3.amazonaws.com/dynatable-docs-assets/css/bootstrap-2.3.2.min.css" />
        <!-- Plugin styles -->
        <link rel="stylesheet" media="all" href="https://s3.amazonaws.com/dynatable-docs-assets/css/jquery.dynatable.css" />
        <!--  jQuery v3.0.0-beta1 -->
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.0.0-beta1/jquery.js"></script>
        <!-- JS Pluging -->
        <script type='text/javascript' src='https://s3.amazonaws.com/dynatable-docs-assets/js/jquery.dynatable.js'></script>
        <script type="text/javascript"> 
                $(document).ready( function() {
                $('#example').dynatable({
                dataset: {
                        records: JSON.parse($('#patients').text())
                }
                });
                });

        </script>               
        </head> 
                <body>
                        <div class = "container"  style="float:left;">
                                <table id="example" class="table table-striped table-bordered" cellspacing="0" width="100%">
                                        <thead>
                                                <tr>
                                                        <th>sample_id</th>
                                                        <th>viral_name</th>
                                                        <th>kraken_nr_reads</th>
                                                        <th>fast_virome_explorer</th>
                                                        <th>metaphlan</th>
                                                        <th>classification</th>
                                                        <th>full_taxon_level:family_genus_species_subspecies</th>
                                                        <th>name_of_pool</th>
                                                        <th>amplification</th>
                                                        <th>sequence_type</th>
                                                        <th>human_material</th>
                                                        <th>total_sequences</th>
                                                </tr>
                                        </thead>

                                        <tfoot>
                                                <tr>
                                                        <th>sample_id</th>
                                                        <th>viral_name</th>
                                                        <th>kraken_nr_reads</th>
                                                        <th>fast_virome_explorer</th>
                                                        <th>metaphlan</th>
                                                        <th>classification</th>
                                                        <th>full_taxon_level:family_genus_species_subspecies</th>
                                                        <th>name_of_pool</th>
                                                        <th>amplification</th>
                                                        <th>sequence_type</th>
                                                        <th>human_material</th>
                                                        <th>total_sequences</th>
                                                </tr>
                                        </tfoot>
                                         <tbody id="tbody">
                                        </tbody>
                                        <tbody>

                                        </tbody>
                                </table>
                         </div>
                <script id="patients">
                
        """
        strHTML_end = """
        
                        </script>
                </body>
         </html>
         """
        #dst = name of json file
        
        if Bool_metadata == 1:
                HTML_table(dst, strHTMLheaderWithMetaData, strHTML_end)
        else:
                HTML_table(dst, strHTMLheader, strHTML_end)
        
        ###### Moving files to output dir ########
        if Bool_filter == 0:
            # Create Dir
            if not os.path.exists("./viral_parser_contigs_output"):
                os.makedirs("./viral_parser_contigs_output")
        
            # Move HTML file to output dir
            os.rename("viral_results.html", "./viral_parser_contigs_output/viral_contigs_results.html")
        
            # Move Json file to output dir
            os.rename("viral_results.json", "./viral_parser_contigs_output/viral_contigs_results.json")
        
            # Move and rename final TSV file to output dir
            os.rename("TMPresults.txt", "./viral_parser_contigs_output/viral_contigs_results.tsv")
            os.rename("Unknown_sequences_from_diamond.txt", "./viral_parser_contigs_output/Unknown_sequences_from_diamond.tsv")
            if Bool_metadata == 1:
                os.rename("TMPresults_fixed_added_metadata.txt", "./viral_parser_contigs_output/TSV_with_metaData_viral_results.txt")
        else:
            if not os.path.exists("./parser_contigs_output"):
                os.makedirs("./parser_contigs_output")
        
            # Move HTML file to output dir
            os.rename("viral_results.html", "./parser_contigs_output/contigs_results.html")
        
            # Move Json file to output dir
            os.rename("viral_results.json", "./parser_contigs_output/contigs_results.json")
        
            # Move and rename final TSV file to output dir
            os.rename("TMPresults.txt", "./parser_contigs_output/contigs_results.tsv")
            os.rename("Unknown_sequences_from_diamond.txt", "./parser_contigs_output/Unknown_sequences_from_diamond.tsv")
            if Bool_metadata == 1:
                os.rename("TMPresults_fixed_added_metadata.txt", "./viral_parser_contigs_output/TSV_with_metaData_viral_results.txt")
        
        
        
        # Removing tmp files
        #os.remove("")
        files_tmp = glob.glob('./Combined_Results_*.txt')
        for tmpFile in files_tmp:
                os.remove(tmpFile)
                
        return None

main()
