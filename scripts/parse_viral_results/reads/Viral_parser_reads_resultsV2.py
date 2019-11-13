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
Bool_metadata = 0

### Boolean for all results or only viral results. 1 for all or 0 for only viral.
Bool_filter = 0

#### PATHWAYS ####
###File pathway for creating table
files_metaphlan2 = glob.glob('/home/amanj/Documents/Results_from_pipeline/T1D_ABIS_Results/T1D_ABIS_discovery_readsOnly_results/discovery/P*/reads/*metaphlan2.tsv')
files_fastviromeexplorer = glob.glob('/home/amanj/Documents/Results_from_pipeline/T1D_ABIS_Results/T1D_ABIS_discovery_readsOnly_results/discovery/P*/reads/*fastviromeexplorer_abundance.tsv')
files_kraken2 = glob.glob('/home/amanj/Documents/Results_from_pipeline/T1D_ABIS_Results/T1D_ABIS_discovery_readsOnly_results/discovery/P*/reads/*kraken2_report.txt')

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
        tmp = str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t"+str(r[3])+"\t"+str(r[4])+"\t"+str(r[5])+"\t"+str(r[6])+"\n"
        f.write(tmp)
    f.close() 
    return None

def ireplace(old, new, text):
    index_l = text.lower().index(old.lower())
    return text[:index_l] + new + text[index_l + len(old):]

def parseMetaphlan2(fileInput):
    #Classification info
    """
    Kingdom[0]: k__, Phylum[1]: p__, Class[2]: c__, Order[3]: o__, Family[4]: f__, Genus[5]: g__, Species[6]: s__, Sub_species/Strain[7]: t__
    
    Since sequence-based profiling is relative and does not provide absolute cellular abundance measures, 
    clades are hierarchically summed. Each level will sum to 100%; that is, 
    the sum of all kindom-level clades is 100%, the sum of all genus-level clades (including unclassified) is also 100%, and so forth. 
    OTU equivalents can be extracted by using only the species-level s__ clades from this file (again, making sure to include clades unclassified at this level).
    
    Source: https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2
    
    Output from function:
    ['Classification', 'Name', 'Taxon Order', 'Abundance']
    ['Family', 'Flaviviridae', 'k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Flaviviridae', '100.0']  
    """
    #Header from table
    f = open(fileInput)
    line = f.readline().split()
    tmp =  []
    prevTmp = []
    prevTmpContent = []
    container = []
    classType = []
    
    #sample ID from header
    SampleID = line[1]

    #Reading file
    addClass = ['f__','g__','s__','t__']
    ClassNames = ['Family','Genus','Species','Sub_Species']
    line = f.readline().split()
    tmpline = line[0].split('|')
    prevTmp = tmpline
    prevTmpContent = line
    prevTmpContent.insert(0, prevTmp[len(prevTmp)-1])
    run = 1
    control = 1
    while run != 0:
        line = f.readline().split()
        if len(line) == 0:
            for i in range(len(addClass)):
                if addClass[i] in prevTmpContent[0]:
                    prevTmpContent[0] = prevTmpContent[0].replace(addClass[i],'')
                    prevTmpContent.insert(0, ClassNames[i])
            if "virus" in str(prevTmpContent[2].lower()) or "phage" in str(prevTmpContent[2].lower()) or "other" in str(prevTmpContent[2].lower()) or "unclassified" in str(prevTmpContent[2].lower()) or Bool_filter == 1:
                container.append(prevTmpContent)
            break
        tmpline = line[0].split('|')
        
        for n in prevTmp:
            if any(str(it) in str(n) for it in addClass):
                control = 0
                break
        if control != 1: 
            for i in range(len(addClass)):
                if addClass[i] in prevTmpContent[0]:
                    prevTmpContent[0] = prevTmpContent[0].replace(addClass[i],'')
                    prevTmpContent.insert(0, ClassNames[i])
            if "virus" in str(prevTmpContent[2].lower()) or "phage" in str(prevTmpContent[2].lower()) or "other" in str(prevTmpContent[2].lower()) or "unclassified" in str(prevTmpContent[2].lower()) or Bool_filter == 1:
                container.append(prevTmpContent)
        control = 1
        prevTmp = tmpline
        prevTmpContent = line
        prevTmpContent.insert(0, prevTmp[len(prevTmp)-1])
    container.insert(0, [SampleID])
    """
    for n in container:
        print(n)
        print('\n')
    """
    container = replaceOldDB("DictonaryDB.txt", container)
    return container

def fixFVElist(inputList):
    #print(inputList[0])
    newlist = [inputList[0]]
    inputList = inputList[1:len(inputList)]
    check = 1
    virusName = ''
    className = ''
    while check != 0:
        if ";" not in inputList[0]:
            virusName = virusName + inputList[0]
            inputList = inputList[1:len(inputList)]
            if ";" not in inputList[0]:
                virusName = virusName +"_"
        else:
            newlist.append(virusName)
            break
    while check != 0:
        if len(inputList) != 1:
            className = className + inputList[0]
            inputList = inputList[1:len(inputList)]
            if len(inputList) != 1:
                className = className +"_"
        else:
            newlist.append(className)
            newlist.append(inputList[0])
            break
    classStr = newlist[2].split(';')
    index = int(len(classStr));
    label = " "
    if index == 8:
        classStr = classStr[4:8]
        label = "Sub_Species"
    elif index == 7:
        classStr = classStr[4:7]
        label = "Species"
    elif index == 6:
        classStr = classStr[4:6]
        label = "Genus"
    elif index == 5:
        classStr = classStr[4:5]
        label = "Family"
    elif index == 4:
        classStr = classStr[4:4]
        label = "Family"
    else:
        classStr = ['delete']
        
    #Printing for debuggin index
    #print(index)
    #print(len(classStr))
    #print(classStr)
    
    #2 new colomns added index 0 gives taxon level and index 1 gives name of taxon level same as in methaphlan
    newlist.insert(0, classStr[len(classStr)-1])
    newlist.insert(0, label)
    return newlist

def parseFVE(fileInput):
    #Parsing FastViromeExplorer output data
    #Classification info
    """
    kingdom;    phylum; class;  order;  family; genus;  species
    
    Header:
    ['#VirusIdentifier', 'VirusName', 'kingdom;phylum;class;order;family;genus;species', 'EstimatedAbundance']
    
    unfixed
    ['NC_001710.1', 'GB', 'virus', 'C/Hepatitis', 'G', 'virus,', 'complete', 'genome', 'Unclassified;Unclassified;Unclassified;Unclassified;Flaviviridae;Pegivirus;Pegivirus', 'C', '66.0']
    
    fixed output from fixFVElist
    ['NC_001710.1', 'GB_virus_C/Hepatitis_G_virus,_complete_genome', 'Unclassified;Unclassified;Unclassified;Unclassified;Flaviviridae;Pegivirus;Pegivirus_C', '66.0']
    
    new fixed index 0 = taxon level and index 1 = name of taxon level
    ['Species', 'Human_mastadenovirus_C', 'AC_000008.1', 'Human_adenovirus_5,_complete_genome', 'Unclassified;Unclassified;Unclassified;Unclassified;Adenoviridae;Mastadenovirus;Human_mastadenovirus_C', '1248.75']  
    """
    #Skipping Header
    f = open(fileInput)
    line = f.readline().split()
    #Reading file FIXA loop för att läsa filen från FVE
    container = []
    run = 1
    while run != 0:
        line = f.readline().split()
        #print(line)
        if len(line) != 0:
            line = fixFVElist(line)
            container.append(line)
        else:
            break
    return container

def FVEaddSubSpecies(prevCont):
    container = []
    tmp = []
    tmpRow = []
    tmpRow2 = []
    subSpecies = []
    check = 0
    sumValues = 0.0
    for n in prevCont:
        tmpRow = [x for x in n]
        for m in prevCont:
            check = 0
            if m[2] != n[2] and m[1] == n[1]:
                tmpRow2 = [x for x in m]
                #n[len(n)-1] = float(n[len(n)-1]) + float(m[len(m)-1])
                sumValues = sumValues + float(m[len(m)-1])
                if subSpecies == []:
                    tmpRow[0] = "Sub_Species"
                    tmpRow[3] = ireplace(",", "", tmpRow[3])
                    tmpRow[1] = ireplace("_complete_genome", "", tmpRow[3])
                    subSpecies.append(tmpRow)
                else:
                    tmpRow2[0] = "Sub_Species"
                    tmpRow2[3] = ireplace(",", "", tmpRow2[3])
                    tmpRow2[1] = ireplace("_complete_genome", "", tmpRow2[3])
                    for t in subSpecies:
                        if t[2] == tmpRow2[2]:
                            check = 1
                    if check != 1:
                        subSpecies.append(tmpRow2)  
                check = 0
        if container == []:
            n[len(n)-1] = float(n[len(n)-1]) + sumValues
            container.append(n)
            sumValues = 0.0
        else:
            for c in container:
                if c[1] == n[1]:
                    check = 1
                    break
            if check != 1:
                n[len(n)-1] = float(n[len(n)-1]) + sumValues
                container.append(n)
                sumValues = 0.0
    for s in subSpecies:
        container.append(s)
    return container

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
            if "virus" not in str(line[5].lower()) or "phage" not in str(line[5].lower()) or "other" not in str(line[5].lower()) or "unclassified" not in str(line[5].lower()):
                check_level = 1
            if "virus" in str(line[5].lower()) or "phage" in str(line[5].lower()) or "other" in str(line[5].lower()) or "unclassified" in str(line[5].lower()):
                check_level = 0
        if check_level != 1 or Bool_filter == 1:    
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

def combineResults(Kraken, FVE, Meta, ID):
    """
    #Create header 
    # Patient_ID    Name    Kraken  FastViromeExplorer  MetaPhlan   Classification  Full_Taxon_level
    Header = "Patient_ID\tName\tKraken_#n/home/amanj/Documents/190710
    /home/amanj/Documents/190714rReads\tFastViromeExplorer\tMetaPhlan\tClassification\tFull_Taxon_level(Family->Genus->Species->Sub_species)\n"
    CreateFile(ID)
    name = "Combined_Results_"+ID+".txt"
    f=open(name, "a+")
    f.write(Header)
    f.close() 
    """
    #tmp = [ID,"Name","Kraken", "FVE", "Meta", "Class", "level"]
    tmp = [ID,"-","-", "-", "-", "-", "-"]
    container = []
    check = 1
    for k in Kraken:
        #Name
        tmp[1] = k[1]
        #Class
        tmp[5] = k[0]
        #Kraken
        tmp[2] = k[3]
        #Taxon level from Kraken
        tmp[6] = k[7]
        for f in FVE:
            if tmp[1] == f[1] and tmp[5] == f[0]:
                tmp[3] = f[len(f)-1]
                #tmp[6] = f[4].replace("Unclassified;", "")
                #tmp[6] = tmp[6].replace(";","->")
            for m in Meta:
                if tmp[1] == m[1] and tmp[5] == m[0]:
                    tmp[4] = m[len(m)-1]
        container.append(tmp)
        tmp = [ID,"-","-", "-", "-", "-", "-"]
        
    for f in FVE:
        tmp = [ID,"-","-", "-", "-", "-", "-"]
        check = 0
        for k in Kraken:
            if f[1] == k[1] and f[0] == k[0]:
                check = 1
        if check != 1:
            #Name
            tmp[1] = f[1]
            #Class
            tmp[5] = f[0]
            #FVE
            tmp[3] = f[len(f)-1]
            #level
            tmp[6] = f[4].replace("Unclassified;", "")
            tmp[6] = tmp[6].replace(";","->")
            for m in Meta:
                if tmp[1] == m[1] and tmp[5] == m[0]:
                    tmp[4] = m[len(m)-1]
            container.append(tmp)
        
    for m in Meta:
        tmp = [ID,"-","-", "-", "-", "-", "-"]
        check = 0
        for f in FVE:
            if f[1] == m[1] and f[0] == m[0]:
                check = 1
            for k in Kraken:
                if k[1] == m[1] and k[0] == m[0]:
                    check = 1
        if check != 1:
            #Name
            tmp[1] = m[1]
            #Class
            tmp[5] = m[0]
            #Meta
            tmp[4] = m[len(m)-1]
            tmp[6] = taxonLevelString(m[2].split("|"), 2)
            container.append(tmp)
    #for c in container:
    #   print(c)
    WriteToFile(container, ID)
    return None

def CreateInitialTSVresults(files_metaphlan2, files_fastviromeexplorer, files_kraken2):
    
    for i in range(len(files_metaphlan2)):
        fileName = files_metaphlan2[i]
        M_container = parseMetaphlan2(fileName)
        #Fetching Patient ID
        Patient_ID = M_container[0][0]
        M_container.pop(0)
        fileName = files_fastviromeexplorer[i]
        #fileName = "/home/amanj/Documents/190425/pcr_products_180912/discovery/P11463_1001_S1_L001/reads/P11463_1001_S1_L001_fastviromeexplorer_abundance.tsv"
        F_container = parseFVE(fileName)
        F_container = FVEaddSubSpecies(F_container)
        fileName = files_kraken2[i]
        K_container = parseKraken2(fileName)
        combineResults(K_container, F_container, M_container, Patient_ID)
    
    all_results = glob.glob('Combined_Results_*.txt')
    all_results.sort()
    
    #Create header 
    # Patient_ID    Name    Kraken  FastViromeExplorer  MetaPhlan   Classification  Full_Taxon_level
    Header = "sample_id\tviral_name\tkraken_nr_reads\tfast_virome_explorer\tmetaphlan\tclassification\tfull_taxon_level:family_genus_species_subspecies\n"
    f = open("TMPresults.txt", "w")
    f.write(Header)
    for result in all_results:
        f1 = open(result, "r")
        f.write(f1.read())
        f1.close()
    f.close()
    return None

def fixViromeExplorer(TSV_table):
    #Viral results table
    #TSV_table = glob.glob('/home/amanj/Documents/190902/New_Script_all_in_one_parse_viral_data/TMPresults.txt')
    listOfLists_FVE_results = []
    FVE_result_one_sample_id = []
    tmp = []
    #f = open(TSV_table[0])
    f = open(TSV_table)
    lines = f.readlines()
    header = lines[0]
    lines.pop(0)
    genusSum = 0
    familySum = 0
    for line in lines: 
        tmp = [line.split()[0], line.split()[5], line.split()[3]]
        #print(tmp)
        if len(FVE_result_one_sample_id) == 0 and 'family' in tmp[1].lower():
            FVE_result_one_sample_id.append(tmp)
        elif FVE_result_one_sample_id[0][0] == tmp[0] and 'family' not in tmp[1].lower():
            FVE_result_one_sample_id.append(tmp)
        else:
            #print(FVE_result_one_sample_id)
            #print('\n')
            listOfLists_FVE_results.append(FVE_result_one_sample_id)
            FVE_result_one_sample_id = []
            FVE_result_one_sample_id.append(tmp)
    new_FVE_results = []
    for i in listOfLists_FVE_results:
        tmp = i[::-1]
        for t in range(len(tmp)):
            if 'family' in tmp[t][1].lower():
                if  familySum != 0:
                    tmp[t][2] = str(familySum)
                    genusSum = 0
                    familySum = 0
            if 'species' in tmp[t][1].lower():
                if '-' == tmp[t][2]:
                    genusSum = genusSum
                else:
                    genusSum = genusSum + float(tmp[t][2])
                    familySum = familySum + genusSum
            if 'genus' in tmp[t][1].lower():
                if '-' == tmp[t][2]:
                    if genusSum != 0:
                        tmp[t][2] = str(genusSum)
                        genusSum = 0
                else:
                    genusSum = genusSum + float(tmp[t][2])
                    familySum = familySum + genusSum
                    if genusSum != 0:
                        tmp[t][2] = str(genusSum)
                        genusSum = 0
        new_FVE_results.append(tmp[::-1])
        
    #for i in new_FVE_results:
        #print(i)
    
    with open('TMPresultsFixedResults.txt', 'w') as f: 
        f.write(header.lower())
        l = 0 #counter for lines
        for i in new_FVE_results:
            for j in i:
                tmp = lines[l].split()
                #print(tmp)
                tmp[3] = j[2]
                #print(tmp)
                tmp = '\t'.join(tmp).rstrip('\n')
                f.write(tmp + '\n')
                l = l + 1
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
    files_metaphlan2.sort()
    files_fastviromeexplorer.sort()
    files_kraken2.sort()
    CreateInitialTSVresults(files_metaphlan2, files_fastviromeexplorer, files_kraken2)
    #Output file TMPresults.txt
    
    ###################VIROME_EXPLORER_FIX#############################
    #Viral results table. Pathway for initial TSV file
    #TSV_table = glob.glob('/home/amanj/Documents/190902/New_Script_all_in_one_parse_viral_data/TMPresults.txt')
    TSV_table = cur_dir+'/TMPresults.txt'
    fixViromeExplorer(TSV_table)
    #Output file TMPresultsFixedResults.txt
    
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
        src = 'TMPresultsFixedResults.txt'
    dst = 'viral_results.json'
    header = ['sample_id', 'viral_name', 'kraken_nr_reads', 'fast_virome_explorer', 'metaphlan', 'classification', 'full_taxon_level:family_genus_species_subspecies']
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
                            <th>viral_name</th>
                            <th>kraken_nr_reads</th>
                            <th>fast_virome_explorer</th>
                            <th>metaphlan</th>
                            <th>classification</th>
                            <th>full_taxon_level:family_genus_species_subspecies</th>
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
    
    if Bool_filter != 1:
    # Create Dir
        if not os.path.exists("./viral_parser_reads_output"):
            os.makedirs("./viral_parser_reads_output")
    
        # Move HTML file to output dir
        os.rename("viral_results.html", "./viral_parser_reads_output/viral_reads_results.html")
    
        # Move Json file to output dir
        os.rename("viral_results.json", "./viral_parser_reads_output/viral_reads_results.json")
    
        # Move and rename final TSV file to output dir
        os.rename("TMPresultsFixedResults.txt", "./viral_parser_reads_output/viral_reads_results.tsv")
        if Bool_metadata == 1:
            os.rename("TMPresults_fixed_added_metadata.txt", "./viral_parser_reads_output/viral_reads_results_with_metadata.tsv")
    else:
        if not os.path.exists("./parser_reads_output"):
            os.makedirs("./parser_reads_output")
    
        # Move HTML file to output dir
        os.rename("viral_results.html", "./parser_reads_output/reads_results.html")
    
        # Move Json file to output dir
        os.rename("viral_results.json", "./parser_reads_output/reads_results.json")
    
        # Move and rename final TSV file to output dir
        os.rename("TMPresultsFixedResults.txt", "./parser_reads_output/reads_results.tsv")
        if Bool_metadata == 1:
            os.rename("TMPresults_fixed_added_metadata.txt", "./parser_reads_output/reads_results_with_metadata.tsv")
    
    # Removing tmp files
    os.remove("TMPresults.txt")
    files_tmp = glob.glob('./Combined_Results_*.txt')
    for tmpFile in files_tmp:
        os.remove(tmpFile)
    
    return None

main()
