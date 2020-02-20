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

fileTSVname = "P13409_picornavirus_T1D_ABIS"

#Magic-BLAST reads
magic_blast_selected_reads = glob.glob('/home/amanj/Documents/results_picornavirus/P13409_picorna_map_results/P*/reads/magicblast/*selected_reads_magicblast.txt')
magic_blast_sample_results_reads = glob.glob('/home/amanj/Documents/results_picornavirus/P13409_picorna_map_results/P*/reads/magicblast/*sample_results_magicblast_reads.txt')
#BWA reads
bwa_selected_reads = glob.glob('/home/amanj/Documents/results_picornavirus/P13409_picorna_map_results/P*/reads/BWA/*selected_reads_bwaBed.txt')
bwa_sample_results_reads = glob.glob('/home/amanj/Documents/results_picornavirus/P13409_picorna_map_results/P*/reads/BWA/*sample_results_bwaBed_reads.txt')
#Magic-BLAST contigs
magic_blast_selected_contig = glob.glob('/home/amanj/Documents/results_picornavirus/P13409_picorna_map_results/P*/contigs/magicblast/*selected_contig_magicblast.txt')
magic_blast_sample_results_contig = glob.glob('/home/amanj/Documents/results_picornavirus/P13409_picorna_map_results/P*/contigs/magicblast/*sample_results_magicblast_contig.txt')

def magic_blast_contig_parser(selected_contig, sample_result):
    container = []
    tmp = []
    f = open(selected_contig,"r")
    f1 = open(sample_result,"r")
    f_lines = f.readlines()
    f1_lines = f1.readlines()
    f.close()
    f1.close()
    
    for sel_con in f_lines:
        for samp_res in f1_lines:
            if sel_con.split('\t')[0] == samp_res.split('\t')[1]:
                tmp.append(samp_res.split('\t')[0]) 
                tmp.append("contig") 
                tmp.append(sel_con.split('\t')[1]) 
                tmp.append(sel_con.split('\t')[0]) 
                tmp.append(samp_res.split('\t')[2]) 
                tmp.append("-") 
                tmp.append("true") 
                container.append(tmp)
                #print(tmp)
                tmp = []
                break
                
    return container

def magic_blast_read_parser(selected_contig, sample_result):
    container = []
    tmp = []
    f = open(selected_contig,"r")
    f1 = open(sample_result,"r")
    f_lines = f.readlines()
    f1_lines = f1.readlines()
    f.close()
    f1.close()
    
    for sel_con in f_lines:
        for samp_res in f1_lines:
            if sel_con.split('\t')[0] == samp_res.split('\t')[1]:
                tmp.append(samp_res.split('\t')[0]) 
                tmp.append("read") 
                tmp.append(sel_con.split('\t')[1]) 
                tmp.append(sel_con.split('\t')[0]) 
                tmp.append(samp_res.split('\t')[2]) 
                tmp.append("false") 
                tmp.append("true") 
                container.append(tmp)
                #print(tmp)
                tmp = []
                break
                
    return container

def bwa_read_parser(selected_contig, sample_result):
    container = []
    tmp = []
    f = open(selected_contig,"r")
    f1 = open(sample_result,"r")
    f_lines = f.readlines()
    f1_lines = f1.readlines()
    f.close()
    f1.close()
    
    for sel_con in f_lines:
        for samp_res in f1_lines:
            if sel_con.split('\t')[0] == samp_res.split('\t')[1]:
                tmp.append(samp_res.split('\t')[0]) 
                tmp.append("read") 
                tmp.append(sel_con.split('\t')[1][:-1]) 
                tmp.append(sel_con.split('\t')[0]) 
                tmp.append(samp_res.split('\t')[2]) 
                tmp.append("true") 
                tmp.append("false") 
                container.append(tmp)
                #print(tmp)
                tmp = []
                break
                
    return container

def reads_join_containers(magicblast_con, bwa_con):
    list_index = []
    index = -1
    
    for i in range(len(magicblast_con)):
        for j in range(len(bwa_con)):
            index = i
            if bwa_con[j][2] == magicblast_con[i][2]:
                bwa_con[j][6] = "true"
                index = -1
                break
        if index != -1:
            list_index.append(index)
    for x in list_index:
        bwa_con.append(magicblast_con[x])
                
    return bwa_con

def HTML_table(src, header, htmlEndStr, filename):
        f = open(src)
        lines = f.readlines()
        f.close()
        with open(filename, 'w') as f: 
                f.write(header)
                for line in lines:
                        f.write(line)
                f.write(htmlEndStr)
        f.close()
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
        return;

def main():
    magic_blast_selected_reads.sort()
    magic_blast_sample_results_reads.sort()
    bwa_selected_reads.sort()
    bwa_sample_results_reads.sort()
    magic_blast_selected_contig.sort()
    magic_blast_sample_results_contig.sort()
    
    #print(magic_blast_selected_reads)
    #print(magic_blast_sample_results_reads)
    #print(bwa_selected_reads)
    #print(bwa_sample_results_reads)
    #print(magic_blast_selected_contig)
    #print(magic_blast_sample_results_contig)
    
    
    for m in range(len(magic_blast_selected_contig)):
        #contig container
        magicblast_contig_container = magic_blast_contig_parser(magic_blast_selected_contig[m],magic_blast_sample_results_contig[m])
        #reads container
        magicblast_read_container = magic_blast_read_parser(magic_blast_selected_reads[m],magic_blast_sample_results_reads[m])
        bwa_read_container = bwa_read_parser(bwa_selected_reads[m],bwa_sample_results_reads[m])
        magicblast_read_container.sort()
        bwa_read_container.sort()
        reads_container = reads_join_containers(magicblast_read_container,bwa_read_container)
        #reset data
        magicblast_read_container = []
        bwa_read_container = []
        tmp = ""
        fileTSVname1 = fileTSVname + ".tsv"
        fileTSVname2 = fileTSVname + ".json"
        fileTSVname3 = fileTSVname + ".html"
        f=open(fileTSVname1, "a+")
        for r in magicblast_contig_container:
            tmp = str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t"+str(r[3])+"\t"+str(r[4])+"\t"+str(r[5])+"\t"+str(r[6])+"\n"
            f.write(tmp)
        for r in reads_container:
            tmp = str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t"+str(r[3])+"\t"+str(r[4])+"\t"+str(r[5])+"\t"+str(r[6])+"\n"
            f.write(tmp)
    f.close() 
    Header = ['sample_id', 'type_of_data', 'seq_id', 'accession_number', 'ref_name', 'bwa', 'magicblast']
    TSV_file_into_JSON(fileTSVname1, fileTSVname2, Header)
    
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
                                                        <th>type_of_data</th>
                                                        <th>seq_id</th>
                                                        <th>accession_number</th>
                                                        <th>ref_name</th>
                                                        <th>bwa</th>
                                                        <th>magicblast</th>
                                                </tr>
                                        </thead>

                                        <tfoot>
                                                <tr>
                                                        <th>sample_id</th>
                                                        <th>type_of_data</th>
                                                        <th>seq_id</th>
                                                        <th>accession_number</th>
                                                        <th>ref_name</th>
                                                        <th>bwa</th>
                                                        <th>magicblast</th>
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
    HTML_table(fileTSVname2, strHTMLheader, strHTML_end, fileTSVname3)
    dirName = "./" + fileTSVname
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    # Move HTML file to output dir
    os.rename(fileTSVname3, dirName+"/"+fileTSVname3)
    # Move Json file to output dir
    os.rename(fileTSVname2, dirName+"/"+fileTSVname2)
    # Move Tsv fiel to output dir
    os.rename(fileTSVname1, dirName+"/"+fileTSVname1)

    return;
    
main()
    
    
    

