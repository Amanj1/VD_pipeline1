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

thisdict =	{
  "MH118073.1": "Coxsackievirus_A10_strain_10-4214-1",
  "KU866422.1": "Human poliovirus_1_strain_Mahoney_CDC",
  "KP247597.1": "Entrovirus_C_strain_strain_saukett_serotype3",
  "M12197.1": "Poliovirus_type2_completeGenome",
  "KF688606.1": "Rhinovirus C strain human/Australia/SG1/2008"
}

def HTML_table(src, header, htmlEndStr,name):
        f = open(src)
        lines = f.readlines()
        f.close()
        
        with open(name, 'w') as f: 
                f.write(header)
                for line in lines:
                        f.write(line)
                f.write(htmlEndStr)
        f.close()
        return None
        
files_blastn = glob.glob('/home/amanj/Documents/poli_and_coxs/poli_and_coxs_pipeline/poli_and_coxs_results/P*/blastn_reads/*.blastnTabular')
files_blastn.sort()
container = []
tmp = []
read = 0

for i in range(len(files_blastn)):
    Sample_ID = files_blastn[i].split('/').pop()
    ID = Sample_ID.split('_')
    #print(ID)
    Sample_ID = ID[0] + '_' + ID[1]
    
    f = open(files_blastn[i])
    lines = f.readlines()
    print
    for j in range(len(lines)):
        tmp.append(Sample_ID)
        tmp.append(lines[j].split('\t')[1])
        name =  thisdict[str(lines[j].split('\t')[1])]
        tmp.append(name)
        aln_len = lines[j].split('\t')[3]
        tmp.append(str(aln_len))
        tmp.append(lines[j].split('\t')[8])
        tmp.append(lines[j].split('\t')[9])
        seq_id = lines[j].split('\t')[0]
        if "read1" in str(ID[2].split('.')[0]):
            read = 1
        else:
            read = 2
        tmp.append(str(read))
        tmp.append(seq_id)
        tmp.append(lines[j].split('\t')[10])
        tmp.append(str(lines[j].split('\t')[2]))
        #print(tmp)
        container.append(tmp)  
        tmp = []
f.close()
name = "parser_blastn_output.tsv"
f=open(name, "a+")      
for c in container:
    #print(c)
    strline = ""
    for i in range(len(c)):
        if c[i] != c[len(c)-1]:
            strline = strline + c[i] + "\t"
        else:
            strline = strline + c[i] + '\n'
    f.write(strline)
f.close()

header = ['sample_id', 'ref_id', 'name', 'alignment_len', 'start_pos', 'end_pos', 'paired_end_read', 'seq_id', 'evalue', 'percentage_of_identical_matches' ]
data = []
with open("parser_blastn_output.tsv", 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    for row in reader:
        if row[0] == 'sample_id':
            print("\n")
        else:
            if row[0].strip()[0] == '#':  #
                continue
            row = filter(None, row)
            data.append(OrderedDict(zip(header, row)))

    with open("parser_blastn_output.json", 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
        
 
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
                                                        <th>ref_id</th>
                                                        <th>name</th>
                                                        <th>alignment_len</th>
                                                        <th>start_pos</th>
                                                        <th>end_pos</th>
                                                        <th>paired_end_read</th>
                                                        <th>seq_id</th>
                                                        <th>evalue</th>
                                                        <th>percentage_of_identical_matches</th>
                                                </tr>
                                        </thead>

                                        <tfoot>
                                                <tr>
                                                        <th>sample_id</th>
                                                        <th>ref_id</th>
                                                        <th>name</th>
                                                        <th>alignment_len</th>
                                                        <th>start_pos</th>
                                                        <th>end_pos</th>
                                                        <th>paired_end_read</th>
                                                        <th>seq_id</th>
                                                        <th>evalue</th>
                                                        <th>percentage_of_identical_matches</th>
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
HTML_table("parser_blastn_output.json", strHTMLheader, strHTML_end,"parser_blastn_output.html")
