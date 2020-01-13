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
  "KP247597.1": "Entrovirus_C_strain_Saukett_serotype3",
  "M12197.1": "Poliovirus_type2_completeGenome",
  "KF688606.1": "Rhinovirus C strain human/Australia/SG1/2008",
  "GQ183022.1" :"Human parechovirus 1 isolate K129-93, complete genome",
  "k141_3825" : "k141_3825 flag=1 multi=1.0000 len=376",
  "k141_2276" : "k141_2276 flag=1 multi=14.6954 len=988",
  "k141_2271" : "k141_2271 flag=1 multi=2.0000 len=516",
  "k141_2208" : "k141_2208 flag=1 multi=5.0000 len=542",
  "k141_1466" : "k141_1466 flag=1 multi=4.0000 len=960",
  "k141_3265" : "k141_3265 flag=1 multi=6.8986 len=1522"
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
        
files_bed = glob.glob('/home/amanj/Documents/poli_and_coxs/poli_and_coxs_pipeline/poli_and_coxs_results/P*/BWA_reads/*_reads_bwaSamToBed.txt')
files_bed.sort()
container = []
tmp = []
read = 0
#Header = "Patient_ID\tName\tKraken_#n/home/amanj/Documents/190710"

#f=open(name, "a+")
#f.write(Header)
#f.close()
for i in range(len(files_bed)):
    Sample_ID = files_bed[i].split('/').pop()
    ID = Sample_ID.split('_')
    Sample_ID = ID[0] + '_' + ID[1]
    f = open(files_bed[i])
    lines = f.readlines()
    for j in range(len(lines)):
        tmp.append(Sample_ID)
        tmp.append(lines[j].split('\t')[0])
        name =  thisdict[str(lines[j].split('\t')[0])]
        tmp.append(name)
        seq_len = int(lines[j].split('\t')[2]) - int(lines[j].split('\t')[1])
        tmp.append(str(seq_len))
        tmp.append(lines[j].split('\t')[1])
        tmp.append(lines[j].split('\t')[2])
        seq_id = lines[j].split('\t')[3]
        if "/1" in seq_id:
            read = 1
            seq_id = seq_id.replace('/1','')
        else:
            read = 2
            seq_id = seq_id.replace('/2','')
        tmp.append(str(read))
        tmp.append(seq_id)
        tmp.append(lines[j].split('\t')[4])
        tmp.append(lines[j].split('\t')[5])
        container.append(tmp)  
        tmp = []
f.close()
name = "parser_BedTools_output.tsv"
f=open(name, "a+")      
for c in container:
    strline = ""
    for i in range(len(c)):
        if c[i] != c[len(c)-1]:
            strline = strline + c[i] + "\t"
        else:
            strline = strline + c[i]
    f.write(strline)
f.close()

header = ['sample_id', 'ref_id', 'name', 'alignment_len', 'start_pos', 'end_pos', 'paired_end_read', 'seq_id', 'bedtools_score', 'strand' ]
data = []
with open("parser_BedTools_output.tsv", 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    for row in reader:
        if row[0] == 'sample_id':
            print("\n")
        else:
            if row[0].strip()[0] == '#':  #
                continue
            row = filter(None, row)
            data.append(OrderedDict(zip(header, row)))

    with open("parser_BedTools_output.json", 'w') as jsonfile:
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
                                                        <th>bedtools_score</th>
                                                        <th>strand</th>
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
                                                        <th>bedtools_score</th>
                                                        <th>strand</th>
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
HTML_table("parser_BedTools_output.json", strHTMLheader, strHTML_end,"parser_BedTools_output.html")
