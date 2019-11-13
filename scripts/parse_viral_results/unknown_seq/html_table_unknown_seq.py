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

def HTML_table(src, header, htmlEndStr, name):
    
    f = open(src)
    lines = f.readlines()
    f.close()
    
    with open(name+'.html', 'w') as f: 
        f.write(header)
        for line in lines:
            f.write(line)
        f.write(htmlEndStr)
    f.close()
    return None

def main(): 

    
    ###################CONVERT_TSV_file_INTO_JSON#############################
    src = sys.argv[1]
    print(src)
    dst = src.split('.')[0]+'.json'
    print(dst)
    header = [ 'sample_id',	'accession_number',	'sequence_id',	'expected_value',	'percentage_of_identical_matches'	,'title',	'alignment_len', 'sequence_len']  
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
                            <th>accession_number</th>
                            <th>sequence_id</th>
                            <th>expected_value</th>
                            <th>percentage_of_identical_matches</th>
                            <th>title</th>
                            <th>alignment_len</th>
                            <th>sequence_len</th>
                        </tr>
                    </thead>

                    <tfoot>
                        <tr>
                            <th>sample_id</th>
                            <th>accession_number</th>
                            <th>sequence_id</th>
                            <th>expected_value</th>
                            <th>percentage_of_identical_matches</th>
                            <th>title</th>
                            <th>alignment_len</th>
                            <th>sequence_len</th>
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

    HTML_table(dst, strHTMLheader, strHTML_end, src.split('.')[0])
    

    
    return None

main()
