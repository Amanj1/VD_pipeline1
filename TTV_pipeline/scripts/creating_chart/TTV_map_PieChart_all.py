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

contigs_magicblast = glob.glob('/proj/TTV_mapping/TTV_pipeline/TTV_map_results/all/contigs/magicblast/*.tsv')
reads_magicblast = glob.glob('/proj/TTV_mapping/TTV_pipeline/TTV_map_results/all/reads/magicblast/*.tsv')
reads_bwa = glob.glob('/proj/TTV_mapping/TTV_pipeline/TTV_map_results/all/reads/BWA/*.tsv')

File_output = "All_TTV_mapping_PieCharts.html"

def main():     
        contigs_magicblast.sort()
        reads_magicblast.sort()
        reads_bwa.sort()
        
        #print(contigs_magicblast)
        #print(reads_magicblast)
        #print(reads_bwa)


        ###################HTML_CODE_START############################
        strHTMLheader = """
        <!DOCTYPE html>
        <html lang="en-US">
        <body>

        <h1>Pie Chart BWA All Samples Reads</h1>
        <div id="piechart1"></div>
        <div id="piechart2"></div>
        <div id="piechart3"></div>

        <h1>Pie Chart Magicblast All Samples Reads</h1>
        <div id="piechart4"></div>
        <div id="piechart5"></div>
        <div id="piechart6"></div>

        <h1>Pie Chart Magicblast All Samples Contigs</h1>
        <div id="piechart7"></div>
        <div id="piechart8"></div>
        <div id="piechart9"></div>

        <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>

        <script type="text/javascript">
        // Load google charts
        google.charts.load('current', {'packages':['corechart']});
        google.charts.setOnLoadCallback(drawChart1);  //Reads tot
        google.charts.setOnLoadCallback(drawChart2);  //Reads all
        google.charts.setOnLoadCallback(drawChart3);  //Reads comb
        google.charts.setOnLoadCallback(drawChart4);  //Reads all
        google.charts.setOnLoadCallback(drawChart5);  //Reads comb
        google.charts.setOnLoadCallback(drawChart6);  //Contigs all
        google.charts.setOnLoadCallback(drawChart7);  //Contigs comb
        google.charts.setOnLoadCallback(drawChart8);  //Contigs comb
        google.charts.setOnLoadCallback(drawChart9);  //Contigs comb
        """
        
 ###################HTML_CODE_END############################
 
  ###################BWA############################
        strHTML_end = """
        </script>
        </body>
        </html>
         """
        func_bwa1 = """
        function drawChart1() {
        var data1 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_bwa_1_mid = """
        // Optional; add a title and set the width and height of the chart
        var options1 = {'title':'Total TTV hits from all samples (BWA reads)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart1 = new google.visualization.PieChart(document.getElementById('piechart1'));
        chart1.draw(data1, options1);
        }
        """
        func_end = """
        ]);
        """
        with open(reads_bwa[0]) as f:
            content = f.readlines()[1:]
        tmp = []
        tmp_html = []
        for c in content:
            tmp = c.split("\t")
            tmp_str = "['TTV',"+ str(tmp[6]) + "],\n"
            tmp_html.append(tmp_str)
            tot = int(tmp[7]) - int(tmp[6])
            tmp_str = "['Other Sequences',"+ str(tot) + "],"
            tmp_html.append(tmp_str)
            #print(tmp_html)
            break
        outF = open(File_output, "a")
        outF.write(strHTMLheader)
        outF.write(func_bwa1)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_bwa_1_mid)
        
        
        func_bwa2 = """
        function drawChart2() {
        var data2 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_bwa_2_mid = """
        // Optional; add a title and set the width and height of the chart
        var options2 = {'title':'TTV reference occurrences (BWA reads)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart2 = new google.visualization.PieChart(document.getElementById('piechart2'));
        chart2.draw(data2, options2);
        }
        """
        tmp = []
        tmp_html = []
        for c in content:
            tmp = c.split("\t")
            tmp_str = "['"+str(tmp[1])+"', "+ str(tmp[5]) + "],\n"
            tmp_html.append(tmp_str)
        f.close()
        outF.write(func_bwa2)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_bwa_2_mid)
        
        func_bwa3 = """
        function drawChart3() {
        var data3 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_bwa_3_mid = """
        // Optional; add a title and set the width and height of the chart
        var options3 = {'title':'TTV reference occurrences (BWA reads)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart3 = new google.visualization.PieChart(document.getElementById('piechart3'));
        chart3.draw(data3, options3);
        }
        """
        tot_nr_reads = 0
        with open(reads_bwa[1]) as f:
            content = f.readlines()[1:]
        tmp = []
        tmp_html = []
        run = 0
        for c in content:
            tmp = c.split("\t")
            if run == 0:
                sample_id = str(tmp[0])
                ttv_in_sample = str(tmp[7])
                tmp_str = "['"+sample_id+"', "+ ttv_in_sample + "],\n"
                tmp_html.append(tmp_str)
                run = 1
            else:
                if sample_id != str(tmp[0]):
                    sample_id = str(tmp[0])
                    ttv_in_sample = str(tmp[7])
                    tmp_str = "['"+sample_id+"', "+ ttv_in_sample + "],\n"
                    tmp_html.append(tmp_str)


        #tmp_str = "['"+"Other Sequences"+"', "+ str(tot_nr_reads) + "],\n"
        #tmp_html.append(tmp_str)
        f.close()
        outF.write(func_bwa3)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_bwa_3_mid)
        
   ###################MAGICBLAST READS############################       
        func_blast4 = """
        function drawChart4() {
        var data4 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_blast_4_mid = """
        // Optional; add a title and set the width and height of the chart
        var options4 = {'title':'Total TTV hits from all samples (Magicblast reads)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart4 = new google.visualization.PieChart(document.getElementById('piechart4'));
        chart4.draw(data4, options4);
        }
        """
        with open(reads_magicblast[0]) as f:
            content1 = f.readlines()[1:]
        tmp = []
        tmp_html = []
        content = []
        for c in content1:
            tmp = c.split("\t")
            tmp_str = "['TTV',"+ str(tmp[6]) + "],\n"
            tmp_html.append(tmp_str)
            tot = int(tmp[7]) - int(tmp[6])
            tmp_str = "['Other Sequences',"+ str(tot) + "],"
            tmp_html.append(tmp_str)
            #print(tmp_html)
            break
        outF.write(func_blast4)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_blast_4_mid)
        
        
        func_blast5 = """
        function drawChart5() {
        var data5 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_blast_5_mid = """
        // Optional; add a title and set the width and height of the chart
        var options5 = {'title':'TTV reference occurrences (Magicblast reads)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart5 = new google.visualization.PieChart(document.getElementById('piechart5'));
        chart5.draw(data5, options5);
        }
        """
        tmp = []
        tmp_html = []
        for c in content1:
            tmp = c.split("\t")
            tmp_str = "['"+str(tmp[1])+"', "+ str(tmp[5]) + "],\n"
            tmp_html.append(tmp_str)
        f.close()
        outF.write(func_blast5)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_blast_5_mid)
        
        func_blast6 = """
        function drawChart6() {
        var data6 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_blast_6_mid = """
        // Optional; add a title and set the width and height of the chart
        var options6 = {'title':'TTV reference occurrences (Magicblast reads)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart6 = new google.visualization.PieChart(document.getElementById('piechart6'));
        chart6.draw(data6, options6);
        }
        """
        tot_nr_reads = 0
        with open(reads_magicblast[1]) as f:
            content1 = f.readlines()[1:]
        tmp = []
        tmp_html = []
        run = 0
        for c in content1:
            tmp = c.split("\t")
            if run == 0:
                sample_id = str(tmp[0])
                ttv_in_sample = str(tmp[7])
                tmp_str = "['"+sample_id+"', "+ ttv_in_sample + "],\n"
                tmp_html.append(tmp_str)
                run = 1
            else:
                if sample_id != str(tmp[0]):
                    sample_id = str(tmp[0])
                    ttv_in_sample = str(tmp[7])
                    tmp_str = "['"+sample_id+"', "+ ttv_in_sample + "],\n"
                    tmp_html.append(tmp_str)


        #tmp_str = "['"+"Other Sequences"+"', "+ str(tot_nr_reads) + "],\n"
        #tmp_html.append(tmp_str)
        f.close()
        outF.write(func_blast6)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_blast_6_mid)   



###################MAGICBLAST CONTIGS############################       
        func_blast7 = """
        function drawChart7() {
        var data7 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_blast_7_mid = """
        // Optional; add a title and set the width and height of the chart
        var options7 = {'title':'Total TTV hits from all samples (Magicblast contigs)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart7 = new google.visualization.PieChart(document.getElementById('piechart7'));
        chart7.draw(data7, options7);
        }
        """
        with open(contigs_magicblast[0]) as f:
            content2 = f.readlines()[1:]
        tmp = []
        tmp_html = []
        content = []
        for c in content2:
            tmp = c.split("\t")
            tmp_str = "['TTV',"+ str(tmp[6]) + "],\n"
            tmp_html.append(tmp_str)
            tot = int(tmp[7]) - int(tmp[6])
            tmp_str = "['Other Sequences',"+ str(tot) + "],"
            tmp_html.append(tmp_str)
            #print(tmp_html)
            break
        outF.write(func_blast7)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_blast_7_mid)
        
        
        func_blast8 = """
        function drawChart8() {
        var data8 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_blast_8_mid = """
        // Optional; add a title and set the width and height of the chart
        var options8 = {'title':'TTV reference occurrences (Magicblast contigs)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart8 = new google.visualization.PieChart(document.getElementById('piechart8'));
        chart8.draw(data8, options8);
        }
        """
        tmp = []
        tmp_html = []
        for c in content2:
            tmp = c.split("\t")
            tmp_str = "['"+str(tmp[1])+"', "+ str(tmp[5]) + "],\n"
            tmp_html.append(tmp_str)
        f.close()
        outF.write(func_blast8)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_blast_8_mid)
        
        func_blast9 = """
        function drawChart9() {
        var data9 = google.visualization.arrayToDataTable([
        ['Sequence', 'Matches'],
        """
        func_blast_9_mid = """
        // Optional; add a title and set the width and height of the chart
        var options9 = {'title':'TTV reference occurrences (Magicblast contigs)', 'width':1050, 'height':900};

        // Display the chart inside the <div> element with id="piechart"
        var chart9 = new google.visualization.PieChart(document.getElementById('piechart9'));
        chart9.draw(data9, options9);
        }
        """
        tot_nr_reads = 0
        with open(contigs_magicblast[1]) as f:
            content2 = f.readlines()[1:]
        tmp = []
        tmp_html = []
        run = 0
        for c in content2:
            tmp = c.split("\t")
            if run == 0:
                sample_id = str(tmp[0])
                ttv_in_sample = str(tmp[7])
                tmp_str = "['"+sample_id+"', "+ ttv_in_sample + "],\n"
                tmp_html.append(tmp_str)
                run = 1
            else:
                if sample_id != str(tmp[0]):
                    sample_id = str(tmp[0])
                    ttv_in_sample = str(tmp[7])
                    tmp_str = "['"+sample_id+"', "+ ttv_in_sample + "],\n"
                    tmp_html.append(tmp_str)
        f.close()
        outF.write(func_blast9)
        outF.writelines(tmp_html)
        outF.write(func_end)
        outF.write(func_blast_9_mid)        

 ###################HTML_CODE_END############################
        outF.write(strHTML_end)
        outF.close()
        return None

main()
