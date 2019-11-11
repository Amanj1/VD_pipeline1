import sys
import os
import glob
import re
import string
import shutil


files_final_contigs = glob.glob('/home/amanj/Documents/TTV/discorvery_megahit_contigs/discovery_megahit/P*/megahit/1_assembly/final.contigs.fa')
files_filter_contigs = glob.glob('/home/amanj/Documents/TTV/discorvery_megahit_contigs/discovery_megahit/P*/megahit/2_filt_contigs/contigs_filt.fa')

files_final_contigs.sort()
files_filter_contigs.sort()

dstDir_final_contigs = '/home/amanj/Documents/TTV/discorvery_megahit_contigs/megahit_contigs/'
dstDir_filter_contigs = '/home/amanj/Documents/TTV/discorvery_megahit_contigs/megahit_filt_contigs/'

dst_file_final_contigs = '/home/amanj/Documents/TTV/discorvery_megahit_contigs/megahit_contigs/final.contigs.fa'
dst_file_filter_contigs = '/home/amanj/Documents/TTV/discorvery_megahit_contigs/megahit_filt_contigs/contigs_filt.fa'

for i in range(len(files_final_contigs)):
    final = files_final_contigs[i].split('/')[7]
    filter = files_filter_contigs[i].split('/')[7]

    dst_file_name = dstDir_final_contigs + final + "_final_contigs.fa"
    shutil.copy(files_final_contigs[i],dstDir_final_contigs)
    os.rename(dst_file_final_contigs, dst_file_name)
    
    dst_file_name = dstDir_filter_contigs + filter + "_filt_contigs.fa"
    shutil.copy(files_filter_contigs[i],dstDir_filter_contigs)
    os.rename(dst_file_filter_contigs, dst_file_name)