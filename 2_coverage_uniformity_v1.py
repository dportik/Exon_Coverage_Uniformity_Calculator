import sys
import os
import subprocess as sp
import shutil
import numpy as np
from random import shuffle
import time
'''
Usage: python 2_coverage_uniformity_v1.py [directory with bam files] [path to new bed files output directory]


The 'Sample_indexXX__sorted.bam' files are located in the 'intarget_assemblies' folder resulting 
from script 6 (6-TransExonCapPhylo 'contig' option) of the 
https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow.
You should consider moving these to a new directory before beginning here.

The script currently uses 50 random samples to calculate the end bin coverages across 
all exons passing short or long length requirements.  This takes a while.  You can edit 
the script near the ******* below to include all samples, or decrease this number.  

The 'All_out.txt' output has the average coverage for each exon edge bin across all samples
included.  This can be used to generate box plots in R using the bin name as the factor.

##############
DEPENDENCIES:
numpy - Numerical Python
samtools - needs to be in path to call from command line
##############
------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
October 2015
------------------------
'''

bam_dir = sys.argv[1]
os.chdir(bam_dir)

prefix_set = set()

for fl in os.listdir('.'):
    if fl.endswith('_sorted.bam'):
        flnames = fl.split('_')
        prefix = flnames[0]+'_'+flnames[1]
        prefix_set.add(prefix)
        
prefix_list = list(prefix_set)
#shuffle list to randomize sampling to draw from
shuffle(prefix_list)

bed_dir = sys.argv[2]
os.chdir(bed_dir)

cov_list = []

x = int(1)

#************take 50 random samples from the list, which is out of order
#************to take less samples, change 51 below to lower number (ex. 21 = 20 samples)
for pref in prefix_list[:51]:
#************OR USE THE COMMAND BELOW TO INCLUDE ALL SAMPLES (this could take a long time...)
#for pref in prefix_list:
    print '\n'
    for fl in os.listdir('.'):
        flsplit = fl.split('.')
        if flsplit[0] == pref:
            out_name = flsplit[0]+'.'+flsplit[1]+'.coverage.txt'
            cov_list.append(flsplit[1])

            bam_name = bam_dir+'/'+pref+'_sorted.bam'
            #samtools depth -b  bed  sorted.bam > per-site-coverage.txt
            print "Sample {0} - {1}: Now calculating coverage for {2}".format(x,pref,fl)
            cov_string = "samtools depth -b {0} {1} > {2}".format(fl, bam_name, out_name)
            proc_cov = sp.call(cov_string, shell=True)
    x+=1
    print time.asctime( time.localtime(time.time()) )

    
def quickstats(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = str(np.around(x_avg, decimals = 1))
    return x_avg

out_file = 'All_out.txt'
fh_out = open(out_file, 'a')

long_forward0_10bp = []
long_forward10_20bp = []
long_forward20_30bp = []
long_forward30_40bp = []
long_forward40_50bp = []

long_reverse0_10bp = []
long_reverse10_20bp = []
long_reverse20_30bp = []
long_reverse30_40bp = []
long_reverse40_50bp = []
        
short_forward0_10bp = []
short_forward10_20bp = []
short_forward20_30bp = []

short_reverse0_10bp = []
short_reverse10_20bp = []
short_reverse20_30bp = []

            
for fl in os.listdir('.'):
    if fl.endswith('long_forward0_10bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_forward0_10bp.append(contig_avg)
        fh_temp.close()
        
    elif fl.endswith('long_forward10_20bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_forward10_20bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_forward20_30bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_forward20_30bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_forward30_40bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_forward30_40bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_forward40_50bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_forward40_50bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_reverse0_10bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_reverse0_10bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_reverse10_20bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_reverse10_20bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_reverse20_30bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_reverse20_30bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_reverse30_40bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_reverse30_40bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('long_reverse40_50bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            long_reverse40_50bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('short_forward0_10bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            short_forward0_10bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('short_forward10_20bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            short_forward10_20bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('short_forward20_30bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            short_forward20_30bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('short_reverse0_10bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            short_reverse0_10bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('short_reverse10_20bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            short_reverse10_20bp.append(contig_avg)
        fh_temp.close()

    elif fl.endswith('short_reverse20_30bp.coverage.txt'):
        print "Working to calculate averages for {}...".format(fl)
        temp_set = set()
        fh_temp = open(fl,'r')
        lines = fh_temp.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')
            temp_set.add(line[0])
        temp_list = list(temp_set)
        for contig in temp_list:
            concovlist = []
            for line in lines:
                line = line.strip()
                line = line.split('\t')
                if line[0] == contig:
                    concovlist.append(line[2])
            contig_avg = quickstats(concovlist)
            short_reverse20_30bp.append(contig_avg)
        fh_temp.close()

print "Writing averages set 1 (of 16) to output..."
for item in long_forward0_10bp:
    fh_out.write('long_forward0_10bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 2 (of 16) to output..."
for item in long_forward10_20bp:
    fh_out.write('long_forward10_20bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 3 (of 16) to output..."
for item in long_forward20_30bp:
    fh_out.write('long_forward20_30bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 4 (of 16) to output..."
for item in long_forward30_40bp:
    fh_out.write('long_forward30_40bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 5 (of 16) to output..."
for item in long_forward40_50bp:
    fh_out.write('long_forward40_50bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 6 (of 16) to output..."
for item in long_reverse0_10bp:
    fh_out.write('long_reverse0_10bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 7 (of 16) to output..."
for item in long_reverse10_20bp:
    fh_out.write('long_reverse10_20bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 8 (of 16) to output..."
for item in long_reverse20_30bp:
    fh_out.write('long_reverse20_30bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 9 (of 16) to output..."
for item in long_reverse30_40bp:
    fh_out.write('long_reverse30_40bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 10 (of 16) to output..."
for item in long_reverse40_50bp:
    fh_out.write('long_reverse40_50bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 11 (of 16) to output..."
for item in short_forward0_10bp:
    fh_out.write('short_forward0_10bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 12 (of 16) to output..."
for item in short_forward10_20bp:
    fh_out.write('short_forward10_20bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 13 (of 16) to output..."
for item in short_forward20_30bp:
    fh_out.write('short_forward20_30bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 14 (of 16) to output..."
for item in short_reverse0_10bp:
    fh_out.write('short_reverse0_10bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 15 (of 16) to output..."
for item in short_reverse10_20bp:
    fh_out.write('short_reverse10_20bp'+'\t'+'{}'.format(item)+'\n')
print "Writing averages set 16 (of 16) to output..."
for item in short_reverse20_30bp:
    fh_out.write('short_reverse20_30bp'+'\t'+'{}'.format(item)+'\n')

fh_out.close()

print '\n', "Process complete! That took a while...", '\n'
