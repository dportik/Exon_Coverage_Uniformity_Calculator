import sys
import os
import subprocess as sp
import shutil
import numpy as np
'''
Usage: python 1_exon_ends_bed_generator_v1.py [directory with bed files] [path to an output directory]

This script will iterate across all the 'Sample_indexXX_coding.bed' files (initially found in the 
'assemblies/In_target' folder resulting from script 5 [5-FindingTargetsPhylo] from the 
https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow), which should be in a 
new directory of your choosing.

There are many outputs that will be created and moved to the output directory for use in the next script.
These are all bed files that contain information about the exon end bins.

When finished, run the second script: '2_coverage_uniformity_v1.py'


------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
October 2015
------------------------
'''

bed_dir = sys.argv[1]
os.chdir(bed_dir)

output_dir = sys.argv[2]
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

prefix_list = []

#look through s5 target folder for coding.bed files
for fl in os.listdir('.'):
    if fl.endswith('_coding.bed'):
        temp_fh = open(fl, 'r')
        print "Beginning short and long fragment outputs for {}...".format(fl)
        flnames = fl.split('_')
        prefix = flnames[0]+'_'+flnames[1]

        prefix_list.append(prefix)
        
        fh_1 = open('{}.long_forward0_10bp.bed'.format(prefix), 'a')
        fh_2 = open('{}.long_forward10_20bp.bed'.format(prefix), 'a')
        fh_3 = open('{}.long_forward20_30bp.bed'.format(prefix), 'a')
        fh_4 = open('{}.long_forward30_40bp.bed'.format(prefix), 'a')
        fh_5 = open('{}.long_forward40_50bp.bed'.format(prefix), 'a')
        fh_6 = open('{}.long_reverse0_10bp.bed'.format(prefix), 'a')
        fh_7 = open('{}.long_reverse10_20bp.bed'.format(prefix), 'a')
        fh_8 = open('{}.long_reverse20_30bp.bed'.format(prefix), 'a')
        fh_9 = open('{}.long_reverse30_40bp.bed'.format(prefix), 'a')
        fh_10 = open('{}.long_reverse40_50bp.bed'.format(prefix), 'a')
        
        fh_11 = open('{}.short_forward0_10bp.bed'.format(prefix), 'a')
        fh_12 = open('{}.short_forward10_20bp.bed'.format(prefix), 'a')
        fh_13 = open('{}.short_forward20_30bp.bed'.format(prefix), 'a')
        fh_14 = open('{}.short_reverse0_10bp.bed'.format(prefix), 'a')
        fh_15 = open('{}.short_reverse10_20bp.bed'.format(prefix), 'a')
        fh_16 = open('{}.short_reverse20_30bp.bed'.format(prefix), 'a')
        
        for line in temp_fh:
            line = line.strip()
            line = line.split('\t')
            sample_name = line[0]
            start = int(line[1])
            end = int(line[2])
            length = end - start
            start10 = start+10
            start20 = start+20
            start30 = start+30
            start40 = start+40
            start50 = start+50
            end10 = end-10
            end20 = end-20
            end30 = end-30
            end40 = end-40
            end50 = end-50
            #create outputs of 10 bins (5 forward, 5 reverse) for every exon over 200 bp
            if length >= int(200):
            	print '\t', "exon is > 200..."
                fh_1.write(sample_name+'\t'+'{}'.format(start)+'\t'+'{}'.format(start10)+'\n')
                fh_2.write(sample_name+'\t'+'{}'.format(start10)+'\t'+'{}'.format(start20)+'\n')
                fh_3.write(sample_name+'\t'+'{}'.format(start20)+'\t'+'{}'.format(start30)+'\n')
                fh_4.write(sample_name+'\t'+'{}'.format(start30)+'\t'+'{}'.format(start40)+'\n')
                fh_5.write(sample_name+'\t'+'{}'.format(start40)+'\t'+'{}'.format(start50)+'\n')
                fh_6.write(sample_name+'\t'+'{}'.format(end10)+'\t'+'{}'.format(end)+'\n')
                fh_7.write(sample_name+'\t'+'{}'.format(end20)+'\t'+'{}'.format(end10)+'\n')
                fh_8.write(sample_name+'\t'+'{}'.format(end30)+'\t'+'{}'.format(end20)+'\n')
                fh_9.write(sample_name+'\t'+'{}'.format(end40)+'\t'+'{}'.format(end30)+'\n')
                fh_10.write(sample_name+'\t'+'{}'.format(end50)+'\t'+'{}'.format(end40)+'\n')

            #create outputs of 6 bins (3 forward, 4 reverse) for every exon between 60 and 200 bp
            if length > int(60) and length <= int(100):
            	print '\t', "exon is > 60 and < 100 bp...."
                fh_11.write(sample_name+'\t'+'{}'.format(start)+'\t'+'{}'.format(start10)+'\n')
                fh_12.write(sample_name+'\t'+'{}'.format(start10)+'\t'+'{}'.format(start20)+'\n')
                fh_13.write(sample_name+'\t'+'{}'.format(start20)+'\t'+'{}'.format(start30)+'\n')
                fh_14.write(sample_name+'\t'+'{}'.format(end10)+'\t'+'{}'.format(end)+'\n')
                fh_15.write(sample_name+'\t'+'{}'.format(end20)+'\t'+'{}'.format(end10)+'\n')
                fh_16.write(sample_name+'\t'+'{}'.format(end30)+'\t'+'{}'.format(end20)+'\n')
                
        fh_1.close()
        fh_2.close()
        fh_3.close()
        fh_4.close()
        fh_5.close()
        fh_6.close()
        fh_7.close()
        fh_8.close()
        fh_9.close()
        fh_10.close()
        fh_11.close()
        fh_12.close()
        fh_13.close()
        fh_14.close()
        fh_15.close()
        fh_16.close()
        
        temp_fh.close()
        
        print "Finished writing outputs for file {}.".format(fl), '\n'

print "Now moving all files to output directory...", '\n'
for fl in os.listdir('.'):
    if fl.endswith('_10bp.bed'):
        shutil.move(fl, output_dir)
    elif fl.endswith('_20bp.bed'):
        shutil.move(fl, output_dir)
    elif fl.endswith('_30bp.bed'):
        shutil.move(fl, output_dir)
    elif fl.endswith('_40bp.bed'):
        shutil.move(fl, output_dir)
    elif fl.endswith('_50bp.bed'):
        shutil.move(fl, output_dir)
print "All finished!", '\n'

        
#now all files are generated, need to 
#samtools depth -b  bed  sorted.bam > per-site-coverage.txt for each of the 10 bed files

#JMPD002_index4_sorted.bam
#are in: /Volumes/Portik_Storage/S6A_Contig/S6_L2_1/intarget_assemblies
#All bed files moved to output folder


