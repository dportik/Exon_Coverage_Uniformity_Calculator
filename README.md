# Exon_Coverage_Uniformity_Calculator

The uniformity of coverage across the length of exons can be examined using both longer (> 200 bp) 
and shorter (61–100 bp) exons.  Exons matching these criteria are filtered out from bed files 
containing exon coordinates independently for each sample.  For longer exons, five bins of 
10 bp increments are created for both the 5’ and 3’ ends, resulting in the generation of ten 
additional bed files per sample.  For shorter exons, three bins of 10 bp increments are created 
for both the 5’ and 3’ ends, resulting in the generation of six additional bed files per sample.  
Each bed file is used to calculate the per base pair coverage for a specific end bin with the 
samtools depth function using the appropriate sample bam file.  These per base pair coverage values 
are averaged within exons, and all averages of exons for a particular bin are subsequently 
combined across 50 randomly chosen samples (can be changed). The values across bins can be
visualized using the boxplot function in R to assess the median, upper and lower quartiles, 
and range of coverage estimates.  

Bed file generation scheme:
if exon diff between start/stop is >= 200 bp:

	1. start to (start+10) -> SAMPLE_long_forward0_10bp.bed
	2. (start+10) to (start+20) -> SAMPLE_long_forward10_20bp.bed
	3. (start+20) to (start+30) -> SAMPLE_long_forward20_30bp.bed
	4. (start+30) to (start+40) -> SAMPLE_long_forward30_40bp.bed
	5. (start+40) to (start+50) -> SAMPLE_long_forward40_50bp.bed

	1. end to (end-10)  -> SAMPLE_long_reverse0_10bp.bed
	2. (end-10) to (end-20) -> SAMPLE_long_reverse10_20bp.bed
	3. (end-20) to (end-30) -> SAMPLE_long_reverse20_30bp.bed
	4. (end-30) to (end-40) -> SAMPLE_long_reverse30_40bp.bed
	5. (end-40) to (end-50) -> SAMPLE_long_reverse40_50bp.bed
 

For short exons 61-100 bp:
if exon diff between start/stop is >60 and <=100:
	1. start to (start+10) -> short_forward0_10bp.bed
	2. (start+10) to (start+20) -> short_forward10_20bp.bed
	3. (start+20) to (start+30) -> short_forward20_30bp.bed

	1. end to (end-10)  -> short_reverse0_10bp.bed
	2. (end-10) to (end-20) -> short_reverse10_20bp.bed
	3. (end-20) to (end-30) -> short_reverse20_30bp.bed


------------------------------------------------------------------------------------------
STEP 1 - 1_exon_ends_bed_generator_v1.py


Usage: python 1_exon_ends_bed_generator_v1.py [directory with bed files] [path to an output directory]

This script will iterate across all the 'Sample_indexXX_coding.bed' files (initially found in the 
'assemblies/In_target' folder resulting from script 5 [5-FindingTargetsPhylo] from the 
https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow), which should be in a 
new directory of your choosing.

There are many outputs that will be created and moved to the output directory for use in the next script.
These are all bed files that contain information about the exon end bins.

When finished, run the second script: '2_coverage_uniformity_v1.py'

------------------------------------------------------------------------------------------
STEP 2 - 2_coverage_uniformity_v1.py

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
