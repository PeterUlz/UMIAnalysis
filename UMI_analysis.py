#! /usr/bin/python

# Analyze Experiments involving molecular barcodes

# Assumptions:
#  1) FastQ File is 4 lines per read (sequence does not wrap into two lines)
#  2) Reads from same read group are VERY similar
#  3) InDel (errors) are very rare (InDels would lead to strand lagging in consensus building step)

# version 0.2: Allow adaptor trimming (usnig cutadapt)
# version 0.1: Initial set up
version="0.2"
import argparse
import gzip
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from UMI_Reads import UMIRead
from UMI_Reads import UMIReadGroup
###################################################################################

script_dir = os.path.dirname(os.path.realpath(__file__))
###################################################################################
def saveSettings(args,filename):
   SETTINGS = open(filename,"w")
   SETTINGS.write(time.strftime("%d/%m/%Y:  %H:%M:%S  : Starting Analysis"))
   SETTINGS.write("  Version: "+version+"\n")
   SETTINGS.write("  Fastq file 1: "+args.fastq_file1+"\n")
   SETTINGS.write("  Fastq file 2: "+args.fastq_file2+"\n")
   SETTINGS.write("  Length UMI Read 1: "+str(args.umi1)+"\n")
   SETTINGS.write("  Length UMI Read 2: "+str(args.umi2)+"\n")
   SETTINGS.write("  Length Offset Read 1: "+str(args.off1)+"\n")
   SETTINGS.write("  Length Offset Read 2: "+str(args.off2)+"\n")
   SETTINGS.write("  Member Threshold: "+str(args.member_threshold)+"\n")
   SETTINGS.write("  Consensus Threshold: "+str(args.consensus)+"\n")
   SETTINGS.write("  Threads: "+str(args.threads)+"\n")
   SETTINGS.write("  Sample name file: "+args.name+"\n")
   SETTINGS.write("  Output Directory: "+args.outdir+"\n")
   SETTINGS.close()

###################################################################################
# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription end')
parser.add_argument('-fq1','--fastq1', dest='fastq_file1', 
                   help='Fastq file R1',required=True)
parser.add_argument('-fq2','--fastq2', dest='fastq_file2', 
                   help='Fastq file R2',required=True)
parser.add_argument('-umi1','--umi-read1', dest='umi1', 
                   help='Length of UMI in Read 1 (starting from 5\' end)',required=True,type=int,default=6)
parser.add_argument('-umi2','--umi-read2', dest='umi2', 
                   help='Length of UMI in Read 2 (starting from 5\' end)',type=int,default=0)
parser.add_argument('-off1','--offset1', dest='off1', 
                   help='Length of offset after UMI 1 to trim away from read',type=int,default = 0)
parser.add_argument('-off2','--offset2', dest='off2', 
                   help='Length of offset after UMI 2 to trim away from read',type=int,default = 0)
parser.add_argument('-members','--members-threshold', dest='member_threshold', 
                   help='Minimum amount of UMI Members of Read Group in order to output this sequence [default 1]',type=int,default = 1)
parser.add_argument('-cons','--consensus-threshold', dest='consensus', 
                   help='Minimum consensus of bases in read group for assignment [default 0.5]',type=float,default = 0.5)
parser.add_argument('-s','--sample-name', dest='name',
                   help='Sample name to be used subsequently',required=True)
parser.add_argument('-o','--out-dir', dest='outdir',
                   help='Output Directory [default .]',default=".")
parser.add_argument('-t','--threads', dest='threads',
                   help='No. threads for alignment [default: 1]',type=int,default=1)
parser.add_argument('-mutect','--mutect-regions', dest='target_regions',
                   help='BED file containing regions for SNP-calling with MuTect [default: ""]',default="")
parser.add_argument('-adaptor1','--adaptor1', dest='adaptor1',
                   help='Adaptor sequence to trim in collapsed reads [default: ""]',default="")
parser.add_argument('-adaptor2','--adaptor2', dest='adaptor2',
                   help='Adaptor sequence to trim in collapsed reads [default: ""]',default="")
parser.add_argument('-cons_type','--consensus-type', dest='cons',
                   help='Generate exact match consensus or Multiple Alignment (msa) [default: exact]',default="exact",choices=["exact","msa"])
parser.add_argument('-skip_aln','--skip-aligment', dest='skip_aln',
                   help='Stop analysis after Collapsed FastQ Generation',action="store_true")
parser.add_argument('-skip_snp','--skip-snp-calling', dest='skip_snp',
                   help='Stop analysis after Alignment',action="store_true")

args = parser.parse_args()
print time.strftime("%d/%m/%Y:  %H:%M:%S  : Starting Analysis")
print "  Fastq file 1: ",args.fastq_file1
print "  Fastq file 2: ",args.fastq_file2
print "  Length UMI Read 1: ",args.umi1
print "  Length UMI Read 2: ",args.umi2
print "  Sample name file: ",args.name
print "  Output Directory: ",args.outdir

proj_dir = args.outdir+"/"+args.name

umi_reads=dict()
filenames=dict()
filenames["Settings"]=proj_dir+"/Umi_AnalysisSettings.txt"
filenames["Fastq1"]=args.fastq_file1
filenames["Fastq2"]=args.fastq_file2
filenames["CollapsedR1FastQ"]=proj_dir+"/"+args.name+".collapsedR1.fastq.gz"
filenames["CollapsedR2FastQ"]=proj_dir+"/"+args.name+".collapsedR2.fastq.gz"
filenames["TrimCollapsedR1FastQ"]=proj_dir+"/"+args.name+".trimmed.collapsedR1.fastq.gz"
filenames["TrimCollapsedR2FastQ"]=proj_dir+"/"+args.name+".trimmed.collapsedR2.fastq.gz"
filenames["ConcatenateFastQ"]=proj_dir+"/"+args.name+".extendedFrags.fastq.gz"
filenames["collapseSAM"]=proj_dir+"/"+args.name+".collapsed.sam"
filenames["collapseBAM"]=proj_dir+"/"+args.name+".collapsed.bam"
filenames["uncollapseSAM"]=proj_dir+"/"+args.name+".uncollapsed.sam"
filenames["uncollapseBAM"]=proj_dir+"/"+args.name+".uncollapsed.bam"
filenames["collapseVCF"]=proj_dir+"/"+args.name+".collapsed.vcf"
filenames["uncollapseVCF"]=proj_dir+"/"+args.name+".uncollapsed.vcf"
####################################################################################################################
# Prepare working directory and output
if not os.path.exists(proj_dir):
    os.makedirs(proj_dir)
saveSettings(args,filenames["Settings"])
####################################################################################################################
# Step1 Import Reads and extract UMI
count = 0
# Fastq 1 file
print "Step 1 a) Read FastQ1 file"
with gzip.open(filenames["Fastq1"], 'r') as f:
    col_count = 0
    for line in f:
        if col_count % 4 == 0:
            line1 = line
        elif col_count % 4 == 1:
            line2 = line 
        elif col_count % 4 == 2:
            line3 = line
        elif col_count % 4 == 3:
            line4 = line
            count += 1
            read = UMIRead()
            read.createFromRead1(line1[1:].split(" ")[0],line2,line4,args.umi1,args.off1)
            umi_reads[read.getName()] = read
        col_count += 1

# Fastq 2 file
print "Step 1 b) Read FastQ2 file"
with gzip.open(filenames["Fastq2"], 'r') as f:
    col_count = 0
    for line in f:
        if col_count % 4 == 0:
            line1 = line
        elif col_count % 4 == 1:
            line2 = line 
        elif col_count % 4 == 2:
            line3 = line
        elif col_count % 4 == 3:
            line4 = line
            count += 1
            read = umi_reads[line1[1:].split(" ")[0]]
            read.createFromRead2(line1[1:].split(" ")[0],line2,line4,args.umi2,args.off2)
            read.combineUMIs()
        col_count += 1

print "  "+str(len(umi_reads))+" reads found"

####################################################################################################################
# Step2 create Read groups
print "Step 2) Create Read Groups"

umis=list(UMIRead.umi_list)

read_groups = dict.fromkeys(umis)
read_groups.update(dict(read_groups))

for read in umi_reads.values():
    umi = read.getUMI()
    read_groups[umi]=UMIReadGroup(umi,read)

####################################################################################################################
# Step3 Get ReadGroup Statistics
print "Step 3) Get Read Group Distribution"
read_group_members = list()
read_group_count = 0
read_group_single = 0
read_group_2mem = 0
read_group_3mem = 0
read_group_4mem = 0
read_group_5mem = 0
read_group_6_10mem = 0
read_group_10_30mem = 0
read_group_over30mem = 0
read_group_over100mem = 0
read_group_threshold = 0
for key in read_groups.keys():
    read_group_count += 1
    count = read_groups[key].getCount()
    read_group_members.append(count)
    if count == 1:
        read_group_single += 1
    if count == 2:
        read_group_2mem += 1
    if count == 3:
        read_group_3mem += 1
    if count == 4:
        read_group_4mem += 1
    if count == 5:
        read_group_5mem += 1
    if count > 5 and count <= 10:
        read_group_6_10mem += 1
    if count > 10 and count <= 30:
        read_group_10_30mem += 1
    if count > 30:
        read_group_over30mem += 1
    if count > 100:
        read_group_over100mem += 1
    if count >= args.member_threshold:
        read_group_threshold += 1


plt.hist(read_group_members,facecolor='green',range=[1,1000])
plt.savefig(proj_dir+'/ReadGroup_Hist_large.pdf')
plt.clf()
plt.hist(read_group_members,facecolor='green',bins=20,range=[1,100])
plt.savefig(proj_dir+'/ReadGroup_Hist_small.pdf')

STATS=open(proj_dir+"/ReadGroupStats.txt","w")
STATS.write("ReadGroup Statistics\n")
STATS.write("--------------------\n\n")
STATS.write("Reads: "+str(len(umi_reads.keys()))+"\n")
STATS.write("Read Groups: "+str(read_group_count)+"\n")
STATS.write("Member Threshold: "+str(args.member_threshold)+"\n")
STATS.write("Read groups with more members than threshold: "+str(read_group_threshold)+"\n")
STATS.write("Read Groups with 1 member: "+str(read_group_single)+"\n")
STATS.write("Read Groups with 2 members: "+str(read_group_2mem)+"\n")
STATS.write("Read Groups with 3 members: "+str(read_group_3mem)+"\n")
STATS.write("Read Groups with 4 members: "+str(read_group_4mem)+"\n")
STATS.write("Read Groups with 5 members: "+str(read_group_5mem)+"\n")
STATS.write("Read Groups with >5 but < 11 members: "+str(read_group_6_10mem)+"\n")
STATS.write("Read Groups with >10 but < 31 members: "+str(read_group_10_30mem)+"\n")
STATS.write("Read Groups with >30 members: "+str(read_group_over30mem)+"\n")
STATS.write("Read Groups with >100 members: "+str(read_group_5mem)+"\n")
STATS.write("Read Groups mean member count: "+str(np.mean(read_group_members))+"\n")
STATS.write("Read Groups median member count: "+str(np.median(read_group_members))+"\n")
STATS.write("Read Groups minimum member count: "+str(np.amin(read_group_members))+"\n")
STATS.write("Read Groups maximum member count: "+str(np.amax(read_group_members))+"\n")
STATS.close()
print "  "+str(len(read_groups))+" read groups found"

read_groups[read_groups.keys()[10]].printFASTA("test.fa")

####################################################################################################################
# Step4 Build Consensus Sequence of Read Group
print "Step 4) Build Consensus sequences"
for key in read_groups.keys():
    read_groups[key].createConsensusR1(args.consensus)
    read_groups[key].createConsensusR2(args.consensus)


####################################################################################################################
# Step5 Output consensus sequences as Fastq Files
print "Step 5) Output consensus sequences as Fastq Files"
with gzip.open(filenames["CollapsedR1FastQ"], 'wb') as f:
    for key in read_groups.keys():
        if read_groups[key].getCount < args.member_threshold:
            continue
        line1=read_groups[key].getName()+"_"+str(read_groups[key].getCount())
        line2=read_groups[key].getConsSequenceR1()
        line3="+"
        line4="I"*len(read_groups[key].getConsSequenceR1())
        f.write(line1+"\n"+line2+"\n"+line3+"\n"+line4+"\n")
with gzip.open(filenames["CollapsedR2FastQ"], 'wb') as f:
    for key in read_groups.keys():
        if read_groups[key].getCount < args.member_threshold:
            continue
        line1=read_groups[key].getName()+"_"+str(read_groups[key].getCount())
        line2=read_groups[key].getConsSequenceR2()
        line3="+"
        line4="I"*len(read_groups[key].getConsSequenceR2())
        f.write(line1+"\n"+line2+"\n"+line3+"\n"+line4+"\n")

####################################################################################################################
# Trim adaptors if specified
if args.adaptor1 != "":
    if args.adaptor2 != "":
        subprocess.call([script_dir+"/Software/cutadapt","-a",args.adaptor1,"-A",args.adaptor2,"-o",filenames["TrimCollapsedR1FastQ"],"-p",
              filenames["TrimCollapsedR2FastQ"],filenames["CollapsedR1FastQ"],filenames["CollapsedR2FastQ"]])
    else:
        subprocess.call([script_dir+"/Software/cutadapt","-a",args.adaptor1,"-o",filenames["TrimCollapsedR1FastQ"],filenames["CollapsedR1FastQ"]])
####################################################################################################################
# Step6 Concatenate Forward and reverse reads using FLASH
print "Step 6) Concatenate Forward and reverse reads using FLASH"
if args.adaptor1 != "":
    if args.adaptor2 != "":
        subprocess.call([script_dir+"/Software/flash","-z","-M","160","-p","33","-x","0.5","-o",proj_dir+"/"+args.name,filenames["TrimCollapsedR1FastQ"],
                filenames["TrimCollapsedR2FastQ"]])
    else:
        subprocess.call([script_dir+"/Software/flash","-z","-M","160","-p","33","-x","0.5","-o",proj_dir+"/"+args.name,filenames["TrimCollapsedR1FastQ"],
                filenames["CollapsedR2FastQ"]])
else:
    subprocess.call([script_dir+"/Software/flash","-z","-M","160","-p","33","-x","0.5","-o",proj_dir+"/"+args.name,filenames["CollapsedR1FastQ"],
                filenames["CollapsedR2FastQ"]])
####################################################################################################################
# Step7 Align to Human Reference Genome
if args.skip_aln:
    print "Stop Analysis now"
    sys.exit()
print "Step 7) Alignment to Human Reference Genome"
SAM = open(filenames["collapseSAM"],"w")
subprocess.call([script_dir+"/Software/bwa","mem","-M","-R","@RG\tID:"+args.name+"\tSM:"+args.name+"\tPL:illumina","-t",str(args.threads),script_dir+"/Ref/hg19",filenames["ConcatenateFastQ"]],stdout=SAM)
SAM.close()
bam_conversion=subprocess.call(["java","-jar",script_dir+"/Software/picard.jar","SortSam","SO=coordinate","INPUT="+filenames["collapseSAM"],"OUTPUT="+filenames["collapseBAM"],"VALIDATION_STRINGENCY=LENIENT","CREATE_INDEX=true"])
subprocess.call(["rm",filenames["collapseSAM"]])

#get Statistics on BAM file
environ=os.environ.copy()
environ["LD_LIBRARY_PATH"]=script_dir+"/Software"
BAMSTAT = open(filenames["collapseBAM"]+".stats","w")
subprocess.call([script_dir+"/Software/bamtools","stats","-in",filenames["collapseBAM"]],stdout=BAMSTAT,env=environ)
BAMSTAT.close()

####################################################################################################################
# Step8 Align Uncollapsed FastQs to Human Reference Genome
print "Step 8) Align Uncollapsed FastQs to Human Reference Genome"
SAM = open(filenames["uncollapseSAM"],"w")
subprocess.call([script_dir+"/Software/bwa","mem","-M","-R","@RG\tID:"+args.name+"\tSM:"+args.name+"\tPL:illumina","-t",str(args.threads),script_dir+"/Ref/hg19",filenames["Fastq1"],filenames["Fastq2"]],stdout=SAM)
SAM.close()
bam_conversion=subprocess.call(["java","-jar",script_dir+"/Software/picard.jar","SortSam","SO=coordinate","INPUT="+filenames["uncollapseSAM"],"OUTPUT="+filenames["uncollapseBAM"],"VALIDATION_STRINGENCY=LENIENT","CREATE_INDEX=true"])
subprocess.call(["rm",filenames["uncollapseSAM"]])

#get Statistics on BAM file
BAMSTAT = open(filenames["uncollapseBAM"]+".stats","w")
subprocess.call([script_dir+"/Software/bamtools","stats","-insert","-in",filenames["uncollapseBAM"]],stdout=BAMSTAT,env=environ)
BAMSTAT.close()

####################################################################################################################
# Step9 SNP calling on collapsed BAM
if args.skip_snp:
    print "Stop Analysis now"
    sys.exit()
print "Step 9) SNP calling MuTect on collapsed BAM"
if args.target_regions != "":
    subprocess.call(["java","-jar",script_dir+"/Software/GenomeAnalysisTK.jar","--analysis_type","MuTect2","--reference_sequence",
           script_dir+"/Ref/hg19.fa","--input_file:tumor",filenames["collapseBAM"],"--out",filenames["collapseVCF"],"-L",
           args.target_regions,"-rf","BadCigar"])
else:
    print "SNP calling skipped"

####################################################################################################################
# Step10 SNP calling on uncollapsed BAM
print "Step 10) SNP calling MuTect on uncollapsed BAM"
if args.target_regions != "":
    subprocess.call(["java","-jar",script_dir+"/Software/GenomeAnalysisTK.jar","--analysis_type","MuTect2","--reference_sequence",
           script_dir+"/Ref/hg19.fa","--input_file:tumor",filenames["uncollapseBAM"],"--out",filenames["uncollapseVCF"],"-L",
           args.target_regions,"-rf","BadCigar","-nct",str(args.threads)])
else:
    print "SNP calling skipped"
