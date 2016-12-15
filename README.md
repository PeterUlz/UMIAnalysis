# UMI Analysis
Analyze resequencing experiments containing UMIs

Preparation:
-) install java

Basic Usage:
usage: UMI_analysis.py [-h] -fq1 FASTQ_FILE1 -fq2 FASTQ_FILE2 -umi1 UMI1 [-umi2 UMI2] [-off1 OFF1] [-off2 OFF2][-members MEMBER_THRESHOLD] [-cons CONSENSUS] -s NAME [-o OUTDIR] [-t THREADS]

Note:
This calls every step of the below list. UMIs are extracted from FastQ1 and (possibly) from FastQ2. Offsets (Off1) may be specified: This trims away a certain amount of basepairs after the UMI sequence.

*Attention: Large hg19 reference index files (with pseudo-autosomal region masked) for BWA are not in this commit*

CNVs of ctDNA low-coverage WGS sequencing are analyzed in several steps

 1) Import Reads and extract UMI  
 2) Create Read Groups  
 3) Get Read Group Distribution  
 4) Build Consensus sequences  
 5) Output consensus sequences as collapsed Fastq Files  
 6) Concatenate collapsed forward and reverse reads using FLASH  
 7) Alignment (collapsed and uncollapsed reads) to Human Reference Genome  

Settings:
* SafeSeq:
  * UMI1: 14
  * UMI2: 0
  * OFF1: 0
  * OFF2: 0
  * Consensus: 0.5
  * Members: 5

* Thruplex:
  * UMI1: 6
  * UMI2: 6
  * OFF1: 11
  * OFF2: 11
  * Consensus: 0.5
  * Members: 5

* smMIP:
  * UMI1: 5
  * UMI2: 5
  * OFF1: 0
  * OFF2: 0
  * Consensus: 0.5
  * Members: 1
  * Adapter Trimming

* QiaSeq:
  * UMI1: 12
  * UMI2: 12
  * OFF1: 0
  * OFF2: 0
  * Consensus: 0.5
  * Members: 1

* NebNext:
