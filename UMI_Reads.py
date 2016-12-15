# UMI Module

class UMIRead():
   name = ""
   sequence_r1 = ""
   quality_r1 = 0
   sequence_r2 = ""
   quality_r2 = 0
   umi1 = ""
   umi2 = ""
   combined_umi = ""
   
   def createFromRead1(self,r1_name,r1_sequence,r1_quality,umi1_length,offset1):
       self.name = r1_name.rstrip()
       self.umi1 = r1_sequence[:umi1_length]
       self.sequence_r1 = r1_sequence[(umi1_length+offset1):].rstrip()
       self.quality = r1_quality[(umi1_length+offset1):].rstrip()


   def createFromRead2(self,r2_name,r2_sequence,r2_quality,umi2_length,offset2):
       self.name = r2_name.rstrip()
       self.umi2 = r2_sequence[:umi2_length]
       self.sequence_r2= r2_sequence[(umi2_length+offset2):].rstrip()
       self.quality_r2 = r2_quality[(umi2_length+offset2):].rstrip()

   def combineUMIs(self):
       self.combined_umi = self.umi1 + self.umi2

   def printReadSummary(self):
       print self.name+" "+self.sequence_r1
   
   def getName(self):
       return self.name

   def getUMI(self):
       return self.combined_umi

   def getSequenceR1(self):
       return self.sequence_r1
   def getSequenceR2(self):
       return self.sequence_r2





######################################################################################
class UMIReadGroup():
   name = ""
   read_list = list()
   cons_sequence_r1 = ""
   cons_sequence_r2 = ""
   umi = ""
   count = 0
   max_size_r1 = 0
   max_size_r2 = 0

   def __init__(self,umi_comb,read):
       self.name = umi_comb
       self.umi = umi_comb
       self.read_list=list()
       self.read_list.append(read)
       self.count = 1
       self.max_size_r1 = len(read.getSequenceR1())
       self.max_size_r2 = len(read.getSequenceR2())


   def getCount(self):
       return self.count

   def getLenReadList(self):
       return len(self.read_list)
   def getName(self):
       return self.name
   def getConsSequenceR1(self):
       return self.cons_sequence_r1
   def getConsSequenceR2(self):
       return self.cons_sequence_r2
######################################################################################
   def addRead(self,read):
       self.read_list.append(read)
       self.count += 1
       if len(read.getSequenceR1()) > self.max_size_r1:
           self.max_size_r1 = len(read.getSequenceR1())
       if len(read.getSequenceR2()) > self.max_size_r2:
           self.max_size_r2 = len(read.getSequenceR2())
 ######################################################################################
   def getCount(self):
       return self.count
 ######################################################################################
   def printFASTA(self,filename):
       FASTA = open(filename,"w")
       for read in self.read_list:
           FASTA.write(">"+read.getName()+"\n"+read.getSequenceR1()+"\n")
       FASTA.close()
######################################################################################
   def createConsensusR1(self,consensus_thresh):
       consensus = ""
       for i in range (0,self.max_size_r1):
           bases_seen = list()
           for read in self.read_list:
              seq = read.getSequenceR1()
              if i < len(seq):
                  bases_seen.append(seq[i])
           a_count = bases_seen.count('A')
           t_count = bases_seen.count('T')
           g_count = bases_seen.count('G')
           c_count = bases_seen.count('C')
           if a_count > consensus_thresh * float(len(bases_seen)):
               consensus += "A"
           elif t_count > consensus_thresh * float(len(bases_seen)):
               consensus += "T"
           elif g_count > consensus_thresh * float(len(bases_seen)):
               consensus += "G"
           elif c_count > consensus_thresh * float(len(bases_seen)):
               consensus += "C"
           else:
               consensus += "N"
       self.cons_sequence_r1 = consensus
######################################################################################
   def createConsensusR2(self,consensus_thresh):
       consensus = ""
       for i in range (0,self.max_size_r2):
           bases_seen = list()
           for read in self.read_list:
              seq = read.getSequenceR2()
              if i < len(seq):
                  bases_seen.append(seq[i])
           a_count = bases_seen.count('A')
           t_count = bases_seen.count('T')
           g_count = bases_seen.count('G')
           c_count = bases_seen.count('C')
           if a_count > consensus_thresh * float(len(bases_seen)):
               consensus += "A"
           elif t_count > consensus_thresh * float(len(bases_seen)):
               consensus += "T"
           elif g_count > consensus_thresh * float(len(bases_seen)):
               consensus += "G"
           elif c_count > consensus_thresh * float(len(bases_seen)):
               consensus += "C"
           else:
               consensus += "N"
       self.cons_sequence_r2 = consensus
           
       
