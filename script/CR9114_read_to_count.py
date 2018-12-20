#!/usr/bin/python
import os
import sys
import glob
import string
from Bio import SeqIO
from collections import Counter

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def readingFile2SampleID(file2SampleID_file):
  file2SampleID_dict = {}
  infile = open(file2SampleID_file,'r')
  for line in infile.xreadlines():
    if 'R1File' in line: continue
    line = line.rstrip().rsplit("\t")
    file2SampleID_dict[line[0]] = line[1]
  infile.close()
  return file2SampleID_dict

def Processlib(R1file):
  R2file = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file, "fastq")
  R2records = SeqIO.parse(R2file,"fastq") 
  muts = []
  count_record = 0
  for R1record in R1records:
    count_record += 1
    #if count_record == 10000: break
    R2record = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    Mutcodon_dict = {
      'R1resi31': str(R1seq[26:29]),   'R2resi98': str(rc(R2seq[26:29])),
      'R1resi52a': str(R1seq[92:95]),  'R2resi54': str(rc(R2seq[167:170])),
      'R1resi53': str(R1seq[95:98]),   'R2resi53': str(rc(R2seq[170:173])),
      'R1resi54': str(R1seq[98:101]),  'R2resi52a': str(rc(R2seq[173:176])),
      'R1resi98': str(R1seq[239:242]), 'R2resi31': str(rc(R2seq[239:242])),
      }
    R1roi = Mutcodon_dict['R1resi31']+Mutcodon_dict['R1resi52a']+Mutcodon_dict['R1resi53']+ \
            Mutcodon_dict['R1resi54']+Mutcodon_dict['R1resi98']
    R2roi = Mutcodon_dict['R2resi31']+Mutcodon_dict['R2resi52a']+Mutcodon_dict['R2resi53']+ \
            Mutcodon_dict['R2resi54']+Mutcodon_dict['R2resi98']
    if R1roi != R2roi: continue
    if len(R1roi) != 15: continue
    if 'N' in R1roi: continue
    if str(R1seq[23:26]) != 'AAC': continue
    if str(R1seq[29:32]) != 'TAC': continue
    if str(R1seq[89:92]) != 'TCT': continue
    if str(R1seq[101:104]) != 'GGG': continue
    if str(R1seq[236:239]) != 'AAC': continue
    if str(R1seq[242:245]) != 'TAT': continue
    mut = translation(R1roi)
    if mut == 'NPIFY':
      muts.append(mut)
    else:
      if R1roi[2] not in ['G','C']: continue
      if R1roi[5] not in ['G','C']: continue
      if R1roi[8] not in ['G','C']: continue
      if R1roi[11] not in ['G','T']: continue
      if R1roi[14] not in ['G','C']: continue
      muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def Output(mutcount_dict, outfile):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','count'])+"\n")
  for mut in mutcount_dict.keys():
    outfile.write("\t".join(map(str,[mut,mutcount_dict[mut]]))+"\n")
  outfile.close()

def main():
  R1filenames = glob.glob('fastq/*R1*.fastq')
  file2SampleID_dict = readingFile2SampleID('doc/SampleID.tsv')
  print "Analyzing %i samples" % len(R1filenames)
  for R1filename in sorted(R1filenames,key=lambda x:int(x.rsplit('-')[1].rsplit('_')[0])):
    SampleID = file2SampleID_dict[R1filename]
    outfile  = 'count/count_'+SampleID+'.tsv'
    print "Analyzing %s (%s)" % (R1filename, SampleID)
    mutcount_dict = Processlib(R1filename)
    Output(mutcount_dict, outfile)

if __name__ == "__main__":
  main()
