#!/usr/bin/python
import os
import sys
import glob
import string
from collections import Counter

def hashincounttable(filename):
  MutDict = {}
  infile = open(filename,'r')
  countline = 0
  for line in infile:
    if 'count' in line: continue
    line = line.rstrip().rsplit("\t")
    MutDict[line[0]] = int(line[1])
  infile.close()
  return MutDict

def FreqCal(Dict, mut, Total):
  freq = float(Dict[mut])/float(Total)*1000000 if Dict.has_key(mut) else float(0)
  return freq

def main():
  outfile    = 'data/VariantFreqTable.tsv'
  DNADict    = hashincounttable('count/count_DNA.tsv')
  R0Dict     = hashincounttable('count/count_R0.tsv')
  WT_R1Dict   = hashincounttable('count/count_WT-R1.tsv')
  WT_R2Dict   = hashincounttable('count/count_WT-R2.tsv')
  WT_R3Dict   = hashincounttable('count/count_WT-R3.tsv')
  I45M_R1Dict   = hashincounttable('count/count_I45M-R1.tsv')
  I45M_R2Dict   = hashincounttable('count/count_I45M-R2.tsv')
  I45M_R3Dict   = hashincounttable('count/count_I45M-R3.tsv')
  I45T_R1Dict   = hashincounttable('count/count_I45T-R1.tsv')
  I45T_R2Dict   = hashincounttable('count/count_I45T-R2.tsv')
  I45T_R3Dict   = hashincounttable('count/count_I45T-R3.tsv')
  I45F_R1Dict   = hashincounttable('count/count_I45F-R1.tsv')
  I45F_R2Dict   = hashincounttable('count/count_I45F-R2.tsv')
  I45F_R3Dict   = hashincounttable('count/count_I45F-R3.tsv')
  muts       = list(set(R0Dict.keys()+WT_R1Dict.keys()+WT_R2Dict.keys()+WT_R3Dict.keys()+
                        I45M_R1Dict.keys()+I45M_R2Dict.keys()+I45M_R3Dict.keys()+
                        I45T_R1Dict.keys()+I45T_R2Dict.keys()+I45T_R3Dict.keys()+
                        I45F_R1Dict.keys()+I45F_R2Dict.keys()+I45F_R3Dict.keys()))
  R0Total   = sum(R0Dict.values())
  WT_R1Total = sum(WT_R1Dict.values())
  WT_R2Total = sum(WT_R2Dict.values())
  WT_R3Total = sum(WT_R3Dict.values())
  I45M_R1Total = sum(I45M_R1Dict.values())
  I45M_R2Total = sum(I45M_R2Dict.values())
  I45M_R3Total = sum(I45M_R3Dict.values())
  I45T_R1Total = sum(I45T_R1Dict.values())
  I45T_R2Total = sum(I45T_R2Dict.values())
  I45T_R3Total = sum(I45T_R3Dict.values())
  I45F_R1Total = sum(I45F_R1Dict.values())
  I45F_R2Total = sum(I45F_R2Dict.values())
  I45F_R3Total = sum(I45F_R3Dict.values())

  R0min  = float(1)/float(R0Total)*1000000
  outfile     = open(outfile,'w')
  outfile.write("\t".join(['Variant','R0','WT_R1','WT_R2','WT_R3',
                           'I45M_R1','I45M_R2','I45M_R3',
                           'I45T_R1','I45T_R2','I45T_R3',
                           'I45F_R1','I45F_R2','I45F_R3'])+"\n")
  mutcount   = 0
  for mut in muts:
    mutcount += 1
    if mutcount%100000 == 0: print 'Processed %i variants' % mutcount
    R0freq    = FreqCal(R0Dict, mut, R0Total)
    WT_R1freq  = FreqCal(WT_R1Dict, mut, WT_R1Total)
    WT_R2freq  = FreqCal(WT_R2Dict, mut, WT_R2Total)
    WT_R3freq  = FreqCal(WT_R3Dict, mut, WT_R3Total)
    I45M_R1freq  = FreqCal(I45M_R1Dict, mut, I45M_R1Total)
    I45M_R2freq  = FreqCal(I45M_R2Dict, mut, I45M_R2Total)
    I45M_R3freq  = FreqCal(I45M_R3Dict, mut, I45M_R3Total)
    I45T_R1freq  = FreqCal(I45T_R1Dict, mut, I45T_R1Total)
    I45T_R2freq  = FreqCal(I45T_R2Dict, mut, I45T_R2Total)
    I45T_R3freq  = FreqCal(I45T_R3Dict, mut, I45T_R3Total)
    I45F_R1freq  = FreqCal(I45F_R1Dict, mut, I45F_R1Total)
    I45F_R2freq  = FreqCal(I45F_R2Dict, mut, I45F_R2Total)
    I45F_R3freq  = FreqCal(I45F_R3Dict, mut, I45F_R3Total)
    outfile.write("\t".join(map(str,[mut,R0freq,WT_R1freq,WT_R2freq,WT_R3freq,
                                     I45M_R1freq,I45M_R2freq,I45M_R3freq,
                                     I45T_R1freq,I45T_R2freq,I45T_R3freq,
                                     I45F_R1freq,I45F_R2freq,I45F_R3freq,]))+"\n")
  outfile.close()

if __name__ == "__main__":
  main()
