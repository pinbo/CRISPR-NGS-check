#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import sys
#import csv
from local_align import *

# added by Junli Zhang on 11/30/2019
# function to extract sequences from a fasta file


def get_fasta(infile):
    fasta = {}  # dictionary for alignment
    with open(infile) as file_one:
        for line in file_one:
            line = line.strip()
            if line:  # skip blank lines
                if line.startswith(">"):
                    # left strip > or space, so " > abc edf" will be "abc edf", then split by space to get "abc"
                    sequence_name = line.lstrip("> ").split()[0]
                    fasta[sequence_name] = ""
                else:
                    # remove spaces in case
                    fasta[sequence_name] += line.replace(" ", "")
    return fasta

# reverse complement a sequence


def ReverseComplement(seq):
    s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
    s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
    seq_dict = {s1[i]: s2[i] for i in range(len(s1))}
    return "".join([seq_dict[base] for base in reversed(seq)])

# interleave R1 and R2
def interleave(f1, f2):
    out = {}  # a dict of reads count
    n = 0 # number of reads
    while True:
        n += 1
        line = f1.readline()
        if line.strip() == "":
            break
        # print line.strip()
        for i in range(3):
            read = f1.readline().strip()
            if i == 0:
                R1 = read
        for i in range(4):
            read = f2.readline().strip()
            if i == 1:
                R2 = ReverseComplement(read) # RC of read2
        # put into the dict
        if (R1, R2) in out:
            out[(R1, R2)] += 1
        else:
            out[(R1, R2)] = 1
    return [out, n]
# check two sequences with only 1 SNP
def checkSNP(seq1,seq2):
    l1 = len(seq1)
    l2 = len(seq2)
    if l1 != l2:
        return 100 # big difference
    else:
        ndiff = 0
        for i in range(l1):
            if seq1[i] != seq2[i]:
                ndiff += 1
        return ndiff

# check PE fastq
def checkFastq(prefix, wtSeq): # wtSeq is the PCR amplicon of the unedited template 
    R1file = glob.glob(prefix + "*R1*.fastq")[0]
    R2file = glob.glob(prefix + "*R2*.fastq")[0]
    with open(R1file) as f1:
        with open(R2file) as f2:
            R1R2, nreads = interleave(f1, f2)
    # check the 
    seqList = []
    countList = []
    for k in sorted(R1R2, key=R1R2.get, reverse=True):
        (R1, R2) = k
        if seqList:
            new = 1
            for i in range(len(seqList)):
                (r1, r2) = seqList[i]
                if checkSNP(r1, R1) < 3 and checkSNP(r2, R2) < 3: # less than 3 SNPs, treat as the same
                    countList[i] += 1
                    new = 0
                    break
            if new:
                seqList.append(k)
                countList.append(R1R2[k])
        else:
            seqList.append(k)
            countList.append(R1R2[k])
    # now check the indels
    indelPosList = []
    # algnList = [] # alignment of r1 and r2
    algnList2 = [] # alignment of wt and r1, wt and r2
    for i in range(len(seqList)):
        if countList[i] * 100 / nreads < 10:
            break
        (r1, r2) = seqList[i]
        # check whether there is overlap between R1 and R2
        # a, b = local_align(r1, r2, ScoreParam(3, -3, -6))
        # algnList.append((a,b))
        # if a and ("-" not in a and "-" not in b): # overlap
        #     print(i, a)
        #     p1 = r1.find(a)
        #     p2 = r2.find(b)
        #     merged = r1[0:p1] + r2[(p2+len(b)):]
        #     c, d = local_align(wtSeq, merged, ScoreParam(3, -3, -5))
        #     if "-" in c: # indel
        #         xx = c.find("-")
        #         indelPos = wtSeq.find(c[:xx]) + xx
        #     elif "-" in d: # indel
        #         xx = d.find("-")
        #         indelPos = wtSeq.find(c[:xx]) + xx
        #     else:
        #         indelPos = -1
        #     indelPosList.append(indelPos)
        #     algnList2.append((c,d))
        # else: # no overlap, need to check r1 and r2 separately
        # r1 first
        c, d = local_align(wtSeq, r1, ScoreParam(3, -3, -5))
        if "-" in c: # indel
            xx = c.find("-")
            indelPos1 = wtSeq.find(c[:xx]) + xx
        elif "-" in d: # indel
            xx = d.find("-")
            indelPos1 = wtSeq.find(c[:xx]) + xx
        else:
            indelPos1 = -1
        # r2
        e, f = local_align(wtSeq, r2, ScoreParam(3, -3, -5))
        if "-" in e: # indel
            xx = e.find("-")
            indelPos2 = wtSeq.find(e[:xx]) + xx
        elif "-" in f: # indel
            xx = f.find("-")
            indelPos2 = wtSeq.find(e[:xx]) + xx
        else:
            indelPos2 = -1
        indelPosList.append((indelPos1, indelPos2))
        algnList2.append((c,d,e,f))

    return [nreads, seqList, algnList2, countList, indelPosList]




def main():
    ref_lib = get_fasta(os.path.abspath(sys.argv[1]))
    ID = sys.argv[2]
    ref_seq = ref_lib[ID].upper()
    prefix = "3_"
    nread, seqList, algnList2, countList, indelPosList = checkFastq(prefix, ref_seq)
    print("Total reads is ", nread)
    for i in range(len(algnList2)):
        print(seqList[i], algnList2[i], countList[i], indelPosList[i])

if __name__ == "__main__":
    main()
