#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import sys
#import csv
from local_align import local_align, ScoreParam

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
            # print(read)
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
    leftseq = wtSeq[15:30]
    rightseq = wtSeq[-30:-15]
    file1= glob.glob(prefix + "*R1*")[0]
    file2 = glob.glob(prefix + "*R2*")[0]
    print("files are ", file1, file2)
    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1, mode='rt') as f1:
            with gzip.open(file2, mode='rt') as f2:
                R1R2, nreads = interleave(f1, f2)
    else:
        with open(file1) as f1:
            with open(file2) as f2:
                R1R2, nreads = interleave(f1, f2)
    # check the 
    seqList = []
    countList = []
    for k in sorted(R1R2, key=R1R2.get, reverse=True):
        (R1, R2) = k
        if leftseq not in R1 or rightseq not in R2:
            continue
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
    algnList = [] # alignment of r1 and r2
    algnList1 = [] # alignment of wt and r1
    algnList2 = [] # alignment of wt and r2
    indexList = [] # list of index in seqList used
    nread2 = 0 # reads for on target amplicons
    for i in range(len(seqList)):
        nread2 += countList[i]
        # if countList[i] * 100 / nreads < 5:
        #     break # seqList is sorted
        (r1, r2) = seqList[i]
        # if leftseq not in r1 or rightseq not in r2:
        #     continue
        # add i to index list
        indexList.append(i)
        # check whether there is overlap between R1 and R2
        a, b = local_align(r1, r2, ScoreParam(3, -4, -5))
        algnList.append((a,b))
        # if a and ("-" not in a and "-" not in b): # overlap
        #     # print(i, a)
        #     p1 = r1.find(a)
        #     p2 = r2.find(b)
        #     merged = r1[0:p1] + r2[(p2+len(b)):]
        #     c, d = local_align(wtSeq, merged, ScoreParam(3, -4, -5))
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
        c, d = local_align(wtSeq, r1, ScoreParam(3, -4, -5))
        if "-" in c: # indel
            xx = c.find("-")
            indelPos1 = wtSeq.find(c[:xx]) + xx
        elif "-" in d: # indel
            xx = d.find("-")
            indelPos1 = wtSeq.find(c[:xx]) + xx
        else:
            indelPos1 = -1
        algnList1.append((c,d))
        # r2
        e, f = local_align(wtSeq, r2, ScoreParam(3, -4, -5))
        if "-" in e: # indel
            xx = e.find("-")
            indelPos2 = wtSeq.find(e[:xx]) + xx
        elif "-" in f: # indel
            xx = f.find("-")
            indelPos2 = wtSeq.find(e[:xx]) + xx
        else:
            indelPos2 = -1
        indelPosList.append((indelPos1, indelPos2))
        algnList2.append((e,f))

    return [nreads, nread2, seqList, algnList, algnList1, algnList2, countList, indelPosList, indexList]




def main():
    ref_lib = get_fasta(os.path.abspath(sys.argv[1]))
    ID = sys.argv[2]
    ref_seq = ref_lib[ID].upper()
    prefix = sys.argv[3]
    nreads, nread2, seqList, algnList, algnList1, algnList2, countList, indelPosList, indexList = checkFastq(prefix, ref_seq)
    print("Total reads is ", nreads)
    print("Reads on target are ", nread2)
    print("Length of index, algnList, algnList2, seqList, countList, indelPosList", len(indexList), len(algnList), len(algnList2), len(seqList), len(countList), len(indelPosList))
    print("PCR amplicon length is ", len(ref_seq))
    print("PE reads\talignment between R1 and R2\talignment length between R1 and R2\talignment between WT and R1\talignment length between WT and R1\talignment between WT and R2\talignment length between WT and R2\tNumber of reads\tIndel position in R1 and R2")
    for i in range(len(algnList2)):
        print(seqList[indexList[i]], algnList[i], 0 if "-" in algnList[i][0] else len(algnList[i][0]), algnList1[i], len(algnList1[i][0]), algnList2[i], len(algnList2[i][0]), countList[indexList[i]], indelPosList[i], sep='\t')

if __name__ == "__main__":
    main()
