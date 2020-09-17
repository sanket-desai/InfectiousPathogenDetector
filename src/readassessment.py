#!/usr/bin/python3

'''

'''
import pysam
from globals import *

#Assess read based on quality, alignment, sequence complexity,
class ShortReadAssessment(object):
    def __init__(self, samread):
        self.alignedread_=samread
        self.meanphread_=sum(samread.query_qualities)/len(samread.query_qualities)
        self.nuccounts_={}
        atc=0
        gcc=0
        for i in range(len(samread.query_sequence)-1):
            if samread.query_sequence[i:i+2] == "AT":
                atc+=1
            elif samread.query_sequence[i:i+2] == "GC":
                gcc+=1
        self.readlength_=len(samread.query_sequence)
        self.nuccounts_["AT"]=atc
        self.nuccounts_["GC"]=gcc
        self.nuccounts_["A"]=samread.query_sequence.count("A")
        self.nuccounts_["T"]=samread.query_sequence.count("T")
        self.nuccounts_["C"]=samread.query_sequence.count("C")
        self.nuccounts_["G"]=samread.query_sequence.count("G")
    def is_low_complexity(self):
        lc=False
        if (self.nuccounts_["GC"] / self.readlength_)  >= 0.89  :
            lc=True
        elif (self.nuccounts_["AT"] / self.readlength_) >= 0.87 :
            lc=True
        elif (self.nuccounts_["A"] / self.readlength_) >= 0.6 or (self.nuccounts_["T"]/self.readlength_) >= 0.6 or (self.nuccounts_["G"] / self.readlength_) >= 0.6 or (self.nuccounts_["C"]/self.readlength_) >= 0.6:
            lc=True
        return lc
    def is_low_length(self, minlen=40):
        return self.readlength_ < minlen
    def is_low_confidence_alignment(self): #use match / mismatch ratio instead // hard cut off for the mismatch could be 10
        lca=False
        t=self.alignedread_.cigartuples
        #Modify CIGAR tuples by only including the non-clipped portion of the read to calculate the mismatch / read length ratio
        if t[0][0]==4 and t[0][1] <= 15:
            t=t[1:]
        if t[-1][0]==4 and t[-1][1] <= 15:
            t=t[:-1]
        curlen=0
        matchlen=0
        mismatch=0
        for ti in t:
            curlen+=ti[1]
            if ti[0]==0:
                matchlen+=ti[1]
            elif ti[0] == 8:
                mismatch+=ti[1]
#        if (matchlen / curlen) < 0.7 : #not even 70 percent alignment
#            lca=True
        if mismatch > 10: #hard rule, 10 mismatches
            lca=True
        if mismatch/ curlen > 0.1: #10 percent of read is mismatching
            lca=True
        return lca
    def is_low_phread(self, minp=20):
        return self.meanphread_ <= minp
