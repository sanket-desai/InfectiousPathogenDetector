#!/usr/bin/python3
'''
Author      : Sanket Desai
Date        : 27/08/2020
Description : Create IPD Database
'''
import sys
import os
import datetime
from Bio import SeqIO
from Bio import Entrez
from globals import *
import subprocess
import os
import argparse
import time
import gzip
'''
create fasta given ncbi IDS
'''
class IPDUpdateException(Exception):
    def __init__(self, exstr=""):
        Exception.__init__(self, exstr)

class IPDDatabaseCreator(object):
    def __init__(self):
        #Format type good and file exists
        #GlobalVar.initialize()
        self.failedgb_=[] #gbk ids which failed
        self.downloadcount_=0
        self.pathogbkfname_=GlobalVar.pathoids_+".gbk"
        self.secondarygbkfname_=GlobalVar.secondarydbids_+".gbk"
        self.primarypathoids_=[]
        pathoidfh=open(GlobalVar.pathoids_)
        processedsecdbids=[]
        for p in pathoidfh:
            if not p.startswith("#"):
                self.primarypathoids_.append(p.strip())
        pathoidfh.close()
        #self.secondarydbids_=[]
        #secdbidfh=open(GlobalVar.secondarydbids_)
        #for s in secdbidfh:
        #    if not s.startswith("#"):
        #        self.secondarydbids_.append(s.strip())
        #secdbidfh.close()
        print("Read %d primary IDS.." %(len(self.primarypathoids_)) )
        #Open up the files for writing
        pathofafi=open(GlobalVar.pathofa_,'w')
        secdbfafi=open(GlobalVar.blastndb_,'w')
        pathogbkfi=open(self.pathogbkfname_,'a')
        #secdbgbkfi=open(self.secondarygbkfname_,'a')
        pathogfffi=open(GlobalVar.pathogff_,'a')
        pathoannfi=open("../data/annotation/primaryannotation.tsv",'a')
        secdbannfi=open("../data/annotation/secondaryannotation.tsv",'w')
        self.idlist_to_localgb(self.primarypathoids_, pathogbkfi)
        while len(self.failedgb_)>0:
            self.idlist_to_localgb(self.failedgb_, pathogbkfi)
        self.failedgb_=[]
        pathogbkfi.close()
        pathogiter=SeqIO.parse(self.pathogbkfname_, "genbank")
        #Write pathogff, annotation and fasta for primary pathogen
        for pgenrec in pathogiter:
            anno=self.gbtoannotationtsvrecord(pgenrec)
            pathogfffi.write("%s\n" %(self.gbtogffrecord(pgenrec)) )
            pathoannfi.write("%s\n" %(anno))
            secdbannfi.write("%s\n" %(anno))
            pgenrec.description=""
            SeqIO.write(pgenrec, pathofafi,"fasta")
            SeqIO.write(pgenrec, secdbfafi, "fasta")
            processedsecdbids.append(pgenrec.name)
        pathogfffi.close()
        pathoannfi.close()
        pathofafi.close()
        print("Primary Database printed into IPD..!")
        try:
            cmd="ncbi-genome-download -s refseq --flat-output -l complete -F genbank -r 10 -v -R all -o ../data/annotation/gbktemp archaea"
            cprocess=subprocess.run(cmd, shell=True)
            cprocess.check_returncode()
            print("Completed genome download for Archaea.")
        except:
            print("ncbi-genome-download failed. Please re-install the ncbi-genome-download and re-run IPD installation!!")
            sys.exit(0)
        try:
            cmd="ncbi-genome-download -s refseq --flat-output -l complete -F genbank -r 10 -v -R all -o ../data/annotation/gbktemp viral"
            cprocess=subprocess.run(cmd, shell=True)
            cprocess.check_returncode()
            print("Completed genome download for Viral genomes.")
        except:
            print("ncbi-genome-download failed. Please re-install the ncbi-genome-download and re-run IPD installation!!")
            sys.exit(0)
        try:
            cmd="ncbi-genome-download -s refseq --flat-output -l complete -F genbank -r 10 -v -R all -o ../data/annotation/gbktemp bacteria"
            cprocess=subprocess.run(cmd, shell=True)
            cprocess.check_returncode()
            print("Completed genome download for Bacteria.")
        except:
            print("ncbi-genome-download failed. Please re-install the ncbi-genome-download and re-run IPD installation!!")
            sys.exit(0)
        #write to secondary database annotation and fasta
        #Combine with primary and secondary
        #cmd="zcat ../data/annotation/gbktemp/*.gz >> ../data/annotation/secondarydb.ids.gbk"
        #cprocess=subprocess.run(cmd, shell=True)
        #cprocess.check_returncode()
        gbktemp="../data/annotation/gbktemp/"
        secdbgbkout=open("../data/annotation/secondarydb.ids.gbk",'w')
        gind=0
        for g in os.listdir(gbktemp):
            try:
                gfi=gzip.open(gbktemp+g)
                for i in gfi:
                    secdbgbkout.write(i.decode('utf-8'))
                gfi.close()
                print("%d %s" %(gind+1, g))
                gind+=1
            except Exception as e:
                print("Could not open the gbk file : %s" %(g) )
                continue
        #pathogbkfi=open(self.pathogbkfname_)
        #for p in pathogbkfi:
        #    secdbgbkout.write(p)
        #    #this may need a \n
        pathogbkfi.close()
        secdbgbkout.close()
        try:
            secdbgiter=SeqIO.parse(self.secondarygbkfname_,"genbank")
            print("Secondary database creation started.")
            #processed secondary ids
            try:
                for sgenrec in secdbgiter:
                    if not sgenrec.name in processedsecdbids:
                        secdbannfi.write("%s\n" %(self.gbtoannotationtsvrecord(sgenrec)))
                        sgenrec.description=""
                        SeqIO.write(sgenrec, secdbfafi,"fasta")
                        processedsecdbids.append(sgenrec.name)
            except ValueError as v:
                print(v)
                pass
            print("Secondary database annotations created.")
            secdbannfi.close()
        except Exception as e:
            print(e)
            sys.exit(0)
        secdbfafi.close()

    def idlist_to_localgb(self, l, fo): #Takes in the list of genbank ids
        Entrez.email="abc@xyz.com"
        max_tries=5 #try downloading 5 times
        tries=0
        downloadfailedgb=[]
        for gid in l:
            downloadsuccessfull=False
            while tries < max_tries:
                try:
                    handle = Entrez.efetch(db = "nucleotide", id =gid, rettype = "gbwithparts", retmode="full")
                    fo.write(handle.read())
                    downloadsuccessfull=True
                    print("%d %s" %(self.downloadcount_+1, gid))
                    self.downloadcount_+=1
                    break
                except:
                    tries += 1
                    time.sleep(2)
                    continue
            if not downloadsuccessfull:
                downloadfailedgb.append(gid)
        self.failedgb_=downloadfailedgb
    def gbtoannotationtsvrecord(self, rec):
        # ID, genus, name, length, description, taxonomy
        annotations=[]
        try:
            annotations.append(rec.id)
            annotations.append(rec.annotations['taxonomy'][-1])
            annotations.append(rec.description)
            annotations.append(str(len(rec)))
            annotations.append(rec.description)
            annotations.append(",".join(rec.annotations['taxonomy']))
        except:
            print("Record ommitted, as the annotations are incomplete for the record : %s" %(rec.name))
        return "\t".join(annotations)
    #grec is the genbank record as generated by next function of seqio
    def gbtogffrecord(self, grec):
        gffrec=[]
        gffrec.append(grec.id)
        gffrec.append("GenBank")
        gffrec.append("gene")
        gffrec.append("1")
        gffrec.append(str(len(grec)))
        gffrec.append(".")
        gffrec.append("+")
        gffrec.append("1")
        attributes=[]
        attributes.append("ID="+grec.name)
        attributes.append("Dbxref="+",".join(grec.dbxrefs))
        attributes.append("Name="+grec.name)
        attributes.append("gene_id="+grec.description)
        gffrec.append(";".join(attributes))
        return "\t".join(gffrec)

def main():
    #Since all the input and output file names are mapped, just create
    i=IPDDatabaseCreator()
    print("Primary and Secondary database successfully created..")
if __name__=="__main__":
    main()
