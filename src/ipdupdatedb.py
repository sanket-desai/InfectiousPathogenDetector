#!/usr/bin/python3
'''
Author      : Sanket Desai
Date        : 18/08/2020
Description : Update the IPD database, index sequences and annotation
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
'''
Algorithm
Read the genbank file/s
Check if ID exists in the database
Create a fasta sequence file with >ID as the header and append to the current Pathogen sequence (1060) and human pathogen sequence fasta
Create GFF records and append to the human pathogen and pathogen annotation file (GFF)
Append the IDs to the existing IDs file
Append annotations to the primary and secondary annotation TSV files in data directory
Hisat indexing hspathofa
minimap2 indexing pathofa
samtools indexing of pathofa
Create annotation for snpEff
'''
import shutil
# remove_files_from_path source: https://stackoverflow.com/questions/185936/how-to-delete-the-contents-of-a-folder
def remove_files_from_path(fpath):
    folder = fpath
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

class IPDUpdateException(Exception):
    def __init__(self, exstr=""):
        Exception.__init__(self, exstr)

class IPDUpdate(object):
    def __init__(self, inmap): #Takes in the genbank file format file (multi gbk is also acceptable)
        if not inmap['format'].strip()=="gb" and not inmap['format'].strip()=="txt":
            raise IPDUpdateException("Please specify the right format type (acceptable inputs are \"gb\" or \"txt\")")
            sys.exit(0)
        if not os.path.isfile( inmap['input'].strip()):
            raise IPDUpdateException("File %s does not exist or is unreadable!!" %(inmap['input'].strip()) )
            sys.exit(0)
        self.gbkfname=inmap['input']
        #Format type good and file exists
        GlobalVar.initialize()
        #print(gbkids)
        self.failedgb_=[]
        self.downloadcount_=0
        self.primaryneedsindexing_=False
        self.secondaryneedsindexing_=False
        #load primary pathogen ids for reading
        self.primarypathoids_=[]
        pathoidfh=open(GlobalVar.pathoids_)
        for p in pathoidfh:
            if not p.startswith("#"):
                if "." in p:
                    p=p[:p.find(".")]
                self.primarypathoids_.append(p.strip())
        pathoidfh.close()
        self.secondarydbids_=[]
        secdbidfh=open(GlobalVar.secondarydbids_)
        for s in secdbidfh:
            if not s.startswith("#"):
                if "." in s:
                    s=s[:s.find(".")]
                self.secondarydbids_.append(s.strip())
        secdbidfh.close()
        gbtodownload=[]
        if inmap['format'].strip()=="txt":
            self.gbkfname=inmap['input']+".gbk"
            infi=open(inmap['input'])
            gbkids=[]
            for i in infi:
                if "." in i:
                    i=i[:i.find(".")]
                gbkids.append(i.strip())
            print("Read %d number of Genkbank IDS for the input file." %(len(gbkids)) )
            infi.close()
            for gi in gbkids:
                if not self.exists_in_primary(gi):
                    self.primaryneedsindexing_=True
                    gbtodownload.append(gi)
                if not self.exists_in_secondary(gi):
                    self.secondaryneedsindexing_=True
        elif inmap['format'].strip()=="gb":
            tempiter=SeqIO.parse(self.gbkfname, "genbank")
            for t in tempiter:
                tid=t.id
                if "." in tid:
                    tid=tid[:tid.find(".")]
                if not self.exists_in_primary(tid):
                    self.primaryneedsindexing_=True
                    gbtodownload.append(tid)
                if not self.exists_in_secondary(tid):
                    self.secondaryneedsindexing_=True
        if self.primaryneedsindexing_:
            #load fa for appending
            pathofafi=open(GlobalVar.pathofa_,'a')
            hspathofafi=open(GlobalVar.hspathofa_,'a')
            secdbfafi=open(GlobalVar.blastndb_,'a')
            #load gff for appending
            pathogfffi=open(GlobalVar.pathogff_,'a')
            hspathogfffi=open(GlobalVar.hspathogff_,'a')
            snpeffgbkfi=open(GlobalVar.snpeffgbk_,'a')
            #load annotation tsv for appending
            primannfi=open(GlobalVar.primaryannotationtsv_,'a')
            secannfi=open(GlobalVar.secondaryannotationtsv_,'a')
            #load primary pathogen ids for appending
            pathoidfh=open(GlobalVar.pathoids_,'a')
            secdbidfh=open(GlobalVar.secondarydbids_,'a')
            outf=open(self.gbkfname,'w')
            self.idlist_to_localgb(gbtodownload, outf)
            outf.close()
            giter=SeqIO.parse(self.gbkfname, "genbank")
            for genrec in giter:
                pathoidfh.write("#Record added at : %s\n" %(str(datetime.datetime.now())))
                #if "." in genrec.id:
                #    genrec.id=genrec.id[:genrec.id.find(".")]
                ss=genrec.id
                if "." in genrec.id:
                    ss=ss[:ss.find(".")]
                pathoidfh.write( "%s\n" %(ss))
                #Add in GFF
                pathogfffi.write("\n%s" %(self.gbtogffrecord(genrec)))
                hspathogfffi.write("\n%s" %(self.gbtogffrecord(genrec)))
                #Add in annotation tsv
                primannfi.write("%s\n" %(self.gbtoannotationtsvrecord(genrec)))
                secannfi.write("%s\n" %(self.gbtoannotationtsvrecord(genrec)))
                #snpEff annotation database building
                SeqIO.write(genrec, snpeffgbkfi,"genbank")
                #Add in fasta
                genrec.description=""
                SeqIO.write(genrec, hspathofafi,"fasta")
                SeqIO.write(genrec,pathofafi, "fasta")
                print( "Record added to IPD: %s\n" %(genrec.id))
            else:
                print("%s genome exists in IPD!" %(genrec.id))
                pass
            if self.secondaryneedsindexing_:
                secdbidfh.write("#Record added at : %s\n" %(str(datetime.datetime.now())))
                secdbidfh.write("%s\n" %(genrec.name))
                SeqIO.write(genrec, secdbfafi, 'fasta')
            else:
                print("%s genome exists in IPD!" %(genrec.name))
                pass
            print("Proceeding with indexing")
            pathofafi.close()
            hspathofafi.close()
            hspathogfffi.close()
            pathogfffi.close()
            pathoidfh.close()
            secdbidfh.close()
            secdbfafi.close()
            primannfi.close()
            secannfi.close()
            snpeffgbkfi.close()
            #Start indexing - if primary annotation, secondary annotation check
            #Index hspatho
            #hisat build
            cmd=GlobalVar.hisat2_build_+" "+GlobalVar.hspathofa_+" "+GlobalVar.hisat2hspathoindex_
            cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
            cprocess.check_returncode()
            #minimap2 index
            cmd=GlobalVar.minimap2_+" -d "+GlobalVar.minimap2index_+" "+GlobalVar.pathofa_
            cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
            cprocess.check_returncode()
            #samtools patho fa index
            cmd=GlobalVar.samtools_+" faidx "+GlobalVar.pathofa_
            cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
            cprocess.check_returncode()
            #snpEff database building
            cmd=GlobalVar.snpeff_+" build -genbank -v ipd1060"
            cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
            cprocess.check_returncode()
            if self.secondaryneedsindexing_:
                cmd=GlobalVar.makeblastdb_+" -dbtype nucl -in "+GlobalVar.blastndb_
                cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
                cprocess.check_returncode()

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
    def exists_in_primary(self, gbid): #Given a id  check if it exist in IPD
        return gbid in self.primarypathoids_
    def exists_in_secondary(self, gbid):
        return gbid in self.secondarydbids_
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
    aparser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description= '''Infectious Pathogen Detector Database Updater: Module to update the primary pathogen database.
Version v0.1.0''')
    aparser.add_argument('-f', action='store', dest='format', help='File format (expected input is either "gb" or "txt" ', required=True)
    aparser.add_argument('-i', action='store', dest='input', help='Input file to read (either genbank file or text file with Genbank IDS one per line)', required=True)
    pargs=aparser.parse_args()
    try:
        inmap={}
        inmap['format']=pargs.format
        inmap['input']=pargs.input
        i=IPDUpdate(inmap)
    except TypeError as t:
        aparser.print_help()
        print(t)
        print("\n")
        sys.exit(0)

if __name__=="__main__":
    main()
