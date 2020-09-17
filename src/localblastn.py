import sys
import os
import subprocess
from subprocess import PIPE

#perform blastN

class LocalBlastN():
    def __init__(self, blastnpath, blastdb, inputfastafile, numalignments, threads, outputtsvfile, outputformat="6"):
        self.blastnpath_=blastnpath
        self.blastndb_=blastdb
        self.inputfasta_=inputfastafile
        self.numalignments_=numalignments
        self.outputfile_=outputtsvfile
        self.threads_=threads
        self.outputformat_=outputformat
        cmd=self.blastnpath_+" -db "+self.blastndb_ + " -query " + self.inputfasta_ + " -outfmt "+ self.outputformat_+" -num_alignments "+str(numalignments)+" -num_threads "+str(self.threads_) + " -out " + self.outputfile_
        cprocess=subprocess.run(cmd, shell=True)
        cprocess.check_returncode()
        self.blastoutfilehandle_=open(self.outputfile_, 'r')

        #self.blastoutputheader_=self.blastoutfilehandle_.readline().strip().split("\t")
    def get_blastoutputformat6record(self):
        l=self.blastoutfilehandle_.readline()
        return BlastOutputFormat6Record(self.blastoutfilehandle_.readline())
    def close_connection(self):
        self.blastoutfilehandle_.close()
    def executeBlastnOutfmt6Wrapper(self, sequence , numofalignments):
        blastcmd="echo -e \">"+sequence.get_name()+"\\n"+sequence.get_sequence().upper()+"\" | "+self.blastnpath_+" -db "+self.blastndb_+" -outfmt 6 -num_alignments "+str(numofalignments)
        popenobj=subprocess.Popen(blastcmd, shell=True, stdout=PIPE, stderr=PIPE)
        blastnformat6recs=[]
        for bl in popenobj.stdout:
            blastnformat6recs.append(BlastOutputFormat6Record(bl.decode("utf-8")))
        return blastnformat6recs
#Single record for blast output format 6 - a 12 column record
class BlastOutputFormat6Record():
    def __init__(self, rec):
        if rec == '':
            raise ValueError
        self.srec_=rec.strip().split("\t")
        self.qseqid_=self.srec_[0]
        self.sseqid_=self.srec_[1]
        self.pident_=float(self.srec_[2])
        self.length_=int(self.srec_[3])
        self.mismatch_=int(self.srec_[4])
        self.gapopen_=int(self.srec_[5])
        self.qstart_=int(self.srec_[6])
        self.qend_=int(self.srec_[7])
        self.sstart_=int(self.srec_[8])
        self.send_=int(self.srec_[9])
        self.evalue_=float(self.srec_[10])
        self.bitscore_=float(self.srec_[11])
    def __str__(self):
        self.to_string()
    def to_string(self):
        return "\t".join(self.srec_)
    def is_identical(self, thresh=99.99):
        return self.pident_ >= thresh
    def is_significant(self, thresh=0.0001):
        return self.evalue_ < thresh
    def as_list(self):
        return self.srec_
