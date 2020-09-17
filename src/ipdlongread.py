#!/usr/bin/python3
'''

Author      : Sanket Desai
Date        : 25/06/2020
Description : IPD Long Read main classes

'''
import sys
import os
import subprocess
import pysam
from globals import *
from seq import *
import ntpath
import multiprocessing
from readassessment import *
from localblastn import *
import argparse
import datetime
from functools import reduce
'''
Input for the short read module
prefix - prject name
paired/single end - librarytype
r1 - read 1
r2 - read 2
output directory - outdir
threads - threads
molecule - dna/rna
datatype - human/microbe

A map with following keys
prefix
r1
r2
outdir
threads
molecule
datatype
'''

class FeatureCountsRecord(object):
    def __init__(self, line, treads):
        self.treads_=treads
        sline=line.strip().split("\t")
        if len(sline) < 7:
            raise GlobalFormatException("Featurecount file format exception at : %s\n" %(line))
        self.Geneid_=sline[0].strip()
        self.Chr_=sline[1].strip()
        self.Start_= 0
        self.End_=1
        self.Length_=1
        try:
            self.Start_=int(sline[2].strip())
        except ValueError:
            self.Start_=int(sline[2].strip().split(";")[0])
        try:
            self.End_=int(sline[3].strip())
        except ValueError:
            self.End_=int(sline[3].strip().split(";")[0])
        self.Strand_=sline[4].strip()
        self.Length_=int(sline[5].strip())
        self.count_=float(sline[6].strip())
    def get_count(self):
        return self.count_
    def get_length(self):
        return self.Length_
    def get_feature(self):
        return self.Geneid_
    def get_id(self):
        return self.Chr_
    def get_fpm(self):
        return (self.count_ / self.treads_ ) * 1000000
    def get_fpkm(self):
        return ( (self.count_ / (self.Length_ / 1000) ) / self.treads_ ) *1000000

class IPDLongRead(object):
    def __init__(self, inputmap ):
        self.inputr1_isgzip_=False #[could be gzip or text ]
        #self.inputr2_isgzip_=False
        self.inputr1_extn_=""
        #self.inputr2_extn_=""
        self.inputmap_=inputmap
        #self.data_ispaired_=False
        #if "r2" in self.inputmap_:
        #    self.data_ispaired_=True
        #first sam
        r1basename=ntpath.basename(self.inputmap_['r1'])
        #r2basename=""
        self.r1fnameraw_=self.inputmap_['r1']
        #self.r2fnameraw_=""
        #if self.data_ispaired_:
        #    self.r2fnameraw_=self.inputmap_['r2']
        #    r2basename=ntpath.basename(self.inputmap_['r2'])
        if not self.inputmap_["outdir"].endswith('/'):
            self.inputmap_["outdir"]=self.inputmap_["outdir"]+'/'
        ty1=subprocess.run( "file -b "+self.r1fnameraw_ , shell=True, stdout=subprocess.PIPE)
        if "gzip" in ty1.stdout.decode('utf-8').strip():
            self.inputr1_isgzip_=True
        if self.inputr1_isgzip_:
            if self.r1fnameraw_.endswith(".fq.gz"):
                self.inputr1_extn_=".fq.gz"
                #self.inputmap_['r1']= self.inputmap_['outdir']+r1basename+"_clean.fq.gz"
            elif self.r1fnameraw_.endswith(".fastq.gz"):
                self.inputr1_extn_=".fastq.gz"
                #self.inputmap_['r1']= self.inputmap_['outdir']+r1basename+"_clean.fastq.gz"
            else:
                raise InputFileTypeException("Unusual input read 1 file input type: %s \nAccepted file extensions are 'fq.gz', 'fastq.gz'!!" %(self.inputr1_extn_))
        else:
            if self.r1fnameraw_.endswith(".fq"):
                self.inputr1_extn_=".fq"
                #self.inputmap_['r1']= self.inputmap_['outdir']+r1basename+"_clean.fq"
            elif self.r1fnameraw_.endswith(".fastq"):
                self.inputr1_extn_=".fastq"
                #self.inputmap_['r1']= self.inputmap_['outdir']+r1basename+"_clean.fastq"
            else:
                raise InputFileTypeException("Unusual input read 1 file input type: %s \nAccepted file extensions are 'fq', 'fastq'!!" %(self.inputr1_extn_))
        r1basename=r1basename.replace(self.inputr1_extn_,"")
        #Initialization
        GlobalVar.initialize()
        print("Data loaded!! Starting data filtration..")
        #Filtration with prechop and NanoFilt
        self.primarysortedbam_=""
        self.filteredfastq_=""
        self.samplebasename_=r1basename
        self.filtration()
        print("Filtration completed")
        self.primarymetadata_=PrimaryMetadata(GlobalVar.primaryannotationtsv_)
        #Post alignment processing
        self.genomecoveragemap_={} #Need to change the strategy of genome coverage calculation
        #pathoid_len_map=self.primarymetadata_.get_pathogenid_length_map()
        #print(pathoid_len_map)
        #for p in pathoid_len_map:
        #    self.genomecoveragemap_[p]=[0]*pathoid_len_map[p]
        #print("Genome cov map created -- " + str(len(self.genomecoveragemap_)))
        #sys.exit(0)
        #self.primarycountsam_=self.primarysam_.replace(".sam", "_count.bam")
        self.primarycountbam_=self.inputmap_["outdir"]+self.inputmap_["prefix"]+"_"+r1basename+"_count.bam"
        #file to keep only the filtered pathogen sam records
        #self.assembledcontigdirectory_=self.primarysam_.replace(".sam","_assembledcontigs")
        #self.assembledcontigsfinalfile_=self.assembledcontigdirectory_+"/final.contigs.fa"
        #self.assembledcontigsblastoutput_=self.primarysam_.replace(".sam", "_assembledcontigs_blastout.tsv")
        #minimap2 alignment
        self.primary_alignment()
        print("Primary alignment completed")
        self.featurecountsoutputfile_=""
        self.ipdfinalcountsfile_=""
        self.ipdfinalvcf_=self.primarycountbam_.replace(".bam", "_final.vcf")
        self.ipdfinalannotatedvcf_=self.primarycountbam_.replace(".bam", "_final_annotated.vcf")
        self.total_reads_=0
        sfile=pysam.AlignmentFile(self.primarycountbam_,'rb')
        for read in sfile.fetch(until_eof=True):
            self.total_reads_=self.total_reads_+1
        sfile.close()
        #self.totalreads_=reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(self.primarycountbam_) ])

        #Code can fork here -- creating two different processes (multiprocess) / one for counting and other for variant call
        #align the pathofasta file with secondary reference database and perform filtration
        print("Splitting processes")
        cproc=multiprocessing.Process(target=self.counting_process)
        vproc=multiprocessing.Process(target=self.variant_process)
        cproc.start()
        vproc.start()
        cproc.join()
        vproc.join()

    def variant_process(self):
        #print("variant process started..")
        #currenttime=datetime.datetime.now()
        haltvarproc=False
        try:
            self.variantpreprocessing()
        except subprocess.CalledProcessError:
            print("Variant preprocessing failed. PICARD error!!")
            haltvarproc=True
            #sys.exit(0)
        if not haltvarproc:
            try:
                self.variantcalling()
            except subprocess.CalledProcessError:
                print("Variant preprocessing failed. Error in the variant calling tools / annotation !!")
                pass

    def variantpreprocessing(self):
        #sam format conversion
        #sort sam
        bam=self.primarycountbam_
        sortedbam=bam.replace(".bam","_sorted.bam")
        cmd=GlobalVar.picard_+" SortSam INPUT="+ bam +" OUTPUT="+sortedbam+" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        #bam indexing
        cmd=GlobalVar.picard_+" BuildBamIndex INPUT="+ sortedbam +" OUTPUT="+sortedbam+".bai VALIDATION_STRINGENCY=SILENT  TMP_DIR="+self.inputmap_['outdir']+"tmp"
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        self.primarysortedbam_=sortedbam

    def variantcalling(self):
        #medaka
        #cmd="conda init bash"
        #cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        #cprocess.check_returncode()
        #cmd='bash -c \"source activate medaka; medaka_variant -i "+self.primarysortedbam_ +" -f "+GlobalVar.pathofa_+" -p -t "+ str(self.inputmap_['threads'])+" -o "+ self.inputmap_['outdir']+self.samplebasename_+"_medaka_variant" #self.primarysortedbam_.replace(".bam", "_medaka.vcf\"')
        #medakacmd="medaka_variant -i "+self.primarysortedbam_ +" -f "+GlobalVar.pathofa_+" -p -t "+ str(self.inputmap_["threads"])+" -o "+ self.inputmap_["outdir"]+self.samplebasename_+"_medaka_variant"
        #cmd="source activate medaka && 'medaka_variant -i "+self.primarysortedbam_ +" -f "+GlobalVar.pathofa_+" -p -t "+ str(self.inputmap_["threads"])+" -o "+ self.inputmap_["outdir"]+self.samplebasename_+"_medaka_variant"' && source deactivate"
        #cmd="source activate medaka ; "+medakacmd
        #my_env = os.environ.copy()
        #cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)#, env=my_env)
        #cprocess=subprocess.run(medakacmd, shell=True, stdout=subprocess.DEVNULL)#, env=my_env)
        #cprocess.check_returncode()
        #cmd="conda activate medaka"
        #cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        #cprocess.check_returncode()
        #cmd="medaka_variant -i "+self.primarysortedbam_ +" -f "+GlobalVar.pathofa_+" -p -t "+ str(self.inputmap_['threads'])+" -o "+ self.inputmap_['outdir']+self.samplebasename_+"_medaka_variant" #self.primarysortedbam_.replace(".bam", "_medaka.vcf")
        cmd="bash -i conda_activate.sh "+self.primarysortedbam_ +" "+GlobalVar.pathofa_+" "+self.inputmap_['outdir']+self.samplebasename_+"_medaka_variant"+" "+str(self.inputmap_['threads']) #self.primarysortedbam_.replace(".bam", "_medaka.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.vcffilter_+ " -f \"QUAL > 20\" " + self.inputmap_['outdir']+self.samplebasename_+"_medaka_variant/round_1.vcf > "+self.primarysortedbam_.replace(".bam","_medaka.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.samtools_+" mpileup -q 1 -f "+GlobalVar.pathofa_+" -o "+self.primarysortedbam_.replace(".bam",".mpileup ")+self.primarysortedbam_
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.varscan_+" mpileup2cns "+self.primarysortedbam_.replace(".bam",".mpileup")+" --output-vcf 1 --variants > "+self.primarysortedbam_.replace(".bam","_varscan.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.lofreq_+" call -f "+GlobalVar.pathofa_+" -o "+self.primarysortedbam_.replace(".bam","_lofreq1.vcf")+" "+self.primarysortedbam_
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.vcffilter_+" -f \"DP > 10 & AF > 0.75\" "+ self.primarysortedbam_.replace(".bam","_lofreq1.vcf") +" > "+self.primarysortedbam_.replace(".bam","_lofreq.vcf")
        #print(cmd)
        cprocess=subprocess.run(cmd, shell=True)
        cprocess.check_returncode()
        #freebayes two-pass variant calling
        cmd=GlobalVar.freebayes_+" -f "+GlobalVar.pathofa_+"  "+self.primarysortedbam_+" | "+ GlobalVar.vcffilter_+" -f \"QUAL > 20\" > "  +self.primarysortedbam_.replace(".bam","_freebayes_r1.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bgzip_+" "+self.primarysortedbam_.replace(".bam","_freebayes_r1.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bcftools_+" index "+self.primarysortedbam_.replace(".bam","_freebayes_r1.vcf.gz")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()

        cmd=GlobalVar.freebayes_+" -f "+GlobalVar.pathofa_+"  --min-repeat-entropy 1 --haplotype-length 500 --haplotype-basis-alleles "+ self.primarysortedbam_.replace(".bam","_freebayes_r1.vcf.gz")+" "+self.primarysortedbam_+" | "+ GlobalVar.vcffilter_+" -f \"QUAL > 20 & DP > 10 & AF > 0.75\" > "  +self.primarysortedbam_.replace(".bam","_freebayes.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()

        cmd=GlobalVar.bgzip_+" "+self.primarysortedbam_.replace(".bam","_medaka.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bgzip_+" "+self.primarysortedbam_.replace(".bam","_freebayes.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bgzip_+" "+self.primarysortedbam_.replace(".bam","_lofreq.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bgzip_+" "+self.primarysortedbam_.replace(".bam","_varscan.vcf")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bcftools_+" index "+self.primarysortedbam_.replace(".bam","_medaka.vcf.gz")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bcftools_+" index "+self.primarysortedbam_.replace(".bam","_lofreq.vcf.gz")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bcftools_+" index "+self.primarysortedbam_.replace(".bam","_varscan.vcf.gz")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        cmd=GlobalVar.bcftools_+" index "+ self.primarysortedbam_.replace(".bam","_freebayes.vcf.gz")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        #consensus call - sr.vcf
        cmd=GlobalVar.bcftools_+" isec -n +2 -o "+self.primarysortedbam_.replace(".bam","_sr.vcf.gz")+" "+ self.primarysortedbam_.replace(".bam","_lofreq.vcf.gz ")+ self.primarysortedbam_.replace(".bam","_varscan.vcf.gz ")+self.primarysortedbam_.replace(".bam","_freebayes.vcf.gz -f PASS -O z -w 2")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        #Filter with DP 10 and
        #cmd=GlobalVar.bgzip_+" "+ self.primarysortedbam_.replace(".bam","_sr.vcf")
        #cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        #cprocess.check_returncode()
        cmd=GlobalVar.bcftools_+" index "+ self.primarysortedbam_.replace(".bam","_sr.vcf.gz")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()

        #common between consensus and medaka
        cmd=GlobalVar.bcftools_+" isec -n +2 -o "+self.ipdfinalvcf_+" "+self.primarysortedbam_.replace(".bam","_sr.vcf.gz ")+ self.primarysortedbam_.replace(".bam","_medaka.vcf.gz -f PASS -O v -w 2")
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()
        print("Final variant calls made")
        #Where it fails
        #cmd=GlobalVar.vcffilter_+" -f \"DP > 10 & AF > 0.75\" "+ self.ipdfinalvcf_.replace(".vcf","n1.vcf") +" > "+self.ipdfinalvcf_
        #print(cmd)
        #cprocess=subprocess.run(cmd, shell=True)
        #cprocess.check_returncode()

#        cmd=GlobalVar.bcftools_+" view -v indels "+self.primarysortedbam_.replace(".bam","_sr.vcf.gz -o ")+ self.primarysortedbam_.replace(".bam", "_sr_indel.vcf.gz -O z")
#        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
#        cprocess.check_returncode()
#        cmd=GlobalVar.bcftools_+" index "+self.primarysortedbam_.replace(".bam","_snp.vcf.gz")
#        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
#        cprocess.check_returncode()
#        cmd=GlobalVar.bcftools_+" index "+self.primarysortedbam_.replace(".bam","_sr_indel.vcf.gz")
#        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
#        cprocess.check_returncode()

#        cmd=GlobalVar.bcftools_+" isec -n 1 -o "+self.ipdfinalvcf_+ " "+self.primarysortedbam_.replace(".bam","_snp.vcf.gz ")+ self.primarysortedbam_.replace(".bam","_sr_indel.vcf.gz -f PASS -O v -w 1")
#        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
#        cprocess.check_returncode()
        cmd=GlobalVar.snpeff_+" -dataDir "+GlobalVar.snpeffdatadir_+" -nodownload ipd1060 "+self.ipdfinalvcf_+" > "+self.ipdfinalannotatedvcf_
        #Commented only to be run on Param / otherwise uncomment this and RUN
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        cprocess.check_returncode()

    def counting_process(self):
        #self.primarysamfilehandle_.close()
        #self.assemblyfastqinput1filehandle_.close()
        self.featurecountsoutputfile_=self.primarycountbam_.replace(".bam", "_featurecounts.tsv")
        self.ipdfinalcountsfile_=self.primarycountbam_.replace(".bam", "_finalcounts.tsv")
        self.counting()

    def filtration(self):
        chopped=self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_chopped"+self.inputr1_extn_
        cmd=GlobalVar.porechop_+" -i "+ self.inputmap_['r1'] +" -o "+ chopped
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        self.filteredfastq_=chopped.replace("_chopped"+self.inputr1_extn_,"_filtered"+self.inputr1_extn_)
        if not self.filteredfastq_.endswith(".gz"):
            self.filteredfastq_=self.filteredfastq_+".gz"
        try:
            cprocess.check_returncode()
        except subprocess.CalledProcessError:
            print("Filtration failed !!")
            sys.exit(0)
        if self.inputr1_isgzip_:
            cmd="gunzip -c "+chopped+ " | NanoFilt -q 13 -l 500 | gzip > "+self.filteredfastq_
            cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
            cprocess.check_returncode()
        else:
            #cmd="NanoFilt -q 13 -l 500 "+chopped+"| gzip > "+chopped.replace("_chopped"+self.inputr1_extn_,"_filtered"+self.inputr1_extn_+".gz")
            cmd="NanoFilt -q 13 -l 500 "+chopped+"| gzip > "+self.filteredfastq_
            cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
            cprocess.check_returncode()
    #argument: complete output path for sam file
    def primary_alignment(self):
        #cmd=GlobalVar.minimap2_+" -ax map-ont -k 13 -t "+str(self.inputmap_["threads"])+" --splice --MD -Y "+GlobalVar.minimap2index_+ " "+self.filteredfastq_+ " | "+GlobalVar.samtools_+" view -bh > "+self.primarycountbam_
        cmd=GlobalVar.minimap2_+" -ax map-pb -k 13 -t "+str(self.inputmap_["threads"])+" --splice --MD -Y "+GlobalVar.minimap2index_+ " "+self.filteredfastq_+ " | "+GlobalVar.samtools_+" view -bh > "+self.primarycountbam_
        cprocess=subprocess.run(cmd, shell=True)
        #raises subprocess.CalledProcessError if returncode is non-zero
        try:
            cprocess.check_returncode()
        except Exception:
            print("Primary alignment failed (cmd): %s" %(cmd))
            sys.exit(0)

    def counting(self):
        #construct cmd for featurecounts and execute using run
        cmd=GlobalVar.featurecounts_+" -a "+GlobalVar.pathogff_+" -t gene -g gene_id -M -O -T " + str(self.inputmap_['threads']) +" --fraction -o "+self.featurecountsoutputfile_+"  "+self.primarycountbam_
        cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
        #raises subprocess.CalledProcessError if returncode is non-zero
        cprocess.check_returncode()
        #Need to postprocess the counts created
        featcountfi=open(self.featurecountsoutputfile_)
        ipdfinalcountfo=open(self.ipdfinalcountsfile_,'w')
        if not self.total_reads_>0:
            raise GlobalFormatException("Total reads found to be '0' in the input file!! Aborting...")
        ipdfinalcountfo.write("%s\t%s\t%s\t%s\t%s\t%s\n" %("Feature","ID","Length","Fragments", "FPM", "FPKM") )
        for f in featcountfi:
            if not f.startswith("#") and not f.startswith("Geneid"):
                frec=FeatureCountsRecord(f, self.total_reads_)
                ipdfinalcountfo.write("%s\t%s\t%d\t%f\t%f\t%f\n" %(frec.get_feature(), frec.get_id(), frec.get_length(), frec.get_count(), frec.get_fpm(), frec.get_fpkm() ) )
        featcountfi.close()
        ipdfinalcountfo.close()

def main():
    aparser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description= '''Infectious Pathogen Detector (Long read): Quantification and variant analysis of pathongens from long-read nanopore sequencing data
Developed by Dutt lab
Version v0.2.0''')
    #aparser.add_argument('-m', action='store', dest='mode', help='Multisample or single sample mode, options [multi or single]')
    aparser.add_argument('-p', action='store', dest='prefix', help='Project name')
    aparser.add_argument('-1', action='store', dest='fastq', help='FASTQ file (with complete path)')
    aparser.add_argument('-o', action='store', dest='outdir', help='Output directory path')
    aparser.add_argument('-t', action='store', dest='threads', required=False, default="4" , help='Number of threads (default = 4)')
    aparser.add_argument('-s', action='store', dest='species', required=False, default="microbe", help='Species sequenced, option can be human / microbe (default = microbe)')
    pargs=aparser.parse_args()
    try:
        inmap={}
        inmap['prefix']=pargs.prefix
        inmap['r1']=pargs.fastq
        inmap['outdir']=pargs.outdir
        inmap['threads']=pargs.threads
        inmap['datatype']=pargs.species
        i=IPDLongRead(inmap)

    except TypeError:
        aparser.print_help()
        print("\n")
        sys.exit(0)

if __name__=="__main__":
    main()
'''
prefix - prject name
paired/single end - librarytype
r1 - read 1
r2 - read 2
output directory - outdir
threads - threads
molecule - dna/rna
datatype - human/microbe
'''
