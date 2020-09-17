#!/usr/bin/python3
'''

Author	  : Sanket Desai
Date		: 25/04/2020
Description : IPD Short Read main classes

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

import datetime
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

class IPDFilterStats(object):
	qualityfiltered_=0 # average phread > 20
	badalignmentfiltered_=0 # 15/15 bases 3',5' soft clipping / rest read should be either M, I, D
	ambigousfiltered_=0 #fitration after secondary alignment
	complexityfiltered_=0 #AT 87% GC 89 % A,T,G,C 60%
	lengthfiltered_=0 # read length < 40
	totalreads_=0

class FeatureCountsRecord(object):
	def __init__(self, line):
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
		return (self.count_ / IPDFilterStats.totalreads_ ) * 1000000
	def get_fpkm(self):
		return ( (self.count_ / (self.Length_ / 1000) ) / IPDFilterStats.totalreads_ ) *1000000

class IPDShortRead(object):
	def __init__(self, inputmap ):
		self.inputr1_isgzip_=False #[could be gzip or text ]
		self.inputr2_isgzip_=False
		self.inputr1_extn_=""
		self.inputr2_extn_=""
		self.inputmap_=inputmap
		self.data_ispaired_=False
		if "r2" in self.inputmap_:
			self.data_ispaired_=True
		#first sam
		r1basename=ntpath.basename(self.inputmap_['r1'])
		r2basename=""
		self.r1fnameraw_=self.inputmap_['r1']
		self.r2fnameraw_=""
		if self.data_ispaired_:
			self.r2fnameraw_=self.inputmap_['r2']
			r2basename=ntpath.basename(self.inputmap_['r2'])
		if not self.inputmap_["outdir"].endswith('/'):
			self.inputmap_["outdir"]=self.inputmap_["outdir"]+'/'
		ty1=subprocess.run( "file -b "+self.r1fnameraw_ , shell=True, stdout=subprocess.PIPE)
		ty2=""
		if self.data_ispaired_:
			ty2=subprocess.run( "file -b "+self.r2fnameraw_ , shell=True, stdout=subprocess.PIPE)
			if "gzip" in ty2.stdout.decode('utf-8').strip():
				self.inputr2_isgzip_=True
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
		if self.data_ispaired_:
			if self.inputr2_isgzip_:
				if self.r2fnameraw_.endswith(".fq.gz"):
					self.inputr2_extn_=".fq.gz"
					#self.inputmap_['r2']= self.inputmap_['outdir']+r2basename+"_clean.fq.gz"
				elif self.r2fnameraw_.endswith(".fastq.gz"):
					self.inputr2_extn_=".fastq.gz"
					#self.inputmap_['r2']= self.inputmap_['outdir']+r2basename+"_clean.fastq.gz"
				else:
					raise InputFileTypeException("Unusual input read 1 file input type: %s \nAccepted file extensions are 'fq.gz', 'fastq.gz'!!" %(self.inputr2_extn_))
			else:
				if self.r2fnameraw_.endswith(".fq"):
					self.inputr2_extn_=".fq"
					#self.inputmap_['r2']= self.inputmap_['outdir']+r2basename+"_clean.fq"
				elif self.r2fnameraw_.endswith(".fastq"):
					self.inputr2_extn_=".fastq"
					#self.inputmap_['r2']= self.inputmap_['outdir']+r2basename+"_clean.fastq"
				else:
					raise InputFileTypeException("Unusual input read 1 file input type: %s \nAccepted file extensions are 'fq.gz', 'fastq.gz'!!" %(self.inputr2_extn_))
		#inputmap r1 and r2 have been changed to clean file files , for running the fastq use r1fnameraw_ and r2fnameraw_
		r1basename=r1basename.replace(self.inputr1_extn_,"")
		self.inputmap_['r1']= self.inputmap_['outdir']+r1basename+"_clean"+self.inputr1_extn_
		if self.data_ispaired_:
			r2basename=r2basename.replace(self.inputr2_extn_, "")
			self.inputmap_['r2']= self.inputmap_['outdir']+r2basename+"_clean"+self.inputr2_extn_
		#Initialization
		GlobalVar.initialize()

		self.samplebasename_=r1basename
		#self.primarysam_=self.inputmap_["outdir"]+self.inputmap_["prefix"]+"_"+r1fname[:r1fname.find(".")]+".sam"
		self.primarysam_=self.inputmap_["outdir"]+self.inputmap_["prefix"]+"_"+r1basename+".sam"
		self.fastphtml_=self.primarysam_.replace(".sam", "_fastp.html")
		self.fastpjson_=self.primarysam_.replace(".sam", "_fastp.json")

		#Filtration with fastp - Commented for Param test
		self.filtration()

		self.primarymetadata_=PrimaryMetadata(GlobalVar.primaryannotationtsv_)
		self.secondarymetadata_=SecondaryMetadata(GlobalVar.secondaryannotationtsv_)
		#Post alignment processing
		self.genomecoveragemap_={} #Need to change the strategy of genome coverage calculation
		pathoid_len_map=self.primarymetadata_.get_pathogenid_length_map()
		#print(pathoid_len_map)
		#for p in pathoid_len_map:
		#	self.genomecoveragemap_[p]=[0]*pathoid_len_map[p]
		#print("Genome cov map created -- " + str(len(self.genomecoveragemap_)))
		#sys.exit(0)
		self.primaryalignmentsummary_=self.primarysam_.replace(".sam", "_alignmentsummary.txt")
		#store the sam id -> reference name pair stored
		self.pathogenaligned_r1_samid_ref_={}
		self.pathogenaligned_r2_samid_ref_={}
		#samid -> not ambiguous (boolean), False means ambigous
		self.samidr1_nonambiguity_={}
		self.samidr2_nonambiguity_={}
		#perform primary alignment
		print("Data loaded, going for primary alignment!")
		#Param testing
		self.primary_alignment()
		if self.inputmap_['datatype']=='dna':
			self.dna_preprocessing()
			#dna_preprocessing adds _dna to the primary sam file name
		print("Primary alignment completed..")
		#sam file to be given for featurecounts - counting sam
		self.primarycountsam_=self.primarysam_.replace(".sam", "_count.bam")
		self.primaryfilteredpathogensam_=self.primarysam_.replace(".sam","_filteredpathogen.sam")
		#Paired pathogen sam reads to fasta file - input for hsblast / secondary alignment
		self.pathoalignedfasta_=self.primarysam_.replace(".sam", "_pathoseqs.fa")
		#fastq files created to be given input to the assembly tool
		self.secondaryalignedtsv_=self.pathoalignedfasta_.replace("_pathoseqs.fa", "_pathosecondaryaligned.tsv")
		self.assemblyfastqinput1_=self.primarysam_.replace(".sam", "_assemblyinput_R1.fq")
		self.assemblyfastqinput2_=self.primarysam_.replace(".sam", "_assemblyinput_R2.fq")
		try:
			self.assemblyfastqinput1filehandle_=open(self.assemblyfastqinput1_, 'w')
			if self.data_ispaired_:
				self.assemblyfastqinput2filehandle_=open(self.assemblyfastqinput2_, 'w')
		except FileNotFoundError:
			print(self.assemblyfastqinput1_)
			sys.exit(0)
		self.assembledcontigdirectory_=self.primarysam_.replace(".sam","_assembledcontigs")
		self.assembledcontigsfinalfile_=self.assembledcontigdirectory_+"/final.contigs.fa"
		self.assembledcontigsblastoutput_=self.primarysam_.replace(".sam", "_assembledcontigs_blastout.tsv")
		#primary alignment sam file / file reader
		self.primarysamfilehandle_=pysam.AlignmentFile(self.primarysam_,'r')
		self.primarycountsamfilehandle_=pysam.AlignmentFile(self.primarycountsam_, 'wb', template=self.primarysamfilehandle_)
		#fastq file(s) for assembly generated, fasta generated (for secondary alignment)
		#currenttime=datetime.datetime.now()
		self.postprimaryalignment_filtration()
		#timediff=datetime.datetime.now()-currenttime
		#print("Filtration completed.. Time taken : ")
		#print(timediff.total_seconds()/60)
		self.featurecountsoutputfile_=""
		self.ipdfinalcountsfile_=""
		self.preprocessedpathogenbam_=""
		self.ipdfinalvcf_=self.primarysam_.replace(".sam", "_final.vcf")
		self.ipdfinalannotatedvcf_=self.primarysam_.replace(".sam", "_final_annotated.vcf")
		#Code can fork here -- creating two different processes (multiprocess) / one for counting and other for variant call
		#align the pathofasta file with secondary reference database and perform filtration
		cproc=multiprocessing.Process(target=self.counting_process)
		vproc=multiprocessing.Process(target=self.variant_process)
		cproc.start()
		vproc.start()
		cproc.join()
		vproc.join()
	def dna_preprocessing(self):
		bam=self.primarysam_.replace(".sam","_primary.bam")
		#sam format conversion
		cmd=GlobalVar.picard_+" SamFormatConverter INPUT="+self.primarysam_+" OUTPUT="+bam+" VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		#raises subprocess.CalledProcessError if returncode is non-zero
		cprocess.check_returncode()
		#Fixmate
		fxdbam=bam.replace(".bam","_fxd.bam")
		cmd=GlobalVar.picard_+" FixMateInformation INPUT="+ bam +" OUTPUT="+fxdbam+" VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#sort sam
		sortedbam=bam.replace(".bam","_fxd_sorted.bam")
		cmd=GlobalVar.picard_+" SortSam INPUT="+ fxdbam +" OUTPUT="+sortedbam+" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#duplicate marking
		rmdupbam=bam.replace(".bam","_fxd_sorted_duprm.bam")
		cmd=GlobalVar.picard_+" MarkDuplicates INPUT="+ sortedbam +" OUTPUT="+rmdupbam+" METRICS_FILE="+rmdupbam.replace(".bam","_info.txt")+ " REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#bam indexing
		cmd=GlobalVar.picard_+" BuildBamIndex INPUT="+ rmdupbam +" OUTPUT="+rmdupbam+".bai VALIDATION_STRINGENCY=SILENT  TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#self.preprocessedpathogenbam_=rmdupbam
		self.primarysam_=self.primarysam_.replace(".sam", "_dna.sam")
		cmd=GlobalVar.samtools_+" view -h -@ "+str(self.inputmap_['threads'])+" -O SAM -o "+ self.primarysam_ + " " +rmdupbam
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()

	def variant_process(self):
		#print("variant process started..")
		#currenttime=datetime.datetime.now()
		try:
			self.variantpreprocessing()
		except subprocess.CalledProcessError:
			print("Variant preprocessing failed. PICARD error!!")
			sys.exit(0)

		try:
			self.variantcalling()
		except subprocess.CalledProcessError:
			print("Variant preprocessing failed. Error in the variant calling tools / annotation !!")
			sys.exit(0)

	def variantpreprocessing(self):
		bam=self.primaryfilteredpathogensam_.replace(".sam",".bam")
		#sam format conversion
		cmd=GlobalVar.picard_+" SamFormatConverter INPUT="+self.primaryfilteredpathogensam_+" OUTPUT="+bam+" VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		#raises subprocess.CalledProcessError if returncode is non-zero
		cprocess.check_returncode()
		#Fixmate
		fxdbam=bam.replace(".bam","_fxd.bam")
		cmd=GlobalVar.picard_+" FixMateInformation INPUT="+ bam +" OUTPUT="+fxdbam+" VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#sort sam
		sortedbam=bam.replace(".bam","_fxd_sorted.bam")
		cmd=GlobalVar.picard_+" SortSam INPUT="+ fxdbam +" OUTPUT="+sortedbam+" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#duplicate marking
		rmdupbam=bam.replace(".bam","_fxd_sorted_duprm.bam")
		cmd=GlobalVar.picard_+" MarkDuplicates INPUT="+ sortedbam +" OUTPUT="+rmdupbam+" METRICS_FILE="+rmdupbam.replace(".bam","_info.txt")+ " REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		#bam indexing
		cmd=GlobalVar.picard_+" BuildBamIndex INPUT="+ rmdupbam +" OUTPUT="+rmdupbam+".bai VALIDATION_STRINGENCY=SILENT  TMP_DIR="+self.inputmap_['outdir']+"tmp"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		self.preprocessedpathogenbam_=rmdupbam

	def variantcalling(self):
		#cmd=GlobalVar.strelka_+" --bam "+self.preprocessedpathogenbam_+" --referenceFasta "+GlobalVar.pathofa_+" --runDir "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_strelka"
		cmd=GlobalVar.freebayes_+" -f "+GlobalVar.pathofa_+"  "+self.preprocessedpathogenbam_+" | "+ GlobalVar.vcffilter_+" -f \"QUAL > 20\" > "  +self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_freebayes.vcf"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.lofreq_+" call -f "+GlobalVar.pathofa_+" -o "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_lofreq.vcf "+self.preprocessedpathogenbam_
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.samtools_+" mpileup -q 1 -f "+GlobalVar.pathofa_+" -o "+ self.preprocessedpathogenbam_.replace(".bam",".mpileup")+" "+self.preprocessedpathogenbam_
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.varscan_+" mpileup2cns "+self.preprocessedpathogenbam_.replace(".bam",".mpileup")+" --output-vcf 1 --variants > "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_varscan.vcf"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.bgzip_+" "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_freebayes.vcf"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.bgzip_+" "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_lofreq.vcf"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.bgzip_+" "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_varscan.vcf"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		currenttime=datetime.datetime.now()
		cmd=GlobalVar.bcftools_+" index "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_lofreq.vcf.gz"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.bcftools_+" index "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_varscan.vcf.gz"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.bcftools_+" index "+ self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_freebayes.vcf.gz"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()
		cmd=GlobalVar.bcftools_+" isec -n +2 -o "+self.ipdfinalvcf_+" "+ self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_lofreq.vcf.gz "+ self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_varscan.vcf.gz "+self.inputmap_['outdir']+self.inputmap_['prefix']+"_"+self.samplebasename_+"_freebayes.vcf.gz -f PASS -O v -w 2 "
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()

		#Commented only to be run on Param / otherwise uncomment this and RUN
		cmd=GlobalVar.snpeff_+" -dataDir "+GlobalVar.snpeffdatadir_+" -nodownload ipd1060 "+self.ipdfinalvcf_+" > "+self.ipdfinalannotatedvcf_+" -noStats"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		cprocess.check_returncode()

	def counting_process(self):
		self.secondaryalignment_filtration()
		self.primarysamfilehandle_.close()
		self.assemblyfastqinput1filehandle_.close()
		self.featurecountsoutputfile_=self.primarysam_.replace(".sam", "_featurecounts.tsv")
		self.ipdfinalcountsfile_=self.primarysam_.replace(".sam", "_finalcounts.tsv")
		#if self.data_ispaired_:
		#	self.assemblyfastqinput2filehandle_.close()
		#counting the genes / pathogens
		self.counting()
		#perform assembly
		self.assembly()
		print("Assembly completed..")
		#Alignment of assembled sequences / fasta
		#self.tertiary_alignment()
		#print("Tertiary alignment completed..")

	def filtration(self):
		cmd=""
		if self.data_ispaired_:
			cmd=GlobalVar.fastp_+" -i "+ self.r1fnameraw_ +" -I "+ self.r2fnameraw_ +" -o "+ self.inputmap_['r1'] +" -O "+ self.inputmap_['r2'] +" -h "+ self.fastphtml_ +" -j "+self.fastpjson_+ " -e 20 -l 40 -y"
		else:
			cmd=GlobalVar.fastp_+" -i "+ self.r1fnameraw_ +" -o "+ self.inputmap_['r1'] +" -h "+ self.fastphtml_ +" -j "+self.fastpjson_+" -e 20 -l 40 -y"
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		#raises subprocess.CalledProcessError if returncode is non-zero
		try:
			cprocess.check_returncode()
		except subprocess.CalledProcessError:
			print("Filtration failed !!")
			sys.exit(0)
	#argument: complete output path for sam file
	def primary_alignment(self):
		cmd=""
		if self.data_ispaired_:
			cmd=GlobalVar.hisat2_+" -x "+GlobalVar.hisat2hspathoindex_+" -1 "+self.inputmap_["r1"]+" -2 "+self.inputmap_["r2"]+" --threads "+str(self.inputmap_["threads"])+  " --summary-file "+ self.primaryalignmentsummary_ +" -S " +self.primarysam_
			#cmd=GlobalVar.hisat2_+" -x "+GlobalVar.hisat2pathoindex_+" -1 "+self.inputmap_["r1"]+" -2 "+self.inputmap_["r2"]+" --threads "+str(self.inputmap_["threads"])+  " --summary-file "+ self.primaryalignmentsummary_ +" -S " +self.primarysam_
		else:
			#cmd=GlobalVar.hisat2_+" -x "+GlobalVar.hisat2hspathoindex_+" -U "+self.inputmap_["r1"]+" --threads "+str(self.inputmap_["threads"])+ " --summary-file "+ self.primaryalignmentsummary_ + " -S " +self.primarysam_
			cmd=GlobalVar.hisat2_+" -x "+GlobalVar.hisat2pathoindex_+" -U "+self.inputmap_["r1"]+" --threads "+str(self.inputmap_["threads"])+ " --summary-file "+ self.primaryalignmentsummary_ + " -S " +self.primarysam_
		#print(cmd)
		cprocess=subprocess.run(cmd, shell=True)
		#raises subprocess.CalledProcessError if returncode is non-zero
		try:
			cprocess.check_returncode()
		except Exception:
			print("Primary alignment failed (cmd): %s" %(cmd))
			sys.exit(0)
	'''
	postprocessing performs multiple jobs:
	if human aligned: write to the countsam
	if unaligned: write into fastq -> assembly
	if pathoaligned: perform 1-4 check; if pass a) put into countsam b) write into fastq c) compute genome coverage
	read quality based filtration
	alignment based filtration
	complexity filtration
	secondary alignment and ambiguous read filtration
	physical genome coverage calculation (pathogens)
	create input for megahit (fastq files from pathogen and unaligned reads)
	'''
	def postprimaryalignment_filtration(self):
		primaryfilteredpathogensamfilehandle=pysam.AlignmentFile(self.primaryfilteredpathogensam_, 'w', template=self.primarysamfilehandle_)
		pathoalignedfafilehandle=open(self.pathoalignedfasta_,'w')
		if self.data_ispaired_:
			#assembler1fo=open(self.assemblyfastqinput1_,'w')
			#assembler2fo=open(self.assemblyfastqinput2_,'w')
			for s1 in self.primarysamfilehandle_.fetch():
				IPDFilterStats.totalreads_+=1
				if s1.is_paired:
					try:
						iter=self.primarysamfilehandle_.fetch()
						s2=next(iter)
						if s1.is_unmapped and s2.is_unmapped: # unaligned fragment - put it for assembler input
							self.assemblyfastqinput1filehandle_.write( "@%s\n" %(s1.query_name) )
							self.assemblyfastqinput1filehandle_.write( "%s\n" %(s1.query_sequence) )
							self.assemblyfastqinput1filehandle_.write("+\n")
							self.assemblyfastqinput1filehandle_.write("%s\n" %(s1.qqual) )
							self.assemblyfastqinput2filehandle_.write( "@%s\n" %(s2.query_name) )
							self.assemblyfastqinput2filehandle_.write( "%s\n" %(s2.query_sequence) )
							self.assemblyfastqinput2filehandle_.write("+\n")
							self.assemblyfastqinput2filehandle_.write("%s\n" %(s2.qqual) )
						#else:
						if not s1.is_unmapped or s2.is_unmapped:
							if self.primarymetadata_.is_human(s1.reference_name) and self.primarymetadata_.is_human(s2.reference_name): #Human fragment
								self.primarycountsamfilehandle_.write(s1)
								self.primarycountsamfilehandle_.write(s2)
								#assembly input write
							else:
								if self.primarymetadata_.is_pathogen(s1.reference_name) and self.primarymetadata_.is_pathogen(s2.reference_name): #Pathogen fragment
									read1assess=ShortReadAssessment(s1)
									read2assess=ShortReadAssessment(s2)
									if read1assess.is_low_length() and read2assess.is_low_length():
										IPDFilterStats.lengthfiltered_+=2
									else:
										if read1assess.is_low_phread() and read2assess.is_low_phread():
											IPDFilterStats.qualityfiltered_+=2
										else:
											self.assemblyfastqinput1filehandle_.write( "@%s\n" %(s1.query_name) )
											self.assemblyfastqinput1filehandle_.write( "%s\n" %(s1.query_sequence) )
											self.assemblyfastqinput1filehandle_.write("+\n")
											self.assemblyfastqinput1filehandle_.write("%s\n" %(s1.qqual) )
											self.assemblyfastqinput2filehandle_.write( "@%s\n" %(s2.query_name) )
											self.assemblyfastqinput2filehandle_.write( "%s\n" %(s2.query_sequence) )
											self.assemblyfastqinput2filehandle_.write("+\n")
											self.assemblyfastqinput2filehandle_.write("%s\n" %(s2.qqual) )
											if read1assess.is_low_complexity() and read2assess.is_low_complexity():
												IPDFilterStats.complexityfiltered_+=2
											else:
												#putative pathogen read -> to be sent to secondary alignment and assembly
												#pathofasta for secondary alignment
												pseq1=Sequence(s1.query_name+"___1", s1.query_sequence) #__1 for read 1
												pseq1.write(pathoalignedfafilehandle)
												pseq2=Sequence(s2.query_name+"___2", s2.query_sequence) #__2 for read 2
												pseq2.write(pathoalignedfafilehandle)
												#pathofilteredsam write
												primaryfilteredpathogensamfilehandle.write(s1)
												primaryfilteredpathogensamfilehandle.write(s2)
												#store the samid -> reference name pair
												#self.pathogenaligned_r1_samid_ref_[s1.query_name]=s1.reference_name
												#self.pathogenaligned_r1_samid_ref_[s2.query_name]=s2.reference_name
												self.pathogenaligned_r1_samid_ref_[pseq1.name_]=s1.reference_name
												self.pathogenaligned_r1_samid_ref_[pseq2.name_]=s2.reference_name
					except StopIteration:
						break
		#	assembler1fo.close()
		#	assembler2fo.close()
		else: #For single end data
			#assembly input writer
			#assembler1fo=open(self.assemblyfastqinput1_,'w')
			for s in self.primarysamfilehandle_.fetch():
				IPDFilterStats.totalreads_+=1
				if s.is_unmapped:
					self.assemblyfastqinput1filehandle_.write( "@%s\n" %(s.query_name) )
					self.assemblyfastqinput1filehandle_.write( "%s\n" %(s.query_sequence) )
					self.assemblyfastqinput1filehandle_.write("+\n")
					self.assemblyfastqinput1filehandle_.write("%s\n" %(s.qqual) )
				#else:
				if not s.is_unmapped:
					if self.primarymetadata_.is_human(s.reference_name):
						self.primarycountsamfilehandle_.write(s)
					elif self.primarymetadata_.is_pathogen(s.reference_name):
						#assessment of read in stages: 1 length, 2 phread, 3 complexity, 4 alignment
						readassess=ShortReadAssessment(s)
						if readassess.is_low_length():
							IPDFilterStats.lengthfiltered_+=1
						else:
							if readassess.is_low_phread():
								IPDFilterStats.qualityfiltered_+=1
							else:
								if readassess.is_low_complexity():
									IPDFilterStats.complexityfiltered_+=1
								else:
									#assembly input writing
									self.assemblyfastqinput1filehandle_.write( "@%s\n" %(s.query_name) )
									self.assemblyfastqinput1filehandle_.write( "%s\n" %(s.query_sequence) )
									self.assemblyfastqinput1filehandle_.write("+\n")
									self.assemblyfastqinput1filehandle_.write("%s\n" %(s.qqual) )
									if readassess.is_low_confidence_alignment():
										IPDFilterStats.badalignmentfiltered_+=1
									else:
										#putative pathogen read -> to be sent to secondary alignment and assembly
										#pathofasta for secondary alignment
										pseq=Sequence(s.query_name, s.query_sequence)
										pseq.write(pathoalignedfafilehandle)
										#store the samid->reference name in map
										#self.pathogenaligned_r1_samid_ref_[s.query_name]=s.reference_name
										self.pathogenaligned_r1_samid_ref_[pseq.name_]=s.reference_name
										#pathosamfiltered write
										primaryfilteredpathogensamfilehandle.write(s)
			#assembler1fo.close()
		pathoalignedfafilehandle.close()
		primaryfilteredpathogensamfilehandle.close()
	def secondaryalignment_filtration(self):
		blastouthandle=LocalBlastN(GlobalVar.blastn_, GlobalVar.blastndb_, self.pathoalignedfasta_, 3, self.inputmap_["threads"], self.secondaryalignedtsv_, "6" )
		primaryfilteredpathogensamfilehandle=pysam.AlignmentFile(self.primaryfilteredpathogensam_, 'r')
		#GO through the fasta aligned file and make a map of sam record id -> ambiguous true / false
		blrecord=blastouthandle.get_blastoutputformat6record()
		if self.data_ispaired_: #for paired end data
			bestquerysubjectpair=(blrecord.qseqid_ , blrecord) #initialize this with the first record of alignment file
			while(True):
				try:
					blrecord=blastouthandle.get_blastoutputformat6record()
					#If a _next hit is found for the same sam record
					if blrecord.qseqid_ == bestquerysubjectpair[0]:
						if blrecord.evalue_ < 0.0001: #check if alignment is significant
							#subject ids differ, length greater than or equal to earlier, percent identity is greater than or equal to earlier
							if blrecord.sseqid_ != bestquerysubjectpair[1].sseqid_ and blrecord.length_ > bestquerysubjectpair[1].length_ and blrecord.pident_ >  bestquerysubjectpair[1].pident_:
								bestquerysubjectpair=(blrecord.qseqid_, blrecord)
					else: #Next sam record alignment has started
						#check if subject id is different genus from the sam; if yes - increament the ambigousfiltered_ couter, else write into samcount file
						if not self.secondarymetadata_.is_genus_equal(bestquerysubjectpair[1].sseqid_, self.pathogenaligned_r1_samid_ref_[bestquerysubjectpair[1].qseqid_]):
							self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=False
						else: #If genus is confirmed as same
							self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=True
						bestquerysubjectpair=(blrecord.qseqid_, blrecord)
				except ValueError as e: #End of file for BlastOutputFormat6Record file
					if self.secondarymetadata_.is_genus_equal(bestquerysubjectpair[1].sseqid_, self.pathogenaligned_r1_samid_ref_[bestquerysubjectpair[1].qseqid_]):
						#If genus is confirmed as same
						self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=True
					else:
						self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=True
					break
			#if read is found to be un ambiguous write in samcount file
			#everything is stored in samidr1 : now while retrieving take out for the pair from samfilereader
			for pr1 in primaryfilteredpathogensamfilehandle.fetch():
				iter=primaryfilteredpathogensamfilehandle.fetch()
				try:
					pr2=next(iter)
					#genus is equal in primary alignment
					if self.secondarymetadata_.is_genus_equal(pr1.reference_name, pr2.reference_name): #If according to SAM genus is equal
						#genus is equal in secondary alignment / both should be non-ambiguous
						try:
							if self.samidr1_nonambiguity_[pr1.query_name+"___1"] and self.samidr1_nonambiguity_[pr2.query_name+"___2"] : #is non ambiguous
								self.primarycountsamfilehandle_.write(pr1)
								self.primarycountsamfilehandle_.write(pr2)
							else:
								#Add these reads to assembly input
								IPDFilterStats.ambigousfiltered_+=2
						except KeyError:
							pass
					else:
						#Add these reads to assembly input
						#You may want to handle chimera in later versions
						IPDFilterStats.ambigousfiltered_+=2
				except Exception as e:
					print(e)
					pass
		else: #for single end data
			bestquerysubjectpair=(blrecord.qseqid_ , blrecord) #initialize this witht the first record of alignment file
			while(True):
				try:
					blrecord=blastouthandle.get_blastoutputformat6record()
					#If a _next hit is found for the same sam record
					if blrecord.qseqid_ == bestquerysubjectpair[0]:
						if blrecord.evalue_ < 0.0001: #check if alignment is significant
							#subject ids differ, length greater than or equal to earlier, percent identity is greater than or equal to earlier
							if blrecord.sseqid_ != bestquerysubjectpair[1].sseqid_ and blrecord.length_ > bestquerysubjectpair[1].length_ and blrecord.pident_ >  bestquerysubjectpair[1].pident_:
								bestquerysubjectpair=(blrecord.qseqid_, blrecord)
					else: #Next sam record alignment has started
						#check if subject id is different genus from the sam; if yes - increament the ambigousfiltered_ couter, else write into samcount file
						if not self.secondarymetadata_.is_genus_equal(bestquerysubjectpair[1].sseqid_, self.pathogenaligned_r1_samid_ref_[bestquerysubjectpair[1].qseqid_]):
							self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=False
						else: #If genus is confirmed as same
							self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=True
						bestquerysubjectpair=(blrecord.qseqid_, blrecord)
				except ValueError as e: #End of file for BlastOutputFormat6Record file
					if self.secondarymetadata_.is_genus_equal(bestquerysubjectpair[1].sseqid_, self.pathogenaligned_r1_samid_ref_[bestquerysubjectpair[1].qseqid_]):
						#If genus is confirmed as same
						self.samidr1_nonambiguity_[bestquerysubjectpair[0]]=True
					break
			#if read is found to be un ambiguous write in samcount file
			for pr in primaryfilteredpathogensamfilehandle.fetch():
				try:
					if self.samidr1_nonambiguity_[pr.query_name]: #is non ambiguous
						self.primarycountsamfilehandle_.write(pr)
					else:
						IPDFilterStats.ambigousfiltered_+=1
				except KeyError:
					pass
		primaryfilteredpathogensamfilehandle.close()
		self.primarycountsamfilehandle_.close()
	def counting(self):
		#construct cmd for featurecounts and execute using run
		cmd=""
		if self.data_ispaired_:
			cmd=GlobalVar.featurecounts_+" -a "+GlobalVar.hspathogff_+" -t gene -g gene_id -p -M -O -T " + str(self.inputmap_['threads']) +" --fraction -o "+self.featurecountsoutputfile_+"  "+self.primarycountsam_
		else:
			cmd=GlobalVar.featurecounts_+" -a "+GlobalVar.hspathogff_+" -t gene -g gene_id -M -O -T " + str(self.inputmap_['threads']) +" --fraction -o "+self.featurecountsoutputfile_+"  "+self.primarycountsam_
		cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
		#raises subprocess.CalledProcessError if returncode is non-zero
		cprocess.check_returncode()
		#Need to postprocess the counts created
		featcountfi=open(self.featurecountsoutputfile_)
		ipdfinalcountfo=open(self.ipdfinalcountsfile_,'w')
		ipdfinalcountfo.write("%s\t%s\t%s\t%s\t%s\t%s\n" %("Feature","ID","Length","Fragments", "FPM", "FPKM") )
		for f in featcountfi:
			if not f.startswith("#") and not f.startswith("Geneid"):
				frec=FeatureCountsRecord(f)
				ipdfinalcountfo.write("%s\t%s\t%d\t%f\t%f\t%f\n" %(frec.get_feature(), frec.get_id(), frec.get_length(), frec.get_count(), frec.get_fpm(), frec.get_fpkm() ) )
		featcountfi.close()
		ipdfinalcountfo.close()
	def assembly(self):
		#construct cmd for megahit and execute using fastq 1/2 from
		cmd=""
		if self.data_ispaired_:
			cmd=GlobalVar.megahit_+" -t "+ str(self.inputmap_['threads'])+" -1 "+ self.assemblyfastqinput1_+ " -2 "+ self.assemblyfastqinput2_ +" -o "+self.assembledcontigdirectory_
		else:
			cmd=GlobalVar.megahit_+" -t "+ str(self.inputmap_['threads'])+" -r "+ self.assemblyfastqinput1_+" -o "+self.assembledcontigdirectory_
		cprocess=subprocess.run(cmd, shell=True)
		cprocess.check_returncode()
	def tertiary_alignment(self):
		#perform blast alignment of the assembled
		#subject (reference pathogen) -> contigs mapping to it
		subject_blastout_map={}
		#query subject map; contig ID -> subject
		contig_subject_map={}
		blastouthandle=LocalBlastN(GlobalVar.blastn_, GlobalVar.blastndb_, self.assembledcontigsfinalfile_, 1, self.inputmap_["threads"], self.assembledcontigsblastoutput_, "6" )
		while(True):
			try:
				brec=blastouthandle.get_blastoutputformat6record()
				contig_subject_map[brec.qseqid_]=brec.sseqid_
				if brec.sseqid_ in subject_blastout_map:
					bout=subject_blastout_map[brec.sseqid_]
					bout.append(brec)
					subject_blastout_map[brec.sseqid_]=bout
				else:
					subject_blastout_map[brec.sseqid_]=[brec]
			except ValueError as v:
				break
				print("Breaking out of the blast output loop!!")
		blastouthandle.close_connection()
'''
		ofile=open(self.assembledcontigsblastoutput_.replace(".tsv","_final.tsv"), 'w')
		for s in subject_blastout_map:
			oblast=subject_blastout_map[s]
			if not oblast == None:
				ofile.write("%s\n" %(s) )
			else:
			#write into the file
			#pass
		ofile.close()
'''

'''
prefix
datatype
r1
r2
outdir
threads
'''
