#!/usr/bin/python3
'''

Author	  : Sanket Desai
Date		: 25/04/2020
Description : Static classes with global variables storing tool and database paths

'''
import sys
import os

class GlobalFormatException(Exception):
	def __init__(self, exstr):
		Exception.__init__(self, exstr)
class GlobalVarException(Exception):
	def __init__(self, exstr=""):
		Exception.__init__(self, exstr)

class InputFileTypeException(Exception):
	def __init__(self, message):
		super().__init__(message)

class GlobalVar(object):
	fastp_="../external/fastp-0.20.1/fastp"
	hisat2_="../external/hisat2-2.1.0/hisat2"
	hisat2_build_="../external/hisat2-2.1.0/hisat2-build"
	blastn_="../external/ncbi-blast-2.10.0+/bin/blastn"
	makeblastdb_="../external/ncbi-blast-2.10.0+/bin/makeblastdb"
	featurecounts_="../external/subread-2.0.0-Linux-x86_64/bin/featureCounts"
	samtools_="../external/samtools-1.10/samtools"
	megahit_="../external/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit"
	picard_="java -Xmx10G -jar ../external/picard-2.23/picard.jar"
	bcftools_="../external/bcftools/bcftools"
	varscan_="java -Xmx10G -jar ../external/varscan/VarScan.v2.4.4.jar"
	bgzip_="../external/tabix/bgzip"
	tabix_="../external/tabix/tabix"
	lofreq_="../external/lofreq_star-2.1.2/bin/lofreq"
	snpeff_="java -Xmx10G -jar ../external/snpEff/snpEff.jar"
	freebayes_="../external/freebayes/bin/freebayes"
	vcffilter_="../external/freebayes/vcflib/build/vcffilter"

	##snpEff variant annotation
	snpeffdatadir_="./data"
	snpeffgbk_="../external/snpEff/data/ipd1060/genes.gbk"

	pathofa_="../data/primaryref/pathoref/patho.fa"
	blastndb_="../data/secondaryref/secondary.fa" #give human + pathogen db name
	primaryannotationtsv_="../data/annotation/primaryannotation.tsv"
	secondaryannotationtsv_="../data/annotation/secondaryannotation.tsv"
	hisat2hspathoindex_="../data/primaryref/hspathoref/hspatho"
	hisat2pathoindex_="../data/primaryref/pathoref/patho"
	hspathogff_="../data/annotation/hspatho.gff"

	pathoids_="../data/annotation/patho.ids"
	secondarydbids_="../data/annotation/secondarydb.ids"
	##Long read analysis
	porechop_="python3 ../external/Porechop/porechop-runner.py"
	minimap2_="../external/minimap2/minimap2"

	minimap2index_="../data/primaryref/pathoref/patho.minimap2.ref"
	pathogff_="../data/annotation/patho.gff"
	#No entry for NanoFilt and longshot, since these will be directly run on CLI

	##cov2 module
	hspathofa_="../data/primaryref/hspathoref/hspatho.fa"
	vcfdatabasefile_="../data/cov2moduleref/gisaidmerged.vcf.gz"
	gisaidmetadatafile_="../data/cov2moduleref/gisaidmetadata.tsv"
	covblastndb_="../data/cov2moduleref/gisaid.fasta"

	@staticmethod
	def initialize():
		#tools
		########if not os.path.isfile(GlobalVar.hisat2_):
		#	raise GlobalVarException("Hisat2 not found in the defined path: %s" % GlobalVar.hisat2_ )
		if not os.path.isfile(GlobalVar.blastn_):
			raise GlobalVarException("BlastN not found in the defined path: %s" % GlobalVar.blastn_ )
		if not os.path.isfile(GlobalVar.samtools_):
			raise GlobalVarException("Samtools not found in the defined path: %s" % GlobalVar.samtools_ )
		#if not os.path.isfile(GlobalVar.megahit_):
		#	raise GlobalVarException("MegaHit not found in the defined path: %s" % GlobalVar.megahit_ )
		#if not os.path.isfile(GlobalVar.subjunc_):
		 #   raise GlobalVarException("Subjunc not found in the defined path: %s" % GlobalVar.subjunc_ )
		#if not os.path.isfile(GlobalVar.subreadalign_):
		 #   raise GlobalVarException("Subread-align not found in the defined path: %s" % GlobalVar.subreadaligner_ )
		if not os.path.isfile(GlobalVar.featurecounts_):
			raise GlobalVarException("FeatureCounts not found in the defined path: %s" % GlobalVar.featurecounts_ )
		#databases
		if not os.path.isfile(GlobalVar.primaryannotationtsv_):
			raise GlobalVarException("Annotation file not found in the defined path: %s" % GlobalVar.primaryannotationtsv_ )
		if not os.path.isfile(GlobalVar.secondaryannotationtsv_):
			raise GlobalVarException("Annotation file not found in the defined path: %s" % GlobalVar.secondaryan )
		if not os.path.isfile(GlobalVar.hspathogff_):
			raise GlobalVarException("Human-pathogen GFF file not found in the defined path: %s" % GlobalVar.hspathogff_ )
		#if not os.path.isfile(GlobalVar.hisat2hspathoindex_+".1.ht2l"):
		#	raise GlobalVarException("Human-pathogen HISAT index file not found in the defined path: %s" % GlobalVar.hisat2hspathoindex_+".1.ht2l" )
		#if not os.path.isfile(GlobalVar.blastndb_+".15.nsq"):
		#	raise GlobalVarException("Secondary database index file not found in the defined path: %s" % (GlobalVar.blastndb_) )
		#if not os.path.isfile(GlobalVar.hisat2pathoindex_+".1.ht2"):
		#	raise GlobalVarException("Pathogen HISAT index file not found in the defined path: %s" % GlobalVar.hisat2pathoindex_+".1.ht2" )
		if not os.path.isfile(GlobalVar.pathoids_):
			raise GlobalVarException("Pathoge ID file not found in the defined path : %s" %(GlobalVar.pathoids_) )
		if not os.path.isfile(GlobalVar.secondarydbids_):
			raise GlobalVarException("Secondary database ID file not found in the defined path: %s" %(GlobalVar.secondarydbids_))
		#check all the is files for all the databases and tools

class PrimaryMetadata(object):
	def __init__(self, anntsv):
		fi=open(anntsv, 'r')
		#self.header_=fi.readline().strip().split("\t")
		self.metadata_={}
		self.humanids_=[]
		self.pathogenids_=[]
		for l in fi:
			sl=l.strip().split("\t")
			self.metadata_[sl[0]]=sl[1:7]
		fi.close()
		self.humanids_=self.get_human_ids()
		self.pathogenids_=self.get_pathogen_ids()
		print("Primary data loaded!")
	def get_primary_ids(self):
		return list(self.metadata_.keys())
	def get_human_ids(self):
		hids=[]
		for k in self.metadata_:
			temp=self.metadata_[k]
			if temp[1]=='Homo sapiens':
				hids.append(k)
		return hids
	def get_pathogen_ids(self):
		pids=[]
		for k in self.metadata_:
			temp=self.metadata_[k]
			if temp[1]!='Homo sapiens':
				pids.append(k)
		return pids
	def get_pathogen_metadata(self):
		pathomap={}
		for k in self.metadata_:
			temp=self.metadata_[k]
			if temp[1]!='Homo sapiens':
				pathomap[k]=self.metadata_[k]
		return pathomap
	def get_human_metadata(self):
		hsmap={}
		for k in self.metadata_:
			temp=self.metadata_[k]
			if temp[1]=='Homo sapiens':
				hsmap[k]=self.metadata_[k]
		return hsmap
	def get_pathogenid_length_map(self):
		pmap=self.get_pathogen_metadata()
		pid_len_map={}
		for p in pmap:
			pid_len_map[p]=int(pmap[p][2])
		return pid_len_map
	def get_humanid_length_map(self):
		pmap=self.get_human_metadata()
		pid_len_map={}
		for p in pmap:
			pid_len_map[p]=int(pmap[p][2])
		return pid_len_map
	def get_genus(self, id_):
		return str(self.metadata_[id_][0]).strip()
	def is_genus_equal(self, id1, id2):
		return self.get_genus(id1) == self.get_genus(id2)
	def get_organism(self, id_):
		return str(self.metadata_[id_][1]).strip()
	def get_sequence_length(self, id_):
		return int(str(self.metadata_[id_][2]).strip())
	def get_description(self, id_):
		return str(self.metadata_[id_][4]).strip()
	def get_taxonomy(self, id_):
		return str(self.metadata_[id_][5]).strip()
	def get_plasmid_flag(self, id_):
		return str(self.metadata_[id_][3]).strip()
	def is_human(self, id_):
		return id_ in self.humanids_
	def is_pathogen(self, id_):
		return id_ in self.pathogenids_

class SecondaryMetadata(PrimaryMetadata):
	def __init__(self, anntsv):
		print("Loading secondary data!!")
		super().__init__(anntsv)
		print("secondary data loaded..")
