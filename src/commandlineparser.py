#!/usr/bin/python3

import sys
import os.path
import re
import glob

class CustomError(Exception):
	def __int__(self,value):
		self.value=value
	def __str__(self):
		return repr(self.value)

class CommandLineInterfaceParser(object):
	def __init__(self, parameterlist):
		if len(parameterlist) == 8:
			self.runmode=parameterlist[0]
			self.outdir = parameterlist[1]
			self.projname = parameterlist[2]
			self.molecular_type = parameterlist[3]
			self.thread = parameterlist[4]		
			self.sampleinfo=parameterlist[5]
			self.fastq_1=parameterlist[6]
			self.fastq_2=parameterlist[7]
			
	def map_paired_sample(self):
		imap_array=[]
		if self.runmode == "PairedEndSample":
			imap={}
			imap['prefix']=self.projname
			imap['datatype']=self.molecular_type
			imap['r1']=self.fastq_1
			imap['r2']=self.fastq_2
			imap['outdir']=self.outdir
			imap['threads']=self.thread
			imap_array.append(imap)
		else:
			print("Incorrect runmode")
		return imap_array
	
	def map_multiple_sample(self):
		imap_array=[]
		if self.runmode == "Multisample":
			fi = open(self.sampleinfo,'r')
			if all(line.count('\t') == 2 for line in fi):
				fi.close()
				fi = open(self.sampleinfo,'r')
				for sample in fi:
					imap={}
					imap['prefix']=self.projname+"_"+sample.split('\t')[0]
					imap['datatype']=self.molecular_type
					imap['outdir']=self.outdir
					imap['threads']=self.thread
					imap['r1']=sample.split('\t')[1]
					imap['r2']=sample.split('\t')[2].rstrip()
					imap_array.append(imap)
			elif all(line.count('\t') == 1 for line in fi):
				fi.close()
				fi = open(self.sampleinfo,'r')
				for sample in fi:
					imap={}
					imap['prefix']=self.projname+"_"+sample.split('\t')[0]
					imap['datatype']=self.molecular_type
					imap['outdir']=self.outdir
					imap['threads']=self.thread
					imap['r1']=sample.split('\t')[1].rstrip()
					imap_array.append(imap)
			else:
				print("Invalid File Format")
		else:
			print("Incorrect runmode")
		return imap_array	
	
	def map_single_end_sample(self):
		imap_array=[]
		if self.runmode == "SingleEndSample":
			imap={}
			imap['prefix']=self.projname
			imap['datatype']=self.molecular_type
			imap['r1']=self.fastq_1
			imap['outdir']=self.outdir
			imap['threads']=self.thread
			imap_array.append(imap)
		else:
			print("Incorrect runmode")
		return imap_array
	
