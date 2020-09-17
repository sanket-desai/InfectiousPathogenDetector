import pysam
import re
import sys
import os
import numpy as np

class CustomError(Exception):
	def __int__(self,value):
		self.value=value
	def __str__(self):
		return repr(self.value)
class GisaidVcfParser(object):
	def __init__(self, VcfFile=None):
		if os.path.exists(VcfFile):
			self.GisaidVCFdatabase = pysam.VariantFile(VcfFile)
			self.key_list=[]
			self.file=list(self.GisaidVCFdatabase.header.samples)
			self.sample_list=[]
			for i in self.file:
				if i.find('_EPI_ISL_') == -1:
					self.sample_list.append(None)
				else:
					self.sample_list.append("EPI_ISL_"+i.split('_EPI_ISL_')[1].split('_')[0])
			self.key_sample_matrix=[]
			for variant in self.GisaidVCFdatabase.fetch():
				i=0
				temp_sample_list = [0] * len(self.sample_list)
				for sample in variant.samples:
					if variant.samples[sample]['DP']:
						temp_sample_list[i] = 1
					i=i+1
				self.key_sample_matrix.append(temp_sample_list)   
				self.key_list.append(str(str(variant.chrom)+"_"+str(variant.pos)+"_"+str(variant.ref)+"_"+str(','.join(str(x) for x in variant.alts))))
			#print(len(self.samples))
			#print(len(self.file))
			#print(len(self.key_sample_matrix))
			#print(len(self.key_sample_matrix[0]))
			
		else:
			sys.exit(0)
	
	def get_genomic_key(self):
		return self.key_list
		
	def get_sample_list(self):
		return self.sample_list
	
	def is_genomic_key_present(self,genomic_key):
		if genomic_key:
			if genomic_key in self.key_list:
				return 1
			else:
				return 0

	def is_gisaid_present(self,gisaid):
		if gisaid:
			if gisaid in self.sample_list:
				return 1
			else:
				return 0

	
	def get_mutational_profile_given_gisaid(self,gisaid):
		mutational_profile=[]
		if gisaid:
			if self.is_gisaid_present(gisaid):
				index=self.sample_list.index(gisaid)
				for i in self.key_sample_matrix:
					mutational_profile.append(i[index])
		return mutational_profile
	
	def get_genomic_key_given_gisaid(self,gisaid):
		key_list=[]
		if gisaid:
			if self.is_gisaid_present(gisaid):
				mutational_profile = self.get_mutational_profile_given_gisaid(gisaid)
				index_pos=list(np.nonzero(mutational_profile)[0])
				for i in index_pos:
					key_list.append(self.key_list[i])
		return key_list

	def get_gisaid_given_genomic_key(self,genomic_key):
		gisaid_key=[]
		if genomic_key:
			if self.is_genomic_key_present:
				index=self.key_list.index(genomic_key)
				gisaid_complete_list=self.key_sample_matrix[index]
				index_gisaid_present=list(np.nonzero(gisaid_complete_list)[0])
				for i in index_gisaid_present:
					if self.sample_list[i]:
						gisaid_key.append(self.sample_list[i])
		return gisaid_key